"""The *NEXAFS edge step* dialog: choose the pre-edge window (and, for the
edge-step method, the post-edge window and the edge energy E0) of a NEXAFS
spectrum and see the result before it is saved.

The numbers come from ``nexafs`` (Tk-free). What the dialog shows first is a
*suggestion* (``nexafs.suggest_edge``: the template's first 4 eV for the
pre-edge, the last tenth of the scan for the post-edge, the steepest rise for
E0); nothing is applied or stated in a report until *Use for this spectrum*
saves it in the annotations (the workbook), as with the ISS beam. Windows are
dragged on the upper plot (choose which with the buttons), E0 is a click.

The dialog talks to the app through ``app.nexafs_regions()``,
``app.nexafs_base(r)`` (the energy and counts before any edge normalisation),
``app.nexafs_edge_get(r)``, ``app.nexafs_edge_set(regions, params)``,
``app.calibration_label(r)`` and ``app.palette``. Not modal. Its behaviour is
in plain methods (``set_mode``, ``set_values``, ``suggest``, ``preview``,
``use_here``, ``use_all``, ``remove``, ``select``) so tests drive it without a
mouse.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk

import nexafs
from workbook_ui import _finish, keep_in_front

ORDERS = ("constant", "linear", "quadratic")


class NexafsEdgeDialog(tk.Toplevel):
    def __init__(self, master, app):
        super().__init__(master)
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure
        self.app = app
        self.regions = app.nexafs_regions()
        self.region = self.regions[0]
        self.saved = False
        self.title("NEXAFS edge step")
        body = ttk.Frame(self, padding=12)
        body.pack(fill="both", expand=True)
        body.columnconfigure(1, weight=1)

        ttk.Label(body, text="Spectrum").grid(row=0, column=0, sticky="w",
                                              padx=(0, 8))
        self.reg_cb = ttk.Combobox(
            body, state="readonly", width=56,
            values=[app.calibration_label(r) for r in self.regions])
        self.reg_cb.grid(row=0, column=1, sticky="ew")
        self.reg_cb.current(0)
        self.reg_cb.bind("<<ComboboxSelected>>",
                         lambda e: self.select(self.reg_cb.current()))

        self.mode = tk.StringVar(value="step")
        modes = ttk.Frame(body)
        modes.grid(row=1, column=0, columnspan=2, sticky="w", pady=(8, 2))
        ttk.Label(modes, text="Method").pack(side="left", padx=(0, 8))
        for m in nexafs.MODES:
            ttk.Radiobutton(modes, text=nexafs.MODE_NAMES[m], value=m,
                            variable=self.mode,
                            command=self._mode_changed).pack(side="left",
                                                             padx=(0, 12))

        grid = ttk.Frame(body)
        grid.grid(row=2, column=0, columnspan=2, sticky="w")
        self.vars = {k: tk.StringVar() for k in
                     ("pre_lo", "pre_hi", "post_lo", "post_hi", "e0")}
        self.orders = {"pre": tk.StringVar(value=ORDERS[1]),
                       "post": tk.StringVar(value=ORDERS[1])}
        self.entries = {}

        def row(r, label, lo, hi, order):
            ttk.Label(grid, text=label).grid(row=r, column=0, sticky="w",
                                             pady=2)
            for c, key in ((1, lo), (3, hi)):
                e = ttk.Entry(grid, textvariable=self.vars[key], width=9)
                e.grid(row=r, column=c, padx=(6, 2))
                e.bind("<Return>", lambda ev: self._redraw())
                e.bind("<FocusOut>", lambda ev: self._redraw())
                self.entries[key] = e
            ttk.Label(grid, text="to").grid(row=r, column=2)
            ttk.Label(grid, text="eV   fit").grid(row=r, column=4, padx=(4, 4))
            cb = ttk.Combobox(grid, state="readonly", width=10, values=ORDERS,
                              textvariable=self.orders[order])
            cb.grid(row=r, column=5)
            cb.bind("<<ComboboxSelected>>", lambda ev: self._redraw())
            self.entries[order] = cb

        row(0, "Pre-edge", "pre_lo", "pre_hi", "pre")
        row(1, "Post-edge", "post_lo", "post_hi", "post")
        ttk.Label(grid, text="Edge energy E0").grid(row=2, column=0,
                                                    sticky="w", pady=2)
        e = ttk.Entry(grid, textvariable=self.vars["e0"], width=9)
        e.grid(row=2, column=1, padx=(6, 2))
        e.bind("<Return>", lambda ev: self._redraw())
        e.bind("<FocusOut>", lambda ev: self._redraw())
        self.entries["e0"] = e
        ttk.Label(grid, text="eV").grid(row=2, column=2, sticky="w")

        self.drag = tk.StringVar(value="pre")
        bar = ttk.Frame(body)
        bar.grid(row=3, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Label(bar, text="Set on the plot:").pack(side="left", padx=(0, 6))
        self.drag_buttons = {}
        for key, text in (("pre", "pre-edge (drag)"),
                          ("post", "post-edge (drag)"), ("e0", "E0 (click)")):
            b = ttk.Radiobutton(bar, text=text, value=key, variable=self.drag,
                                command=self._drag_changed)
            b.pack(side="left", padx=(0, 10))
            self.drag_buttons[key] = b

        self.fig = Figure(figsize=(8.2, 4.6), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=body)
        self.canvas.get_tk_widget().grid(row=4, column=0, columnspan=2,
                                         sticky="nsew", pady=(6, 4))
        body.rowconfigure(4, weight=1)
        self.status = ttk.Label(body, style="Muted.TLabel", wraplength=700,
                                justify="left")
        self.status.grid(row=5, column=0, columnspan=2, sticky="w")

        btns = ttk.Frame(body)
        btns.grid(row=6, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        ttk.Button(btns, text="Suggest", command=self.suggest).pack(
            side="left")
        self.use_btn = ttk.Button(btns, text="Use for this spectrum",
                                  command=self.use_here)
        self.use_btn.pack(side="left", padx=6)
        self.all_btn = ttk.Button(btns, command=self.use_all)
        self.all_btn.pack(side="left")
        ttk.Button(btns, text="Remove from this spectrum",
                   command=self.remove).pack(side="left", padx=6)
        ttk.Button(btns, text="Close", command=self.destroy).pack(
            side="right")
        self.protocol("WM_DELETE_WINDOW", self.destroy)
        self.bind("<Escape>", lambda e: self.destroy())
        _finish(self, app, 920, 760)
        keep_in_front(self, master)
        self.span = None
        self._cid = None
        self.select(0)

    # -- what the dialog edits ---------------------------------------------------
    def select(self, index):
        """Work on spectrum ``index`` of the list: its saved settings, else a
        suggestion (not saved)."""
        self.region = self.regions[index]
        self.reg_cb.current(index)
        energy, counts = self.app.nexafs_base(self.region)
        self._energy, self._counts = energy, counts
        saved = self.app.nexafs_edge_get(self.region)
        self.saved = bool(saved)
        self.all_btn.configure(
            text=f"Use for all {len(self.regions)} listed"
            if len(self.regions) > 1 else "Use for all listed")
        self.all_btn.state(["!disabled"] if len(self.regions) > 1
                           else ["disabled"])
        if saved:
            self._fill(saved)
        else:
            self._fill(nexafs.suggest_edge(energy, counts, self.mode.get())
                       or {"mode": self.mode.get(), "pre": [None, None],
                           "post": [None, None], "e0": None})
        self._redraw()

    def _fill(self, p):
        if p.get("mode") in nexafs.MODES:
            self.mode.set(p["mode"])
        pre, post = p.get("pre") or [None, None], p.get("post") or [None, None]
        for key, v in (("pre_lo", pre[0]), ("pre_hi", pre[1]),
                       ("post_lo", post[0]), ("post_hi", post[1]),
                       ("e0", p.get("e0"))):
            self.vars[key].set("" if v is None else f"{v:.2f}".rstrip("0")
                               .rstrip("."))
        self.orders["pre"].set(ORDERS[nexafs._order(p.get("pre_order"))])
        self.orders["post"].set(ORDERS[nexafs._order(p.get("post_order"))])

    def set_mode(self, mode):
        self.mode.set(mode)
        self._mode_changed()

    def set_values(self, pre=None, post=None, e0=None, pre_order=None,
                   post_order=None):
        """Type values in (windows as ``(lo, hi)``), as the entries would."""
        if pre is not None:
            self.vars["pre_lo"].set(f"{pre[0]:g}")
            self.vars["pre_hi"].set(f"{pre[1]:g}")
        if post is not None:
            self.vars["post_lo"].set(f"{post[0]:g}")
            self.vars["post_hi"].set(f"{post[1]:g}")
        if e0 is not None:
            self.vars["e0"].set(f"{e0:g}")
        if pre_order is not None:
            self.orders["pre"].set(ORDERS[pre_order])
        if post_order is not None:
            self.orders["post"].set(ORDERS[post_order])
        self._redraw()

    def params(self):
        """The settings as typed (a valid dict, or {} while incomplete)."""
        def f(key):
            try:
                return float(self.vars[key].get())
            except ValueError:
                return None
        p = {"mode": self.mode.get(), "pre": [f("pre_lo"), f("pre_hi")],
             "pre_order": ORDERS.index(self.orders["pre"].get())}
        if p["mode"] == "step":
            p.update(post=[f("post_lo"), f("post_hi")], e0=f("e0"),
                     post_order=ORDERS.index(self.orders["post"].get()))
        return nexafs.sanitise_edge(p)

    def suggest(self):
        s = nexafs.suggest_edge(self._energy, self._counts, self.mode.get())
        if not s:
            self.status.configure(text="Too few points to suggest windows.")
            return
        self._fill(s)
        self._redraw()

    def preview(self):
        """``(normalised counts, info)`` for what is typed, or ``(None,
        reasons)``."""
        p = self.params()
        if not p:
            return None, ["fill in every window and E0"]
        return nexafs.edge_normalise(self._energy, self._counts, p)

    def use_here(self):
        p = self.params()
        if not p or nexafs.check_edge(p, self._energy, self._counts):
            self._redraw()
            return False
        self.app.nexafs_edge_set([self.region], p)
        self._after_save()
        return True

    def use_all(self):
        """The same windows for every listed spectrum that they fit (the
        windows are energies, so one edge measured on several channels or
        samples takes one setting); the others are named."""
        p = self.params()
        if not p:
            self._redraw()
            return []
        fits, skipped = [], []
        for r in self.regions:
            energy, counts = self.app.nexafs_base(r)
            (skipped if nexafs.check_edge(p, energy, counts) else fits
             ).append(r)
        self.app.nexafs_edge_set(fits, p)
        self._after_save()
        if skipped:
            self.status.configure(
                text=f"Saved for {len(fits)} spectra; the windows do not fit "
                     "the others: " + ", ".join(self.app.calibration_label(r)
                                                for r in skipped[:4]))
        return fits

    def remove(self):
        self.app.nexafs_edge_set([self.region], {})
        self._after_save()

    def _after_save(self):
        self.saved = bool(self.app.nexafs_edge_get(self.region))
        self._energy, self._counts = self.app.nexafs_base(self.region)
        self._redraw()

    # -- drawing -----------------------------------------------------------------
    def _mode_changed(self):
        step = self.mode.get() == "step"
        for key in ("post", "e0"):
            self.drag_buttons[key].state(["!disabled"] if step
                                         else ["disabled"])
        if not step and self.drag.get() != "pre":
            self.drag.set("pre")
        self._drag_changed()
        self._redraw()

    def _drag_changed(self):
        if self.span is not None:
            self.span.set_active(self.drag.get() in ("pre", "post"))

    def _redraw(self):
        step = self.mode.get() == "step"
        for key in ("post_lo", "post_hi", "e0"):
            self.entries[key].configure(state="normal" if step else "disabled")
        self.entries["post"].configure(state="readonly" if step
                                       else "disabled")
        pal = self.app.palette
        fig = self.fig
        fig.clear()
        fig.set_facecolor(pal["plot_bg"])
        top = fig.add_subplot(211)
        bot = fig.add_subplot(212, sharex=top)
        for ax in (top, bot):
            ax.set_facecolor(pal["plot_bg"])
            ax.tick_params(colors=pal["plot_fg"], labelsize=8)
            for s in ax.spines.values():
                s.set_color(pal["plot_fg"])
        self.top, self.bot = top, bot
        r = self.region
        top.plot(self._energy, self._counts, color=pal["accent"], lw=1.2)
        top.set_ylabel(f"{r.count_label} ({r.count_units})", fontsize=8,
                       color=pal["plot_fg"])
        p = self.params()
        out, info = self.preview()
        if p:
            for w, key in ((p["pre"], "pre"),
                           (p.get("post"), "post")):
                if w:
                    top.axvspan(w[0], w[1], color=pal["accent"], alpha=0.16,
                                lw=0)
            if "e0" in p:
                top.axvline(p["e0"], color=pal["plot_fg"], lw=0.8, ls=":")
        if out is not None:
            for fit, color in ((info["pre"], "#d55e00"),
                               (info.get("post"), "#009e73")):
                if fit:
                    top.plot(self._energy, [nexafs.poly_at(fit, e)
                                            for e in self._energy], color=color,
                             lw=1.0, ls="--")
            bot.plot(self._energy, out, color=pal["accent"], lw=1.2)
            for level in (0.0, 1.0):
                bot.axhline(level, color=pal["plot_fg"], lw=0.6, ls=":")
        bot.set_xlabel("Photon energy (eV)", fontsize=8, color=pal["plot_fg"])
        bot.set_ylabel("Normalised" + (" (step = 1)" if step
                                       else " (0 to 1)"), fontsize=8,
                       color=pal["plot_fg"])
        fig.tight_layout()
        self._selectors()
        self._status(p, info if out is not None else None,
                     None if out is not None else info)
        self.canvas.draw_idle()

    def _selectors(self):
        from matplotlib.widgets import SpanSelector
        pal = self.app.palette
        self.span = SpanSelector(
            self.top, self._on_span, "horizontal", useblit=False,
            interactive=False, props=dict(facecolor=pal["accent"], alpha=0.25))
        self.span.set_active(self.drag.get() in ("pre", "post"))
        if self._cid is not None:
            try:
                self.canvas.mpl_disconnect(self._cid)
            except Exception:                       # noqa: BLE001
                pass
        self._cid = self.canvas.mpl_connect("button_press_event", self._click)

    def _on_span(self, lo, hi):
        if hi - lo <= 0:
            return
        self.set_window(self.drag.get(), lo, hi)

    def _click(self, event):
        if (self.drag.get() == "e0" and event.inaxes is self.top
                and event.xdata is not None and event.button == 1):
            self.vars["e0"].set(f"{event.xdata:.2f}")
            self._redraw()

    def set_window(self, which, lo, hi):
        """What dragging on the plot does."""
        self.set_values(**{which: (round(lo, 2), round(hi, 2))})

    def _status(self, p, info, problems):
        if info is not None:
            bits = [nexafs.edge_text(p) + "."]
            if "step" in info:
                bits.append(f"Edge step {info['step']:.4g} "
                            f"{self.region.count_units}.")
            bits.append("Saved with this file." if self.saved and
                        self.app.nexafs_edge_get(self.region) == p else
                        "Not saved yet: nothing is applied or stated until "
                        "you choose Use for this spectrum.")
            self.status.configure(text=" ".join(bits))
        else:
            why = problems or ["fill in every window and E0"]
            self.status.configure(text="Cannot apply: " + "; ".join(why)
                                  + ".")

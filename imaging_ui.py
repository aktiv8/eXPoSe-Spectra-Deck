"""The viewer for Kratos stigmatic imaging maps (see ``kratosmap``).

Left: one image at a time, with a µm scale and a scale bar. Right: how the
mean counts of the chosen area (the whole image until a box is dragged on it)
change from frame to frame: against stage height for a focus series, against
time for a repeated measurement, else against the frame number; or how sharp
the image is, which peaks at best focus. Step through the maps of the file
with the slider, the arrows or the list, narrow them to the same element or
the same stage position, change the colours and smooth the picture, and save
the picture, the pixel values or the table of frames.

The pixel size is approximate (the file does not record the field of view; see
``kratosmap``) and the picture says so. Smoothing is for the picture only: the
area means, the sharpness and the saved values come from the counts as
recorded.

The dialog talks to the app through ``app.palette``, ``app.cfg``, ``app.status``
and ``app.themes``. Not modal.
"""

from __future__ import annotations

import csv
import os
import tkinter as tk
from tkinter import filedialog, ttk

import kratosmap
import snapmap
import themes
from workbook_ui import _finish

FILTERS = (("All maps", "all"), ("Same energy", "energy"),
           ("Same stage position", "position"))
PLOTS = ("Mean counts per pixel", "Sharpness (focus)")
BAR_CHOICES = (10, 20, 50, 100, 200, 500)        # µm


def _scale_bar(width_um):
    """A round scale-bar length (µm) of about a fifth of ``width_um``."""
    want = width_um / 5
    return min(BAR_CHOICES, key=lambda b: abs(b - want))


class MapSeriesDialog(tk.Toplevel):
    def __init__(self, master, app, parser, region):
        super().__init__(master)
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure
        self.app, self.parser = app, parser
        self.all = kratosmap.frames(parser.regions)
        start = next((f for f in self.all if f.region is region), self.all[0])
        cfg = app.cfg.get("imaging", {})
        self.scale = tk.StringVar(value=cfg.get("scale", "Viridis")
                                  if cfg.get("scale") in themes.SCALE_NAMES
                                  else "Viridis")
        self.sigma = tk.DoubleVar(value=float(cfg.get("sigma", 0.0)))
        self.common = tk.BooleanVar(value=bool(cfg.get("common", False)))
        self.plot_kind = tk.StringVar(
            value=cfg.get("plot", PLOTS[0]) if cfg.get("plot") in PLOTS
            else PLOTS[0])
        auto = kratosmap.default_filter(self.all, start)
        self.how = tk.StringVar(value=next(n for n, k in FILTERS if k == auto))
        self.mask = None                         # boolean (ny, nx); None = whole
        self._sharp = {}                         # id(region) -> sharpness
        self.shown = self._select(start)
        self.i = self.shown.index(start)
        self.title(f"Maps — {os.path.basename(parser.path or '') or region.sample}")

        body = ttk.Frame(self, padding=8)
        body.pack(fill="both", expand=True)
        bar = ttk.Frame(body)
        bar.pack(fill="x")
        ttk.Label(bar, text="Show").pack(side="left")
        cb = ttk.Combobox(bar, state="readonly", width=19, textvariable=self.how,
                          values=[n for n, _k in FILTERS])
        cb.pack(side="left", padx=(4, 12))
        cb.bind("<<ComboboxSelected>>", lambda e: self.set_filter(self.how.get()))
        ttk.Label(bar, text="Colours").pack(side="left")
        cb = ttk.Combobox(bar, state="readonly", width=12, textvariable=self.scale,
                          values=themes.SCALE_NAMES)
        cb.pack(side="left", padx=(4, 12))
        cb.bind("<<ComboboxSelected>>", lambda e: self._settings_changed())
        ttk.Checkbutton(bar, text="Same colour range for all", variable=self.common,
                        command=self._settings_changed).pack(side="left")
        ttk.Label(bar, text="Smooth (px)").pack(side="left", padx=(12, 0))
        ttk.Scale(bar, from_=0.0, to=4.0, variable=self.sigma, length=90,
                  command=lambda v: self._smooth_changed()).pack(side="left",
                                                                 padx=(4, 0))

        nav = ttk.Frame(body)
        nav.pack(fill="x", pady=(4, 0))
        ttk.Button(nav, text="◀", width=3,
                   command=lambda: self.step(-1)).pack(side="left")
        self.slider = ttk.Scale(nav, from_=0, to=max(0, len(self.shown) - 1),
                                command=self._slid)
        self.slider.pack(side="left", fill="x", expand=True, padx=6)
        ttk.Button(nav, text="▶", width=3,
                   command=lambda: self.step(1)).pack(side="left")
        self.frame_cb = ttk.Combobox(nav, state="readonly", width=44)
        self.frame_cb.pack(side="left", padx=(10, 0))
        self.frame_cb.bind("<<ComboboxSelected>>",
                           lambda e: self.set_frame(self.frame_cb.current()))

        self.fig = Figure(figsize=(10.4, 5.0), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=body)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, pady=(6, 4))
        row = ttk.Frame(body)
        row.pack(fill="x")
        self.hint = ttk.Label(row, style="Muted.TLabel", text="", anchor="w")
        self.hint.pack(side="left", fill="x", expand=True)
        ttk.Combobox(row, state="readonly", width=22, textvariable=self.plot_kind,
                     values=PLOTS).pack(side="right")
        self.plot_kind.trace_add("write", lambda *_a: self._settings_changed())
        ttk.Button(row, text="Whole image", command=self.clear_roi).pack(
            side="right", padx=6)
        ttk.Button(row, text="Save…", command=self._save_menu).pack(side="right")
        self.protocol("WM_DELETE_WINDOW", self.destroy)
        self.bind("<Escape>", lambda e: self.destroy())
        _finish(self, app, 1080, 660)
        self._sync_nav()
        self._rebuild()

    # -- what is on show -------------------------------------------------------------
    def _select(self, current):
        key = dict(FILTERS).get(self.how.get(), "all")
        return kratosmap.select(self.all, current, key) or list(self.all)

    @property
    def frame(self):
        return self.shown[self.i]

    def set_filter(self, label):
        """Narrow the frames to those sharing the current one's energy or
        stage position (or all of them); the current frame stays on show."""
        cur = self.frame
        self.how.set(label)
        self.shown = self._select(cur)
        self.i = self.shown.index(cur) if cur in self.shown else 0
        self._sync_nav()
        self._rebuild()

    def set_frame(self, i):
        i = max(0, min(len(self.shown) - 1, int(i)))
        if i == self.i:
            return
        self.i = i
        if (self.mask is not None and
                self.mask.shape != (self.frame.cube.ny, self.frame.cube.nx)):
            self.mask = None
        self._sync_nav()
        self._rebuild()

    def step(self, d):
        self.set_frame(self.i + d)

    def _slid(self, value):
        self.set_frame(round(float(value)))

    def _sync_nav(self):
        n = len(self.shown)
        self.slider.configure(to=max(1, n - 1) if n > 1 else 1)
        self.slider.state(["!disabled"] if n > 1 else ["disabled"])
        self.slider.set(self.i)
        self.frame_cb.configure(values=[f"{k + 1}. {f.label}"
                                        for k, f in enumerate(self.shown)])
        self.frame_cb.current(self.i)

    def set_roi(self, mask):
        self.mask = mask if mask is not None and mask.any() else None
        self._rebuild()

    def clear_roi(self):
        self.set_roi(None)

    # -- settings --------------------------------------------------------------------
    def _save_cfg(self):
        self.app.cfg["imaging"] = {
            "scale": self.scale.get(), "sigma": round(self.sigma.get(), 2),
            "common": self.common.get(), "plot": self.plot_kind.get()}

    def _settings_changed(self):
        self._save_cfg()
        self._rebuild()

    def _smooth_changed(self):
        self._save_cfg()
        self._rebuild()

    def _cmap(self):
        return themes.scale_colourmap(self.scale.get(), False, self.app.palette)

    # -- numbers ---------------------------------------------------------------------
    def _image(self, frame=None):
        return kratosmap.blur(kratosmap.pixels(frame or self.frame),
                              round(self.sigma.get(), 2))

    def _colour_range(self):
        import numpy as np
        if self.common.get():
            return snapmap.colour_range(np.concatenate(
                [self._image(f).ravel() for f in self.shown]))
        return snapmap.colour_range(self._image())

    def _sharpness(self, f):
        v = self._sharp.get(id(f.region))
        if v is None:
            v = self._sharp[id(f.region)] = kratosmap.focus_metric(
                kratosmap.pixels(f))
        return v

    def series(self):
        """``(x label, x values, y label, y values)`` of the plot on the right."""
        xl, xs = kratosmap.series_axis(self.shown)
        if self.plot_kind.get() == PLOTS[1]:
            return xl, xs, "Sharpness (a.u.)", [self._sharpness(f)
                                                for f in self.shown]
        what = ("area" if self.mask is not None else "whole image")
        return (xl, xs, f"Mean counts per pixel, {what}",
                kratosmap.roi_means(self.shown, self.mask))

    def series_table(self):
        """Rows for the CSV of the frames: one per shown frame."""
        mean = kratosmap.roi_means(self.shown, self.mask)
        out = []
        for k, f in enumerate(self.shown):
            out.append({
                "frame": k + 1, "map": f.region.name, "position": f.position,
                "binding_energy_eV": f.be, "kinetic_energy_eV": f.ke,
                "stage_z_um": f.z_um,
                "time": f.when.strftime("%Y-%m-%d %H:%M:%S") if f.when else "",
                "dwell_s": f.region.dwell, "mean_counts_per_pixel": mean[k],
                "sharpness": self._sharpness(f)})
        return out

    # -- drawing ---------------------------------------------------------------------
    def _rebuild(self):
        from matplotlib.widgets import RectangleSelector
        pal = self.app.palette
        fig = self.fig
        fig.clf()
        fig.set_facecolor(pal["plot_bg"])
        gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.1], wspace=0.30,
                              left=0.07, right=0.98, top=0.92, bottom=0.12)
        self.ax_map = fig.add_subplot(gs[0])
        self.ax_ser = fig.add_subplot(gs[1])
        for ax in (self.ax_map, self.ax_ser):
            ax.set_facecolor(pal["plot_bg"])
            ax.tick_params(colors=pal["muted"], labelsize=8)
            for sp in ax.spines.values():
                sp.set_color(pal["muted"])
        f = self.frame
        cube = f.cube
        vmin, vmax = self._colour_range()
        left, right, bottom, top = cube.extent()
        self.im = self.ax_map.imshow(
            self._image(), extent=(left, right, bottom, top), origin="upper",
            cmap=self._cmap(), interpolation="nearest", vmin=vmin, vmax=vmax)
        self.ax_map.set_xlim(left, right)
        self.ax_map.set_ylim(bottom, top)
        self.ax_map.set_xlabel("X (µm, approx.)", color=pal["plot_fg"], fontsize=9)
        self.ax_map.set_ylabel("Y (µm, approx.)", color=pal["plot_fg"], fontsize=9)
        self.ax_map.set_title(f.label, fontsize=9, color=pal["plot_fg"])
        self.cbar = fig.colorbar(self.im, ax=self.ax_map, fraction=0.046, pad=0.02)
        self.cbar.ax.tick_params(colors=pal["muted"], labelsize=8)
        self._draw_scale_bar(left, right, bottom, top)
        self._draw_roi_outline()
        self._draw_series()
        self.rect = RectangleSelector(
            self.ax_map, self._on_roi, useblit=False, button=[1], minspanx=0,
            minspany=0, spancoords="data", interactive=True,
            props=dict(edgecolor=pal["accent"], facecolor="none", linewidth=1.6))
        self.cid = self.canvas.mpl_connect("motion_notify_event", self._hover)
        self.hint.configure(text=(
            "Pixel size approximate (calibrated from stage moves, not recorded "
            "in the file).  Drag a box (or click a pixel) for that area's means."))
        self.canvas.draw_idle()

    def _draw_scale_bar(self, left, right, bottom, top):
        pal = self.app.palette
        width = right - left
        bar = _scale_bar(width)
        x1 = left + 0.05 * width
        y = bottom - 0.05 * (bottom - top)               # y runs downward
        self.ax_map.plot([x1, x1 + bar], [y, y], color="white", lw=3,
                         solid_capstyle="butt")
        self.ax_map.text(x1 + bar / 2, y - 0.045 * (bottom - top),
                         f"{bar} µm (approx.)", color="white", fontsize=8,
                         ha="center", va="bottom")

    def _draw_roi_outline(self):
        if self.mask is None or not self.mask.any():
            return
        import numpy as np
        from matplotlib.patches import Rectangle
        rows, cols = np.flatnonzero(self.mask.any(axis=1)), np.flatnonzero(
            self.mask.any(axis=0))
        c = self.frame.cube
        x0, x1 = c.x_of(cols[0]) - c.dx / 2, c.x_of(cols[-1]) + c.dx / 2
        y0, y1 = c.y_of(rows[0]) - c.dy / 2, c.y_of(rows[-1]) + c.dy / 2
        self.ax_map.add_patch(Rectangle(
            (x0, y0), x1 - x0, y1 - y0, fill=False, lw=1.6, ls="--",
            ec=self.app.palette["accent"], zorder=5))

    def _draw_series(self):
        pal = self.app.palette
        ax = self.ax_ser
        xl, xs, yl, ys = self.series()
        # a focus curve is read along Z, whatever order the stage was swept in
        order = (sorted(range(len(xs)), key=lambda k: xs[k])
                 if xl.startswith("Stage Z") else range(len(xs)))
        ax.plot([xs[k] for k in order], [ys[k] for k in order], "-o",
                color=pal["accent"], lw=1.4, ms=4)
        ax.plot([xs[self.i]], [ys[self.i]], "o", ms=9, mfc="none",
                mec=pal["plot_fg"], mew=1.6)
        ax.set_xlabel(xl, color=pal["plot_fg"], fontsize=9)
        ax.set_ylabel(yl, color=pal["plot_fg"], fontsize=9)
        if len(self.shown) < 2:
            ax.text(0.5, 0.5, "one map: nothing to compare", ha="center",
                    va="center", transform=ax.transAxes, color=pal["muted"],
                    fontsize=9)
        ax.set_title(f"{len(self.shown)} map{'s' if len(self.shown) != 1 else ''}"
                     f" — {self.how.get().lower()}", fontsize=9,
                     color=pal["plot_fg"])

    # -- interaction -----------------------------------------------------------------
    def _on_roi(self, click, release):
        cube = self.frame.cube
        if None in (click.xdata, click.ydata, release.xdata, release.ydata):
            return
        x0, y0, x1, y1 = click.xdata, click.ydata, release.xdata, release.ydata
        if abs(x1 - x0) < abs(cube.dx) / 2 and abs(y1 - y0) < abs(cube.dy) / 2:
            hit = cube.pixel_at(x0, y0)                  # a click: one pixel
            if hit is None:
                return
            import numpy as np
            mask = np.zeros((cube.ny, cube.nx), dtype=bool)
            mask[hit[1], hit[0]] = True
        else:
            mask = cube.rect_mask(x0, y0, x1, y1)
        if mask.any():
            self.set_roi(mask)

    def _hover(self, event):
        cube = self.frame.cube
        if event.inaxes is self.ax_map and event.xdata is not None:
            hit = cube.pixel_at(event.xdata, event.ydata)
            if hit:
                v = kratosmap.pixels(self.frame)[hit[1], hit[0]]
                self.hint.configure(
                    text=f"x {event.xdata:.0f} µm   y {event.ydata:.0f} µm "
                         f"(approx.)    pixel ({hit[0]}, {hit[1]})    counts "
                         f"{v:,.0f}")

    # -- saving ------------------------------------------------------------------------
    def _save_menu(self):
        m = tk.Menu(self, tearoff=0)
        self.app.themes.register_menu(m)
        m.add_command(label="Picture (PNG)…", command=self._save_png)
        m.add_command(label="Pixel values of this map (CSV)…",
                      command=self._save_map)
        m.add_command(label="Table of the maps shown (CSV)…",
                      command=self._save_table)
        try:
            m.tk_popup(*self.winfo_pointerxy())
        finally:
            m.grab_release()

    def _stem(self):
        r = self.frame.region
        return f"{r.sample} {r.name} {self.i + 1}".replace(" ", "_")

    def _ask(self, ext, kind):
        return filedialog.asksaveasfilename(
            parent=self, defaultextension=ext, initialfile=self._stem() + ext,
            filetypes=[(kind, "*" + ext)])

    def _save_png(self):
        path = self._ask(".png", "PNG image")
        if path:
            self.fig.savefig(path, dpi=200, facecolor=self.fig.get_facecolor())
            self.app.status.config(text=f"Saved {os.path.basename(path)}")

    def map_csv(self):
        """The current map's counts as CSV text (as recorded, never smoothed)."""
        f = self.frame
        head = (f"# {f.region.sample} {f.region.name}: {f.label}; counts as "
                "recorded; pixel size approximate (not in the file)\n")
        return head + snapmap.to_csv_grid(f.cube, kratosmap.pixels(f))

    def _save_map(self):
        path = self._ask(".csv", "CSV")
        if path:
            with open(path, "w", encoding="utf-8", newline="") as fh:
                fh.write(self.map_csv())
            self.app.status.config(text=f"Saved {os.path.basename(path)}")

    def _save_table(self):
        path = self._ask(".csv", "CSV")
        if not path:
            return
        rows = self.series_table()
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        self.app.status.config(text=f"Saved {os.path.basename(path)}")

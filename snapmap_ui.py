"""The *SnapMap* viewer: where on the sample a peak comes from.

A SnapMap has a whole spectrum at every pixel (see ``snapmap``). Left: the map,
the counts summed over an energy window; right: the spectrum. Drag across the
spectrum to choose the window (the map redraws), drag a box on the map (or click
one pixel) to see the spectrum of just that area next to the whole map's.
Switch element without losing the box, lay the map over the camera image taken at
the same spot, and save the picture, the map values or the area's spectrum.

The dialog talks to the app through ``app.docs``, ``app.palette``, ``app.cfg``,
``app._photo_image`` and ``app.status``. Not modal.
"""

from __future__ import annotations

import csv
import os
import tkinter as tk
from tkinter import filedialog, ttk

import smoothing
import snapmap
import snapmapfig
import snapshot
import themes
from workbook_ui import _finish, keep_in_front


class SnapMapDialog(tk.Toplevel):
    def __init__(self, master, app, parser, region, view=None, figure=None):
        """``view`` (a ``snapmapfig`` dict) restores a saved look; ``figure``
        is the index in ``app.figures`` it came from."""
        super().__init__(master)
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                        NavigationToolbar2Tk)
        from matplotlib.figure import Figure
        self.app, self.parser = app, parser
        self.entries = [r for r in parser.regions
                        if r.extra.get("cube") is not None
                        and r.sample == region.sample] or [region]
        self.region = region
        self.cube = region.extra["cube"]
        self.rois = []                       # [snapmapfig.new_roi dicts], µm
        self.window = snapmap.default_window(self.cube)
        cfg = app.cfg.get("snapmap", {})
        self.scale = tk.StringVar(value=cfg.get("scale", "Viridis")
                                  if cfg.get("scale") in themes.SCALE_NAMES
                                  else "Viridis")
        self.background = tk.BooleanVar(value=bool(cfg.get("background", False)))
        self.overlay = tk.BooleanVar(value=False)
        self.alpha = tk.DoubleVar(value=float(cfg.get("alpha", 0.65)))
        self.smooth_method = tk.StringVar(
            value=cfg.get("smooth_method", "None")
            if cfg.get("smooth_method") in smoothing.METHODS else "None")
        self.smooth_strength = tk.DoubleVar(
            value=float(cfg.get("smooth_strength", 0.5)))
        self.roi_shape = tk.StringVar(value="rect")
        self.mode = tk.StringVar(value="single")     # snapmapfig.MODES key
        # element, window lo / hi for the rows of the composite and ratio modes
        self.layer_rows = [(tk.StringVar(value=""), tk.StringVar(value=""),
                            tk.StringVar(value="")) for _ in range(3)]
        self.ax_spec = self.span = self.im = None
        self.tiles = []                      # [(parser, region)] laid beside this map
        self.figure = figure
        self.title(f"SnapMap — {region.sample}")
        if view:
            self._apply_view(view)

        body = ttk.Frame(self, padding=8)
        body.pack(fill="both", expand=True)
        bar = ttk.Frame(body)
        bar.pack(fill="x")
        ttk.Label(bar, text="Element").pack(side="left")
        self.el_cb = ttk.Combobox(bar, state="readonly", width=14,
                                  values=[self._label(r) for r in self.entries])
        self.el_cb.current(self.entries.index(region)
                           if region in self.entries else 0)
        self.el_cb.pack(side="left", padx=(4, 12))
        self.el_cb.bind("<<ComboboxSelected>>", lambda e: self._element())
        ttk.Label(bar, text="Colours").pack(side="left")
        cb = ttk.Combobox(bar, state="readonly", width=12,
                          textvariable=self.scale, values=themes.SCALE_NAMES)
        cb.pack(side="left", padx=(4, 12))
        cb.bind("<<ComboboxSelected>>", lambda e: self._settings_changed())
        ttk.Checkbutton(bar, text="Remove background", variable=self.background,
                        command=self._settings_changed).pack(side="left")
        self.photo = self._camera_image()
        self.overlay_cb = ttk.Checkbutton(
            bar, text="On camera image", variable=self.overlay,
            command=self._rebuild,
            state="normal" if self.photo is not None else "disabled")
        self.overlay_cb.pack(side="left", padx=(12, 0))
        if getattr(self, "_want_overlay", False) and self.photo is not None:
            self.overlay.set(True)
        ttk.Scale(bar, from_=0.15, to=1.0, variable=self.alpha, length=70,
                  command=lambda v: self._settings_changed(quiet=True)).pack(
            side="left", padx=(6, 0))

        bar2 = ttk.Frame(body)
        bar2.pack(fill="x", pady=(4, 0))
        ttk.Label(bar2, text="View").pack(side="left")
        self.mode_cb = ttk.Combobox(
            bar2, state="readonly", width=14,
            values=[snapmapfig.MODE_NAMES[m] for m in snapmapfig.MODES if m != "fit"])
        self.mode_cb.set(snapmapfig.MODE_NAMES[self.mode.get()])
        self.mode_cb.pack(side="left", padx=(4, 14))
        self.mode_cb.bind("<<ComboboxSelected>>", lambda e: self.set_mode(
            snapmapfig.MODES[self.mode_cb.current()]))
        ttk.Label(bar2, text="Smoothing").pack(side="left")
        smooth_cb = ttk.Combobox(bar2, state="readonly", width=19,
                                 textvariable=self.smooth_method,
                                 values=smoothing.METHODS)
        smooth_cb.pack(side="left", padx=(4, 12))
        smooth_cb.bind("<<ComboboxSelected>>", lambda e: self._smoothing_changed())
        ttk.Label(bar2, text="Strength").pack(side="left")
        ttk.Scale(bar2, from_=0.0, to=1.0, variable=self.smooth_strength,
                  length=110, command=self._smoothing_changed).pack(
            side="left", padx=(4, 0))

        row = ttk.Frame(body)
        row.pack(side="bottom", fill="x")
        panes = ttk.PanedWindow(body, orient="horizontal")
        panes.pack(fill="both", expand=True, pady=(6, 4))
        left = ttk.Frame(panes)
        panes.add(left, weight=5)
        self.fig = Figure(figsize=(10.4, 4.9), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=left)
        tb = ttk.Frame(left)
        tb.pack(side="top", fill="x")
        self.toolbar = NavigationToolbar2Tk(self.canvas, tb, pack_toolbar=False)
        self.toolbar.pack(side="left")
        for name in ("zoom", "pan"):         # a zoom / pan tool owns the mouse
            orig = getattr(self.toolbar, name)

            def wrapped(*a, _orig=orig, **k):
                r = _orig(*a, **k)
                self._sync_tools()
                return r
            setattr(self.toolbar, name, wrapped)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        panes.add(self._side_panel(panes), weight=1)
        self._show_rows()
        self._prefill_rows()
        self.hint = ttk.Label(row, style="Muted.TLabel", text="", anchor="w")
        self.hint.pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Clear areas", command=self._whole_map).pack(
            side="right")
        ttk.Button(row, text="Mosaic…", command=self._mosaic_dialog).pack(
            side="right", padx=6)
        ttk.Button(row, text="Save…", command=self._save_menu).pack(
            side="right", padx=6)
        self.protocol("WM_DELETE_WINDOW", self.destroy)
        self.bind("<Escape>", lambda e: self.destroy())
        _finish(self, app, 1380, 860)
        self.minsize(900, 560)
        self._rebuild()
        keep_in_front(self, master)

    # -- the view as data (saved figures, picture files) -------------------------------
    def view_spec(self, show_spectrum=True):
        """What is on show, as a ``snapmapfig`` view."""
        return snapmapfig.sanitise({
            "ref": self.app.snapmap_ref(self.parser, self.region),
            "window": list(self.window),
            "background": self.background.get(),
            "scale": self.scale.get(),
            "overlay": {"on": self.overlay.get(), "alpha": self.alpha.get()},
            "rois": [dict(r) for r in self.rois],
            "smooth": {"method": self.smooth_method.get(),
                       "strength": self.smooth_strength.get()},
            "show_spectrum": show_spectrum,
            "tiles": [self.app.snapmap_ref(p, r) for p, r in self.tiles],
            "mode": self.mode.get(),
            "layers": self.layers_spec(),
        })

    def _apply_view(self, view):
        """Set the controls from a saved view (before the widgets exist, so
        only variables and fields are touched)."""
        v = snapmapfig.sanitise(view)
        if v["window"]:
            self.window = tuple(v["window"])
        self.background.set(v["background"])
        if v["scale"] in themes.SCALE_NAMES:
            self.scale.set(v["scale"])
        self.alpha.set(v["overlay"]["alpha"])
        self._want_overlay = v["overlay"]["on"]
        if v["smooth"]["method"] in smoothing.METHODS:
            self.smooth_method.set(v["smooth"]["method"])
        self.smooth_strength.set(v["smooth"]["strength"])
        for ref in v["tiles"]:
            hit = self.app.snapmap_target({"ref": ref})
            if hit is not None and hit[1] is not self.region:
                self.tiles.append(hit)
        if v["mode"] != "single" and not self.tiles:
            self.mode.set(v["mode"])
            rows = {"composite": 3, "ratio": 2}.get(v["mode"], 0)
            for i, ly in enumerate(v["layers"][:rows]):
                hit = self.app.snapmap_target({"ref": ly["ref"]})
                if hit is None or hit[1] not in self.entries:
                    continue
                ev, lo, hi = self.layer_rows[i]
                ev.set(self._label(hit[1]))
                if ly["window"]:
                    lo.set(f"{ly['window'][0]:.2f}")
                    hi.set(f"{ly['window'][1]:.2f}")
        self.rois = [] if self.tiles else [
            r for r in v["rois"] if snapmapfig.roi_mask(self.cube, r) is not None]

    # -- single map, grid, composite, ratio ---------------------------------------------
    def _row_count(self):
        return {"composite": 3, "ratio": 2}.get(self.mode.get(), 0)

    def _row_region(self, i):
        """The region chosen in row ``i`` (None for a blank row)."""
        name = self.layer_rows[i][0].get()
        for r in self.entries:
            if self._label(r) == name:
                return r
        return None

    def _row_window(self, i, region):
        ev, lo, hi = self.layer_rows[i]
        try:
            a, b = float(lo.get()), float(hi.get())
            if a != b:
                return (min(a, b), max(a, b))
        except ValueError:
            pass
        return tuple(snapmap.default_window(region.extra["cube"]))

    def _layers(self):
        """``[(region, window)]`` of the multi-element modes (empty in the
        single map). Elements that are not on this map's pixels are left out."""
        mode = self.mode.get()
        if mode == "single":
            return []
        if mode == "grid":
            regs = list(self.entries)[:snapmapfig.MAX_LAYERS]
            wins = [None] * len(regs)
        else:
            regs, wins = [], []
            for i in range(self._row_count()):
                r = self._row_region(i)
                if r is not None:
                    regs.append(r)
                    wins.append(self._row_window(i, r))
        out = []
        for r, w in zip(regs, wins):
            if snapmapfig.same_grid(self.cube, r.extra["cube"]):
                out.append((r, w or tuple(
                    snapmap.default_window(r.extra["cube"]))))
        return out

    def layers_spec(self):
        mode = self.mode.get()
        if mode == "single":
            return []
        return [{"ref": self.app.snapmap_ref(self.parser, r),
                 "window": list(w) if mode != "grid" else None}
                for r, w in self._layers()]

    def effective_mode(self):
        """The mode actually drawn: a multi-element mode needs two usable
        elements and no mosaic."""
        if self.tiles or len(self._layers()) < 2:
            return "single"
        return self.mode.get()

    def set_mode(self, mode):
        """Single map, element grid, RGB composite or ratio map."""
        if mode not in snapmapfig.MODES:
            return
        self.mode.set(mode)
        if mode != "single":
            self.tiles = []
            self.overlay.set(False)
            self._prefill_rows()
        if hasattr(self, "mode_cb"):
            self.mode_cb.set(snapmapfig.MODE_NAMES[mode])
        self._show_rows()
        self._rebuild()
        if mode != "single" and self.effective_mode() == "single":
            self.hint.configure(text="Choose two or more elements (of this "
                                     "map's pixels) for this view.")

    def _prefill_rows(self):
        """Fill empty rows with the file's elements in order, windows blank
        meaning each element's own."""
        for i in range(self._row_count()):
            ev, lo, hi = self.layer_rows[i]
            if ev.get() in ("", "—") and i < len(self.entries):
                ev.set(self._label(self.entries[i]))
        for i in range(self._row_count()):
            ev, lo, hi = self.layer_rows[i]
            r = self._row_region(i)
            if r is not None and not lo.get():
                a, b = snapmap.default_window(r.extra["cube"])
                lo.set(f"{a:.2f}")
                hi.set(f"{b:.2f}")

    def _layer_changed(self, row=None):
        """An element or window box changed: refill that row's default window
        when its element changed, then redraw."""
        if row is not None:
            r = self._row_region(row)
            if r is not None:
                a, b = snapmap.default_window(r.extra["cube"])
                self.layer_rows[row][1].set(f"{a:.2f}")
                self.layer_rows[row][2].set(f"{b:.2f}")
        self._rebuild(keep_view=True)

    def _show_rows(self):
        if not hasattr(self, "layer_frame"):
            return
        names = {"composite": ("Red", "Green", "Blue"),
                 "ratio": ("Numerator", "Denominator", "")}.get(
            self.mode.get())
        if not names:
            self.layer_frame.pack_forget()
            return
        self.layer_frame.pack(fill="x", pady=(0, 8), before=self.area_head)
        for i, w in enumerate(self._row_widgets):
            show = i < self._row_count()
            for wid in w:
                (wid.grid if show else wid.grid_remove)()
            if show:
                w[0].configure(text=names[i])

    # -- what is on show -------------------------------------------------------------
    @staticmethod
    def _label(r):
        return r.name

    def _camera_image(self):
        c = self.cube
        if c.stage_x_mm is None:
            return None
        images = [b for p in self.app.docs for b in p.images]
        return snapshot.nearest_image(images, c.stage_x_mm, c.stage_y_mm)

    def _cmap(self):
        return themes.scale_colourmap(self.scale.get(), False, self.app.palette)

    def _mosaic(self):
        if not self.tiles:
            return None
        return snapmapfig.mosaic_image(
            [self.cube] + [r.extra["cube"] for _p, r in self.tiles],
            self.window, self.background.get())

    def _image(self):
        lo, hi = self.window
        if self.effective_mode() == "ratio":
            return snapmapfig.mode_image(self.view_spec(False), self.region,
                                         self._layers())
        mos = self._mosaic()
        if mos is not None:
            return mos.image
        return self.cube.image(lo, hi, background=self.background.get())

    def _extent(self):
        mos = self._mosaic()
        return mos.extent if mos is not None else self.cube.extent()

    def _save_cfg(self):
        self.app.cfg["snapmap"] = {
            "scale": self.scale.get(),
            "background": self.background.get(),
            "alpha": round(self.alpha.get(), 2),
            "smooth_method": self.smooth_method.get(),
            "smooth_strength": round(self.smooth_strength.get(), 2),
        }

    def _settings_changed(self, quiet=False):
        self._save_cfg()
        if quiet and self.overlay.get() is False:
            return
        self._update_map()

    def _smoothing_changed(self, *_args):
        self._save_cfg()
        self._update_spectrum()

    def _label_with_smoothing(self, base, extra=None):
        """``base`` (+ ``extra`` in parentheses) with the active smoothing
        method appended, so a smoothed trace is never shown unlabelled."""
        method = self.smooth_method.get()
        if method == "None":
            return base if extra is None else f"{base} ({extra})"
        if extra is None:
            return f"{base} ({method})"
        return f"{base} ({extra}, {method})"

    def _element(self):
        self.region = self.entries[self.el_cb.current()]
        self.cube = self.region.extra["cube"]
        self.tiles = []                      # a mosaic is of one core level
        self.window = snapmap.default_window(self.cube)
        self.photo = self._camera_image()
        self.overlay_cb.configure(
            state="normal" if self.photo is not None else "disabled")
        if self.photo is None:
            self.overlay.set(False)
        self.rois = [r for r in self.rois
                     if snapmapfig.roi_mask(self.cube, r) is not None]
        self._rebuild()

    # -- drawing ---------------------------------------------------------------------
    def _rebuild(self, keep_view=False):
        """Draw everything from scratch (element, overlay or areas changed);
        ``keep_view`` leaves a zoomed map as it was."""
        from matplotlib.widgets import RectangleSelector, SpanSelector
        pal = self.app.palette
        fig = self.fig
        kept = None
        if keep_view and getattr(self, "ax_map", None) is not None:
            kept = (self.ax_map.get_xlim(), self.ax_map.get_ylim(),
                    self._view_sig)
        fig.clf()
        fig.set_facecolor(pal["plot_bg"])
        if self.effective_mode() != "single":
            return self._rebuild_layers(kept)
        gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.28,
                              left=0.07, right=0.98, top=0.93, bottom=0.13)
        self.ax_map = fig.add_subplot(gs[0])
        self.ax_spec = fig.add_subplot(gs[1])
        for ax in (self.ax_map, self.ax_spec):
            ax.set_facecolor(pal["plot_bg"])
            ax.tick_params(colors=pal["muted"], labelsize=8)
            for sp in ax.spines.values():
                sp.set_color(pal["muted"])
        cube = self.cube
        # the photo under the map, in the map's own micrometres
        left, right, bottom, top = self._extent()
        if self.overlay.get() and self.photo is not None and not self.tiles:
            img = self.app._photo_image(self.parser, self.photo)
            if img is not None:
                cal = self.photo.calib
                cx, cy = snapshot.stage_to_pixel(cal, cube.stage_x_mm,
                                                 cube.stage_y_mm)
                ux, uy = cal["um_per_px_x"], cal["um_per_px_y"]
                w, h = img.size
                self.ax_map.imshow(img, extent=[-cx * ux, (w - cx) * ux,
                                                (h - cy) * uy, -cy * uy],
                                   zorder=0)
                padx, pady = (right - left) * 0.3, (bottom - top) * 0.3
                self.ax_map.set_xlim(left - padx, right + padx)
                self.ax_map.set_ylim(bottom + pady, top - pady)
        image = self._image()
        vmin, vmax = snapmap.colour_range(image)
        self.im = self.ax_map.imshow(
            image, extent=self._extent(), origin="upper", cmap=self._cmap(),
            interpolation="nearest", vmin=vmin, vmax=vmax, zorder=1,
            alpha=self.alpha.get() if self.overlay.get() else 1.0)
        if not self.overlay.get():
            self.ax_map.set_xlim(left, right)
            self.ax_map.set_ylim(bottom, top)
        self.ax_map.set_xlabel("X (µm)", color=pal["plot_fg"], fontsize=9)
        self.ax_map.set_ylabel("Y (µm)", color=pal["plot_fg"], fontsize=9)
        self.map_title = self.ax_map.set_title("", fontsize=9,
                                               color=pal["plot_fg"])
        self.cbar = fig.colorbar(self.im, ax=self.ax_map, fraction=0.046,
                                 pad=0.02)
        self.cbar.ax.tick_params(colors=pal["muted"], labelsize=8)
        self._view_sig = (id(self.cube), len(self.tiles),
                          self.overlay.get() and self.photo is not None,
                          "single")
        if kept and kept[2] == self._view_sig:
            self.ax_map.set_xlim(kept[0])
            self.ax_map.set_ylim(kept[1])
        self._draw_roi_outline()
        self._draw_tile_outlines()
        self._draw_spectrum()
        self._update_titles()
        self._refresh_areas()

        self.rect = RectangleSelector(
            self.ax_map, self._on_roi, useblit=False, button=[1],
            minspanx=0, minspany=0, spancoords="data", interactive=True,
            props=dict(edgecolor=pal["accent"], facecolor="none", linewidth=1.6))
        self.span = SpanSelector(
            self.ax_spec, self._on_window, "horizontal", useblit=False,
            interactive=True, props=dict(facecolor=pal["accent"], alpha=0.22))
        self.span.extents = self.window
        self.cid = self.canvas.mpl_connect("motion_notify_event", self._hover)
        self.hint.configure(
            text="Drag across the spectrum to choose the energy window; drag a "
                 "box (or click a pixel) on the map to add an area. The "
                 "toolbar zooms and pans.")
        self._sync_tools()
        try:
            self.toolbar.update()            # the zoom history belongs to the old axes
        except Exception:
            pass
        self.canvas.draw_idle()

    def _rebuild_layers(self, kept):
        """Element grid, RGB composite or ratio map: drawn by
        ``snapmapfig.draw_layers`` (the same drawing a saved figure gets); an
        area can still be dragged on a composite or ratio map."""
        from matplotlib.widgets import RectangleSelector
        pal = self.app.palette
        mode = self.effective_mode()
        self.app.draw_snapmap_view(self.fig, (0.0, 0.0, 1.0, 1.0),
                                   self.view_spec(False), pal, self._cmap())
        self.ax_map = self.fig.axes[0]
        self.ax_spec = self.span = self.im = None
        self.cube = self.region.extra["cube"]
        self._view_sig = (id(self.cube), 0, False, mode)
        if kept and kept[2] == self._view_sig and mode != "grid":
            self.ax_map.set_xlim(kept[0])
            self.ax_map.set_ylim(kept[1])
        self.rect = None
        if mode != "grid":
            self.rect = RectangleSelector(
                self.ax_map, self._on_roi, useblit=False, button=[1],
                minspanx=0, minspany=0, spancoords="data", interactive=True,
                props=dict(edgecolor=pal["accent"], facecolor="none",
                           linewidth=1.6))
        self.cid = self.canvas.mpl_connect("motion_notify_event", self._hover)
        self._refresh_areas()
        self.hint.configure(
            text={"grid": "Every element of this map on its own colour range.",
                  "composite": "Each element on its own 1–99 % range, as "
                               "red, green and blue; drag a box for an area.",
                  "ratio": "First element over the second, blank where the "
                           "second is under 5 % of its own 99th percentile."
                  }[mode])
        self._sync_tools()
        try:
            self.toolbar.update()
        except Exception:
            pass
        self.canvas.draw_idle()

    def _sync_tools(self):
        """The zoom / pan tool and the area / window selectors must not both
        take the mouse."""
        busy = bool(getattr(self.toolbar, "mode", ""))
        for sel in (getattr(self, "rect", None), getattr(self, "span", None)):
            if sel is not None:
                sel.set_active(not busy)

    def _draw_tile_outlines(self):
        if not self.tiles:
            return
        from matplotlib.patches import Rectangle
        cubes = [self.cube] + [r.extra["cube"] for _p, r in self.tiles]
        for _i, l, t, w, h in snapmapfig.mosaic_footprints(cubes):
            self.ax_map.add_patch(Rectangle(
                (l, t), w, h, fill=False, lw=0.8, ls=":",
                ec=self.app.palette["muted"], zorder=4))

    def choose_tiles(self, tiles):
        """Lay these ``[(parser, region)]`` beside the current map (empty =
        back to the single map). Uses one energy window and colour range."""
        self.tiles = [(p, r) for p, r in tiles if r is not self.region]
        if self.tiles:
            self.mode.set("single")
            self.mode_cb.set(snapmapfig.MODE_NAMES["single"])
            self._show_rows()
            self.rois = []
            self.overlay.set(False)
        self._rebuild()
        mos = self._mosaic()
        if mos is not None and mos.notes:
            self.hint.configure(text="Mosaic: " + "; ".join(mos.notes) + ".")

    def _mosaic_dialog(self):
        cands = self.app.snapmap_candidates(self.parser, self.region)
        if not cands:
            from tkinter import messagebox
            messagebox.showinfo(
                "Mosaic", f"No other {self.region.name} SnapMap with a stage "
                          "position is loaded.", parent=self)
            return
        win = tk.Toplevel(self)
        win.title("Mosaic of SnapMaps")
        win.transient(self)
        ttk.Label(win, padding=8, wraplength=380, justify="left",
                  text=f"Maps of {self.region.name} to lay beside this one by "
                       "their stage positions (shared energy window and "
                       "colour range).").pack(anchor="w")
        vars_ = []
        chosen = {id(r) for _p, r in self.tiles}
        for p, r, d in cands:
            v = tk.BooleanVar(value=id(r) in chosen)
            vars_.append((v, p, r))
            ttk.Checkbutton(
                win, variable=v,
                text=f"{r.sample}   ({d:.2f} mm away)").pack(anchor="w", padx=12)

        def ok():
            self.choose_tiles([(p, r) for v, p, r in vars_ if v.get()])
            win.destroy()
        row = ttk.Frame(win, padding=8)
        row.pack(fill="x")
        ttk.Button(row, text="Apply", command=ok).pack(side="right")
        ttk.Button(row, text="Cancel", command=win.destroy).pack(
            side="right", padx=6)
        _finish(win, self.app, 420, 120 + 26 * len(cands))

    def masks(self):
        """One boolean mask (or None) per area, in order."""
        return snapmapfig.roi_masks(self.cube, self.rois)

    def _draw_roi_outline(self):
        """Each area on the map, in its own colour (the selector's own box is
        gone after a redraw)."""
        for k, (roi, m) in enumerate(zip(self.rois, self.masks())):
            if m is not None:
                snapmapfig._outline_roi(self.ax_map, self.cube, roi, m,
                                        snapmapfig.roi_colour(k))

    def _mean_spectra(self):
        """The whole-map and each area's spectrum, smoothed for display when
        a smoothing method is chosen. Never touches ``cube.data`` itself."""
        return snapmapfig.mean_spectra(
            self.cube, self.masks(),
            {"method": self.smooth_method.get(),
             "strength": self.smooth_strength.get()})

    def _draw_spectrum(self):
        pal = self.app.palette
        ax = self.ax_spec
        ax.clear()
        ax.set_facecolor(pal["plot_bg"])
        ax.tick_params(colors=pal["muted"], labelsize=8)
        whole, spectra = self._mean_spectra()
        masks = self.masks()
        e = self.cube.energy
        shown = [(k, sp) for k, sp in enumerate(spectra) if sp is not None]
        ax.plot(e, whole, color=pal["muted"] if shown else pal["accent"], lw=1.4,
                label=self._label_with_smoothing("whole map"))
        for k, sp in shown:
            ax.plot(e, sp, color=snapmapfig.roi_colour(k), lw=1.6,
                    label=self._label_with_smoothing(
                        self.rois[k]["label"], f"{int(masks[k].sum())} px"))
        if shown:
            ax.legend(fontsize=8, frameon=False, labelcolor=pal["plot_fg"])
        r = self.region
        binding = "inding" in (r.energy_label or "")
        ax.set_xlim((max(e), min(e)) if binding else (min(e), max(e)))
        ax.set_xlabel(f"{r.energy_label} ({r.energy_units})", color=pal["plot_fg"],
                      fontsize=9)
        ax.set_ylabel("Counts per pixel", color=pal["plot_fg"], fontsize=9)
        ax.set_title(f"{r.name}", fontsize=9, color=pal["plot_fg"])
        for sp in ax.spines.values():
            sp.set_color(pal["muted"])
        # the window, shown even while nothing is being dragged
        lo, hi = self.window
        self.win_patch = ax.axvspan(lo, hi, color=pal["accent"], alpha=0.14,
                                    lw=0)

    def _update_titles(self):
        lo, hi = self.window
        self.map_title.set_text(
            f"{self.region.name}   {min(lo, hi):.1f}–{max(lo, hi):.1f} "
            f"{self.region.energy_units}")

    def _update_map(self):
        """The window, background or colours changed: redraw just the map."""
        if self.effective_mode() != "single":
            self._rebuild(keep_view=True)
            return
        image = self._image()
        vmin, vmax = snapmap.colour_range(image)
        self.im.set_data(image)
        self.im.set_extent(self._extent())
        self.im.set_clim(vmin, vmax)
        self.im.set_cmap(self._cmap())
        self.im.set_alpha(self.alpha.get() if self.overlay.get() else 1.0)
        self._update_titles()
        try:
            self.win_patch.remove()
        except (ValueError, AttributeError):
            pass
        lo, hi = self.window
        self.win_patch = self.ax_spec.axvspan(
            lo, hi, color=self.app.palette["accent"], alpha=0.14, lw=0)
        self._refresh_areas()
        self.canvas.draw_idle()

    def _update_spectrum(self):
        """Smoothing changed: redraw just the spectrum panel (as ``_on_roi``
        already does for a new ROI), leaving the map untouched."""
        if self.ax_spec is None or self.span is None:
            return
        self._draw_spectrum()
        self.span.extents = self.window
        self.canvas.draw_idle()

    # -- interaction -------------------------------------------------------------------
    def _on_window(self, lo, hi):
        if abs(hi - lo) < 1e-9:
            return
        self.window = (min(lo, hi), max(lo, hi))
        self._update_map()

    def _on_roi(self, click, release):
        if self.tiles:                       # an area belongs to one map
            return
        cube = self.cube
        if None in (click.xdata, click.ydata, release.xdata, release.ydata):
            return
        x0, y0, x1, y1 = click.xdata, click.ydata, release.xdata, release.ydata
        if abs(x1 - x0) < abs(cube.dx) / 2 and abs(y1 - y0) < abs(cube.dy) / 2:
            hit = cube.pixel_at(x0, y0)              # a click: one pixel
            if hit is None:
                return
            rect = snapmapfig.roi_rect_of_pixel(cube, *hit)
        else:
            rect = [x0, y0, x1, y1]
        self.add_roi(rect, self.roi_shape.get())

    def add_roi(self, rect, shape="rect", label=""):
        """Add an area (µm rectangle, rectangle or ellipse); returns its
        index, or None when it holds no pixel or the limit is reached."""
        if len(self.rois) >= snapmapfig.MAX_ROIS:
            self.hint.configure(
                text=f"At most {snapmapfig.MAX_ROIS} areas: delete one first.")
            return None
        roi = snapmapfig.new_roi(rect, shape,
                                 label or f"Area {len(self.rois) + 1}")
        if snapmapfig.roi_mask(self.cube, roi) is None:
            return None
        self.rois.append(roi)
        self._rebuild(keep_view=True)
        return len(self.rois) - 1

    def remove_roi(self, i):
        if 0 <= i < len(self.rois):
            del self.rois[i]
            self._rebuild(keep_view=True)

    def rename_roi(self, i, label):
        if 0 <= i < len(self.rois) and label.strip():
            self.rois[i]["label"] = label.strip()
            self._rebuild(keep_view=True)

    def _blank_mask(self):
        import numpy as np
        return np.zeros((self.cube.ny, self.cube.nx), dtype=bool)

    def _whole_map(self):
        self.rois = []
        self._rebuild(keep_view=True)

    # -- the areas' panel and numbers ---------------------------------------------------
    def _side_panel(self, parent):
        side = ttk.Frame(parent, padding=(8, 0, 0, 0))
        self.layer_frame = ttk.LabelFrame(side, text="Elements", padding=6)
        self._row_widgets = []
        labels = ["—"] + [self._label(r) for r in self.entries]
        for i, (ev, lo, hi) in enumerate(self.layer_rows):
            name = ttk.Label(self.layer_frame, text="")
            cb = ttk.Combobox(self.layer_frame, state="readonly", width=9,
                              textvariable=ev, values=labels)
            cb.bind("<<ComboboxSelected>>",
                    lambda e, i=i: self._layer_changed(i))
            e_lo = ttk.Entry(self.layer_frame, textvariable=lo, width=7)
            e_hi = ttk.Entry(self.layer_frame, textvariable=hi, width=7)
            for e in (e_lo, e_hi):
                e.bind("<Return>", lambda ev_: self._layer_changed())
                e.bind("<FocusOut>", lambda ev_: self._layer_changed())
            for c, wid in enumerate((name, cb, e_lo, e_hi)):
                wid.grid(row=i, column=c, padx=2, pady=1, sticky="w")
            self._row_widgets.append((name, cb, e_lo, e_hi))
        ttk.Label(self.layer_frame, style="Muted.TLabel", wraplength=230,
                  text="Energy window of each element (eV, from / to).").grid(
            row=3, column=0, columnspan=4, sticky="w", pady=(4, 0))
        self.area_head = ttk.Label(side, text="Areas", style="Muted.TLabel")
        self.area_head.pack(anchor="w")
        cols = ("area", "px", "um2", "mean", "ratio")
        self.area_tree = ttk.Treeview(side, columns=cols, show="headings",
                                      height=6, selectmode="browse")
        for c, text, w in (("area", "Area", 80), ("px", "Pixels", 50),
                           ("um2", "µm²", 70), ("mean", "Mean", 70),
                           ("ratio", "vs map", 50)):
            self.area_tree.heading(c, text=text)
            self.area_tree.column(c, width=w, anchor="w" if c == "area" else "e")
        self.area_tree.pack(fill="x", pady=(2, 4))
        for k in range(snapmapfig.MAX_ROIS):
            self.area_tree.tag_configure(f"a{k}",
                                         foreground=snapmapfig.roi_colour(k))
        row = ttk.Frame(side)
        row.pack(fill="x")
        ttk.Label(row, text="New area").pack(side="left")
        ttk.Combobox(row, state="readonly", width=9, textvariable=self.roi_shape,
                     values=snapmapfig.SHAPES).pack(side="left", padx=4)
        row2 = ttk.Frame(side)
        row2.pack(fill="x", pady=(4, 0))
        ttk.Button(row2, text="Rename…", command=self._rename_selected).pack(
            side="left")
        ttk.Button(row2, text="Delete", command=self._delete_selected).pack(
            side="left", padx=4)
        ttk.Button(side, text="Save the numbers (CSV)…",
                   command=self._save_stats).pack(anchor="w", pady=(8, 0))
        self.area_note = ttk.Label(side, style="Muted.TLabel", wraplength=230,
                                   justify="left", text="")
        self.area_note.pack(anchor="w", pady=(8, 0))
        return side

    def stats(self):
        """The areas' numbers over the energy window on show (see
        ``snapmapfig.roi_stats``); empty in a mosaic."""
        mode = self.effective_mode()
        if self.tiles or mode in ("grid", "composite"):
            return []
        image = self._image() if mode == "ratio" else None
        return snapmapfig.roi_stats(self.cube, self.window,
                                    self.background.get(), self.rois,
                                    image=image)

    def _refresh_areas(self):
        tree = getattr(self, "area_tree", None)
        if tree is None:
            return
        tree.delete(*tree.get_children())
        mode = self.effective_mode()
        image = self._image() if mode == "ratio" else None
        for k, roi in enumerate(self.rois):
            got = ([] if self.tiles or mode in ("grid", "composite")
                   else snapmapfig.roi_stats(
                       self.cube, self.window, self.background.get(), [roi],
                       image=image))
            if not got:
                continue
            r = got[0]
            tree.insert("", "end", iid=str(k), tags=(f"a{k}",), values=(
                r["label"], r["px"], f"{r['area_um2']:.0f}", f"{r['mean']:.4g}",
                f"{r['ratio']:.2f}"))
        lo, hi = self.window
        self.area_note.configure(
            text=("An area belongs to one map, so there are none in a mosaic."
                  if self.tiles else
                  "Numbers are given for the single map and the ratio map."
                  if mode in ("grid", "composite") else
                  "Mean is per pixel of the ratio image; vs map is its ratio "
                  "to the whole map's mean." if mode == "ratio" else
                  f"Counts summed over {min(lo, hi):.1f}-{max(lo, hi):.1f} "
                  f"{self.region.energy_units}"
                  f"{', background removed' if self.background.get() else ''}; "
                  "Mean is per pixel, vs map its ratio to the whole map's mean."))

    def _selected_area(self):
        sel = self.area_tree.selection()
        return int(sel[0]) if sel else None

    def _rename_selected(self):
        k = self._selected_area()
        if k is None:
            return
        from tkinter import simpledialog
        name = simpledialog.askstring("Rename area", "Name:", parent=self,
                                      initialvalue=self.rois[k]["label"])
        if name:
            self.rename_roi(k, name)

    def _delete_selected(self):
        k = self._selected_area()
        if k is not None:
            self.remove_roi(k)

    def stats_csv(self):
        """The areas' numbers as CSV text."""
        lo, hi = self.window
        out = [f"# {self.region.sample} {self.region.name}: counts summed over "
               f"{min(lo, hi):.3f}-{max(lo, hi):.3f} {self.region.energy_units}"
               f"{', background removed' if self.background.get() else ''}",
               "Area,Shape,Pixels,Area (um2),Mean per pixel,Sum,SD,Mean vs whole map"]
        for r in self.stats():
            out.append(",".join([
                '"' + r["label"].replace('"', '""') + '"', r["shape"],
                str(r["px"]), f"{r['area_um2']:.6g}", f"{r['mean']:.7g}",
                f"{r['sum']:.7g}", f"{r['sd']:.7g}", f"{r['ratio']:.5g}"]))
        return "\n".join(out) + "\n"

    def _save_stats(self):
        if not self.rois or self.tiles:
            self.hint.configure(text="Draw an area on the map first.")
            return
        path = self._ask(".csv", "CSV")
        if path:
            with open(path, "w", encoding="utf-8-sig", newline="") as fh:
                fh.write(self.stats_csv())
            self.app.status.config(text=f"Saved {os.path.basename(path)}")

    def _hover(self, event):
        cube = self.cube
        if self.im is None:
            return
        if event.inaxes is self.ax_map and event.xdata is not None:
            hit = cube.pixel_at(event.xdata, event.ydata)
            if hit:
                v = self.im.get_array()[hit[1], hit[0]]
                self.hint.configure(
                    text=f"x {event.xdata:.0f} µm   y {event.ydata:.0f} "
                         f"µm    pixel ({hit[0]}, {hit[1]})    value "
                         f"{v:,.1f}")
        elif event.inaxes is self.ax_spec and event.xdata is not None:
            self.hint.configure(text=f"{event.xdata:.2f} "
                                     f"{self.region.energy_units}")

    # -- saving ------------------------------------------------------------------------
    def _save_menu(self):
        m = tk.Menu(self, tearoff=0)
        self.app.themes.register_menu(m)
        m.add_command(label="Picture, map and spectrum (PNG, PDF, SVG)…",
                      command=self._save_png)
        m.add_command(label="Picture, map only (PNG, PDF, SVG)…",
                      command=lambda: self._save_png(spectrum=False))
        m.add_separator()
        m.add_command(label="Save as a figure for the report…",
                      command=self._ask_figure)
        if self.figure is not None and self.figure < len(self.app.figures):
            m.add_command(
                label=f"Update figure “{self.app.figures[self.figure]['name']}”",
                command=self.update_figure)
        m.add_separator()
        m.add_command(label="Map values (CSV)…", command=self._save_map)
        m.add_command(label="Spectra (CSV)…", command=self._save_spectra)
        w = self.winfo_pointerxy()
        try:
            m.tk_popup(*w)
        finally:
            m.grab_release()

    def _stem(self):
        return f"{self.region.sample} {self.region.name}".replace(" ", "_")

    def _ask(self, ext, kind):
        return filedialog.asksaveasfilename(
            parent=self, defaultextension=ext, initialfile=self._stem() + ext,
            filetypes=[(kind, "*" + ext)])

    def save_picture(self, path, spectrum=True):
        """Write the view as a picture (format from the extension), drawn
        for paper: the same drawing a report figure gets."""
        size = (10.4, 4.9) if spectrum else (6.4, 5.2)
        ok = self.app.snapmap_picture(self.view_spec(spectrum), path,
                                      size=size)
        if ok:
            self.app.status.config(text=f"Saved {os.path.basename(path)}")
        return ok

    def _save_png(self, spectrum=True):
        path = filedialog.asksaveasfilename(
            parent=self, defaultextension=".png", initialfile=self._stem() + ".png",
            filetypes=[("PNG image", "*.png"), ("PDF", "*.pdf"),
                       ("SVG", "*.svg")])
        if path:
            self.save_picture(path, spectrum)

    def save_as_figure(self, name):
        """Add the view to the workbook's figures (it appears in the report's
        Figures); returns the figure."""
        fig = self.app.add_figure(name.strip(), state=self.view_spec())
        self.figure = self.app.figures.index(fig)
        self.app.status.config(text=f"Saved figure “{fig['name']}” — caption "
                                    f"it under Workbook ▸ Figures")
        return fig

    def _ask_figure(self):
        from tkinter import simpledialog
        name = simpledialog.askstring(
            "Save as a figure", "Name of the figure:", parent=self,
            initialvalue=f"{self.region.sample} {self.region.name}")
        if name and name.strip():
            self.save_as_figure(name)

    def update_figure(self):
        """Replace the look of the figure this window was opened from."""
        if self.figure is None or self.figure >= len(self.app.figures):
            return False
        self.app.figures[self.figure]["state"] = self.view_spec()
        self.app.wb_touch()
        self.app.status.config(text="Figure updated")
        return True

    def _save_map(self):
        path = self._ask(".csv", "CSV")
        if not path:
            return
        with open(path, "w", encoding="utf-8", newline="") as fh:
            lo, hi = self.window
            fh.write(f"# {self.region.sample} {self.region.name}: counts summed "
                     f"over {min(lo, hi):.3f}-{max(lo, hi):.3f} "
                     f"{self.region.energy_units}"
                     f"{', background removed' if self.background.get() else ''}"
                     "\n")
            fh.write(snapmap.to_csv_grid(self.cube, self._image()))
        self.app.status.config(text=f"Saved {os.path.basename(path)}")

    def _csv_smooth_suffix(self):
        """Appended to a CSV column header so a smoothed export never looks
        like raw data."""
        method = self.smooth_method.get()
        return f", {method}" if method != "None" else ""

    def spectra_csv_rows(self):
        """Header and rows of the spectra export: energy, whole map, then one
        column per area."""
        whole, spectra = self._mean_spectra()
        masks = self.masks()
        suffix = self._csv_smooth_suffix()
        head = [f"{self.region.energy_label} ({self.region.energy_units})",
                f"Whole map (counts per pixel{suffix})"]
        cols = []
        for k, sp in enumerate(spectra):
            if sp is not None:
                head.append(f"{self.rois[k]['label']}, {int(masks[k].sum())} px "
                            f"(counts per pixel{suffix})")
                cols.append(sp)
        rows = [[f"{e:.4f}", f"{whole[i]:.6g}"] + [f"{c[i]:.6g}" for c in cols]
                for i, e in enumerate(self.cube.energy)]
        return head, rows

    def _save_spectra(self):
        path = self._ask(".csv", "CSV")
        if not path:
            return
        head, rows = self.spectra_csv_rows()
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(head)
            w.writerows(rows)
        self.app.status.config(text=f"Saved {os.path.basename(path)}")

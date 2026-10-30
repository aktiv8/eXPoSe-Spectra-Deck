"""The settings file: a Workspace that does not own it (the tests, a script)
neither reads nor writes it, and Help > *Reset settings to defaults* puts it
back to the defaults after a warning.

Run:  python -m unittest discover tests
"""

import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

try:
    import tkinter as tk
    import spectradeck as ee
    HAVE_MPL = ee.HAVE_MPL
except Exception:                                   # pragma: no cover
    tk, ee, HAVE_MPL = None, None, False

import splash  # noqa: E402

SETTINGS = {"iss_axis": "Energy ratio", "reels_axis": "Kinetic",
            "show_splash": False, "geometry": "900x600",
            "external_apps": {"casaxps": r"C:\somewhere\casa.exe"}}


class Files(unittest.TestCase):
    """A temporary home: the module's settings paths point into it."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.dir, True)
        self.config = os.path.join(self.dir, "config.json")
        self.legacy = os.path.join(self.dir, "legacy.json")
        self.calib = os.path.join(self.dir, "calib.json")
        self.saved = (ee.CONFIG_PATH, ee.LEGACY_CONFIG_PATH, ee.CALIB_PATH,
                      ee.LEGACY_CALIB_PATH)
        ee.CONFIG_PATH, ee.LEGACY_CONFIG_PATH = self.config, self.legacy
        ee.CALIB_PATH = self.calib
        ee.LEGACY_CALIB_PATH = os.path.join(self.dir, "legacy_calib.json")
        self.addCleanup(self.restore)

    def restore(self):
        (ee.CONFIG_PATH, ee.LEGACY_CONFIG_PATH, ee.CALIB_PATH,
         ee.LEGACY_CALIB_PATH) = self.saved

    def write(self, path, data):
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(data, fh)

    def read(self, path):
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)


@unittest.skipUnless(ee is not None, "spectradeck not importable")
class TestFiles(Files):
    def test_the_reset_leaves_an_empty_settings_file(self):
        self.write(self.config, SETTINGS)
        self.assertEqual(ee.load_config(), SETTINGS)
        self.assertTrue(ee.reset_settings_files())
        self.assertEqual(self.read(self.config), {})
        self.assertEqual(ee.load_config(), {})

    def test_the_former_applications_settings_do_not_come_back(self):
        # with no new file the legacy one is read: the reset writes an empty
        # file instead of deleting, so the legacy settings are not used again
        self.write(self.legacy, {"theme": "Dark"})
        self.assertEqual(ee.load_config(), {"theme": "Dark"})
        self.assertTrue(ee.reset_settings_files())
        self.assertEqual(ee.load_config(), {})
        self.assertEqual(self.read(self.legacy), {"theme": "Dark"})  # untouched

    def test_a_half_written_file_is_removed_and_the_calibration_is_kept(self):
        self.write(self.config, SETTINGS)
        with open(self.config + ".tmp", "w") as fh:
            fh.write("{ half")
        self.write(self.calib, {"scale": 1.0})
        self.assertTrue(ee.reset_settings_files())
        self.assertFalse(os.path.exists(self.config + ".tmp"))
        self.assertEqual(self.read(self.calib), {"scale": 1.0})

    def test_a_file_that_cannot_be_written_says_so_and_changes_nothing(self):
        os.mkdir(self.config)                   # a folder where the file goes
        self.assertFalse(ee.reset_settings_files())
        self.assertTrue(os.path.isdir(self.config))

    def test_the_splash_reads_the_defaults_after_a_reset(self):
        self.write(self.config, SETTINGS)
        self.assertFalse(splash.wanted([], config_path=self.config))
        ee.reset_settings_files()
        self.assertTrue(splash.wanted([], config_path=self.config))
        self.assertFalse(splash.dpi_aware_wanted([], config_path=self.config))


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestWorkspaceOwnership(Files):
    def workspace(self, **kw):
        import matplotlib
        saved = matplotlib.rcParams.copy()
        try:
            root = tk.Tk()
        except tk.TclError:
            self.skipTest("no display")
        root.withdraw()
        self.addCleanup(lambda: (root.destroy(),
                                 matplotlib.rcParams.update(saved)))
        return ee.Workspace(root, **kw)

    def test_a_workspace_does_not_read_the_settings_by_default(self):
        self.write(self.config, SETTINGS)
        seen = []
        real = ee.load_config
        ee.load_config = lambda: seen.append(1) or real()
        self.addCleanup(setattr, ee, "load_config", real)
        ws = self.workspace()
        self.assertFalse(ws.persist)
        self.assertEqual(seen, [])
        self.assertEqual(ws.cfg, {})
        self.assertEqual(ws.iss_axis_var.get(), "Kinetic")       # the default
        self.assertEqual(ws.reels_axis_var.get(), "Energy loss")  # not Kinetic

    def test_a_workspace_does_not_write_the_settings_by_default(self):
        self.write(self.config, SETTINGS)
        ws = self.workspace()
        before = self.read(self.config)
        ws.set_ident_split(True)
        ws.set_external_hidden("casaxps", True)
        self.assertFalse(ws._write_settings())
        self.assertEqual(self.read(self.config), before)

    def test_nothing_is_created_either(self):
        ws = self.workspace()
        ws.set_ident_split(True)
        self.assertFalse(os.path.exists(self.config))

    def test_the_calibration_seed_is_not_read_or_written_by_default(self):
        self.write(self.calib, {"scale": 9.0, "marker": "from the user"})
        ws = self.workspace()
        self.assertIsNone(ws._seed_calibration())
        writes = []
        real = ee.save_calibration
        ee.save_calibration = lambda c: writes.append(c) or True
        self.addCleanup(setattr, ee, "save_calibration", real)
        ws.calib = {"scale": 2.0}
        ws._hide_calib()                      # "the last one seeds new ones"
        self.assertEqual(writes, [])
        self.assertEqual(self.read(self.calib)["scale"], 9.0)

    def test_the_application_reads_and_writes_the_calibration_seed(self):
        self.write(self.calib, {"scale": 9.0})
        ws = self.workspace(persist=True)
        self.assertEqual(ws._seed_calibration(), {"scale": 9.0})
        writes = []
        real = ee.save_calibration
        ee.save_calibration = lambda c: writes.append(c) or True
        self.addCleanup(setattr, ee, "save_calibration", real)
        ws.calib = {"scale": 2.0}
        ws._hide_calib()
        self.assertEqual(writes, [{"scale": 2.0}])

    def test_the_application_reads_and_writes_them(self):
        self.write(self.config, SETTINGS)
        ws = self.workspace(persist=True)
        self.assertTrue(ws.persist)
        self.assertEqual(ws.cfg["geometry"], "900x600")
        self.assertEqual(ws.iss_axis_var.get(), "Energy ratio")
        ws.set_ident_split(True)
        self.assertTrue(self.read(self.config)["ident_split_doublets"])
        self.assertTrue(ws._write_settings())

    def test_main_builds_the_workspace_with_persistence_on(self):
        import inspect
        src = inspect.getsource(ee.main)
        self.assertIn("Workspace(root, persist=True)", src)


@unittest.skipUnless(HAVE_MPL, "matplotlib / Tk not available")
class TestReset(Files):
    def setUp(self):
        super().setUp()
        import matplotlib
        saved = matplotlib.rcParams.copy()
        try:
            self.root = tk.Tk()
        except tk.TclError:
            self.skipTest("no display")
        self.root.withdraw()
        self.addCleanup(lambda: (self.root.destroy(),
                                 matplotlib.rcParams.update(saved)))
        self.write(self.config, SETTINGS)
        self.ws = ee.Workspace(self.root, persist=True)
        self.asked = []
        self.answers = {"askokcancel": True, "askyesno": False}
        mb = ee.messagebox
        self.saved_boxes = (mb.askokcancel, mb.askyesno, mb.showerror,
                            mb.showinfo)

        def ask(name):
            def f(*a, **k):
                self.asked.append((name, a, k))
                return self.answers[name]
            return f
        mb.askokcancel, mb.askyesno = ask("askokcancel"), ask("askyesno")
        mb.showerror = lambda *a, **k: self.asked.append(("showerror", a, k))
        mb.showinfo = lambda *a, **k: None
        self.addCleanup(self.restore_boxes)
        self.closed = []
        self.ws._on_close = lambda: self.closed.append(1)

    def restore_boxes(self):
        mb = ee.messagebox
        (mb.askokcancel, mb.askyesno, mb.showerror,
         mb.showinfo) = self.saved_boxes

    def test_the_warning_says_what_goes_and_what_stays(self):
        self.ws.reset_settings()
        name, args, kw = self.asked[0]
        self.assertEqual(name, "askokcancel")
        text = args[1]
        for token in ("cannot be undone", "theme", "layout", "recent files",
                      "CasaXPS", "plot styles", "report presets",
                      "data files", "workbooks", "calibration seed",
                      "next time the program starts"):
            self.assertIn(token, text)
        self.assertEqual(kw.get("icon"), "warning")
        self.assertEqual(kw.get("default"), "cancel")     # Enter does not reset

    def test_cancel_changes_nothing(self):
        self.answers["askokcancel"] = False
        self.assertFalse(self.ws.reset_settings())
        self.assertEqual(self.read(self.config), SETTINGS)
        self.assertFalse(self.ws._settings_reset)
        self.assertEqual(self.ws.cfg["geometry"], "900x600")
        self.assertEqual(len(self.asked), 1)              # no second question
        self.assertEqual(self.closed, [])

    def test_ok_resets_the_file_and_stops_the_old_settings_being_written_back(self):
        self.assertTrue(self.ws.reset_settings())
        self.assertEqual(self.read(self.config), {})
        self.assertEqual(self.ws.cfg, {})
        self.assertTrue(self.ws._settings_reset)
        # closing, or any change, must not put the old settings back
        self.assertFalse(self.ws._write_settings({"iss_axis": "Energy ratio"}))
        self.ws.set_ident_split(True)
        self.assertEqual(self.read(self.config), {})
        self.assertEqual(ee.load_config(), {})

    def test_it_offers_to_close_and_does_so_only_on_yes(self):
        self.ws.reset_settings()
        self.assertEqual([n for n, _a, _k in self.asked],
                         ["askokcancel", "askyesno"])
        self.assertEqual(self.closed, [])                 # said no
        self.assertIn("defaults", self.asked[1][1][1])
        self.answers["askyesno"] = True
        self.ws._settings_reset = False
        self.write(self.config, SETTINGS)
        self.assertTrue(self.ws.reset_settings())
        self.assertEqual(self.closed, [1])

    def test_a_failure_is_reported_and_nothing_is_marked_reset(self):
        real = ee.reset_settings_files
        ee.reset_settings_files = lambda: False
        self.addCleanup(setattr, ee, "reset_settings_files", real)
        self.assertFalse(self.ws.reset_settings())
        self.assertEqual(self.read(self.config), SETTINGS)
        self.assertFalse(self.ws._settings_reset)
        self.assertEqual(self.ws.cfg["geometry"], "900x600")
        self.assertEqual(self.asked[-1][0], "showerror")
        self.assertEqual(self.closed, [])

    def test_a_workspace_that_does_not_own_the_file_cannot_reset_it(self):
        ws = ee.Workspace(self.root)                      # persist False
        self.assertFalse(ws.reset_settings())
        self.assertEqual(self.asked, [])                  # not even asked
        self.assertEqual(self.read(self.config), SETTINGS)

    def test_the_defaults_apply_to_the_next_workspace(self):
        self.ws.reset_settings()
        self.ws.root.update_idletasks()
        # (a second window is not built here: the file is what matters)
        self.assertEqual(ee.load_config().get("iss_axis", "Kinetic"), "Kinetic")
        self.assertNotIn("external_apps", ee.load_config())

    def help_menu(self):
        bar = self.ws.menubar
        index = next(i for i in range(bar.index("end") + 1)
                     if bar.type(i) == "cascade"
                     and bar.entrycget(i, "label") == "Help")
        return bar.nametowidget(bar.entrycget(index, "menu"))

    def reset_entry(self, menu):
        """The index of the reset entry (separators count in an index)."""
        return next(i for i in range(menu.index("end") + 1)
                    if menu.type(i) == "command"
                    and menu.entrycget(i, "label") == ee.RESET_SETTINGS_LABEL)

    def test_the_menu_has_the_entry_and_it_asks_first(self):
        helpm = self.help_menu()
        i = self.reset_entry(helpm)
        self.assertTrue(ee.RESET_SETTINGS_LABEL.endswith("…"))   # it asks
        self.answers["askokcancel"] = False
        helpm.invoke(i)
        self.assertEqual(self.asked[0][0], "askokcancel")
        self.assertEqual(self.read(self.config), SETTINGS)       # cancelled

    def test_choosing_it_and_confirming_resets(self):
        helpm = self.help_menu()
        helpm.invoke(self.reset_entry(helpm))
        self.assertEqual(self.read(self.config), {})

    def test_the_entry_is_not_dimmed_without_data(self):
        self.assertEqual(self.ws.docs, [])
        helpm = self.help_menu()
        self.assertEqual(helpm.entrycget(self.reset_entry(helpm), "state"),
                         "normal")


if __name__ == "__main__":
    unittest.main()

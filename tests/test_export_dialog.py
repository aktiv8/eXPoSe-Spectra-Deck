"""The Export dialog's CasaXPS options (VAMAS): built, read and passed on."""
import tkinter as tk
import unittest
from unittest import mock

import spectradeck


class FakeApp:
    def __init__(self):
        self.docs = []
        self.checked = set()
        self.cfg = {}
        self.calls = []

    def _write_export(self, chosen, fmt, **kw):
        self.calls.append((chosen, fmt, kw))
        return True

    def _write_settings(self, cfg=None):      # the Workspace's settings seam
        self.settings_written = getattr(self, "settings_written", 0) + 1
        return False


class TestExportDialog(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.root = tk.Tk()
        except tk.TclError:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()

    def dialog(self, app):
        with mock.patch.object(spectradeck.ExportDialog, "grab_set"):
            dlg = spectradeck.ExportDialog(self.root, app)
        dlg.sections = [mock.Mock(selected_regions=lambda: ["r"])]
        return dlg

    def test_defaults_are_off_and_unchanged(self):
        app = FakeApp()
        dlg = self.dialog(app)
        self.assertFalse(dlg.casa_tf.get())
        self.assertFalse(dlg.mfp_on.get())
        self.assertEqual(dlg.mfp_value.get(), "0.6")
        with mock.patch.object(spectradeck, "save_config"):
            dlg.do_export()
        _c, _f, kw = app.calls[0]
        self.assertFalse(kw["casa_transmission"])
        self.assertIsNone(kw["escape_depth"])

    def test_options_reach_the_export(self):
        app = FakeApp()
        dlg = self.dialog(app)
        dlg.casa_tf.set(True)
        dlg.mfp_on.set(True)
        dlg.mfp_value.set("0,6")
        dlg.mfp_type.set(spectradeck.ESCAPE_TYPES[1])
        with mock.patch.object(spectradeck, "save_config"):
            dlg.do_export()
        kw = app.calls[0][2]
        self.assertTrue(kw["casa_transmission"])
        # typed 0,6 -> written -0.6, the value CasaXPS stores for 1/KE^0.6
        self.assertEqual(kw["escape_depth"], {"exponent": -0.6, "type": 1})
        self.assertEqual(app.cfg["export_mfp_value"], -0.6)

    def test_sign_typed_does_not_matter(self):
        for typed in ("0.6", "-0.6"):
            app = FakeApp()
            dlg = self.dialog(app)
            dlg.mfp_on.set(True)
            dlg.mfp_value.set(typed)
            with mock.patch.object(spectradeck, "save_config"):
                dlg.do_export()
            self.assertEqual(app.calls[0][2]["escape_depth"]["exponent"], -0.6)

    def test_every_casaxps_type_is_offered(self):
        self.assertEqual(sorted(spectradeck.ESCAPE_TYPES), [0, 1, 2, 3, 4, 5, 6])

    def test_thickness_types_write_it_positive(self):
        for kind in (4, 6):
            app = FakeApp()
            dlg = self.dialog(app)
            dlg.mfp_on.set(True)
            dlg.mfp_type.set(spectradeck.ESCAPE_TYPES[kind])
            dlg.mfp_value.set("-1,5")
            with mock.patch.object(spectradeck, "save_config"):
                dlg.do_export()
            self.assertEqual(app.calls[0][2]["escape_depth"],
                             {"exponent": 1.5, "type": kind})
            self.assertEqual(app.cfg["export_escape_thickness"], 1.5)
            self.assertNotIn("export_mfp_value", app.cfg)

    def test_type_without_a_value_keeps_the_exponent_form(self):
        app = FakeApp()
        dlg = self.dialog(app)
        dlg.mfp_on.set(True)
        dlg.mfp_type.set(spectradeck.ESCAPE_TYPES[2])
        self.assertEqual(str(dlg.mfp_entry.cget("state")), "disabled")
        with mock.patch.object(spectradeck, "save_config"):
            dlg.do_export()
        self.assertEqual(app.calls[0][2]["escape_depth"],
                         {"exponent": -0.6, "type": 2})

    def test_bad_exponent_is_refused(self):
        app = FakeApp()
        dlg = self.dialog(app)
        dlg.mfp_on.set(True)
        dlg.mfp_value.set("abc")
        with mock.patch.object(spectradeck.messagebox, "showwarning") as w:
            dlg.do_export()
        self.assertTrue(w.called)
        self.assertEqual(app.calls, [])


if __name__ == "__main__":
    unittest.main()

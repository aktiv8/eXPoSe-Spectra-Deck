"""Uncaught errors are logged and shown, never lost (crashlog.py).

Run:  python -m unittest discover tests
"""

import io
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import appinfo  # noqa: E402
import crashlog  # noqa: E402


def boom():
    try:
        raise ValueError("the spectrum has no counts")
    except ValueError:
        return sys.exc_info()


class TestReport(unittest.TestCase):
    def test_the_report_has_the_version_and_the_traceback(self):
        text = crashlog.format_report(*boom())
        self.assertIn(f"{appinfo.NAME} {appinfo.VERSION}", text)
        self.assertIn("Traceback (most recent call last)", text)
        self.assertIn("ValueError: the spectrum has no counts", text)
        self.assertIn("boom", text)                      # the frame

    def test_the_short_message_is_one_trimmed_line(self):
        msg = crashlog.short_message(ValueError("a\n  b " + "x" * 400))
        self.assertNotIn("\n", msg)
        self.assertTrue(msg.startswith("ValueError: a b"))
        self.assertLessEqual(len(msg), 300)


class TestLog(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.path = os.path.join(self.dir, "crash.log")
        self.addCleanup(self._clean)

    def _clean(self):
        crashlog.teardown()
        for n in os.listdir(self.dir):
            try:
                os.remove(os.path.join(self.dir, n))
            except OSError:
                pass
        os.rmdir(self.dir)

    def test_an_error_is_written_to_the_log(self):
        self.assertTrue(crashlog.setup(self.path))
        self.assertTrue(crashlog.setup(self.path))        # idempotent
        crashlog.record(*boom())
        with open(self.path, encoding="utf-8") as fh:
            text = fh.read()
        self.assertEqual(text.count("uncaught exception"), 1)
        self.assertIn("the spectrum has no counts", text)

    def test_an_unwritable_place_is_not_fatal(self):
        bad = os.path.join(self.dir, "no", "such", "folder", "x.log")
        self.assertFalse(crashlog.setup(bad))
        crashlog.record(*boom())                          # must not raise


class QuietStderr:
    """The handler echoes to stderr (a console still gets it); hide that."""

    def setUp(self):
        self._stderr, sys.stderr = sys.stderr, io.StringIO()
        self.addCleanup(self._restore_stderr)
        # never write to the real ~/.spectradeck.log
        self._logdir = tempfile.mkdtemp()
        patch = mock.patch.object(crashlog, "LOG_PATH",
                                  os.path.join(self._logdir, "t.log"))
        patch.start()
        self.addCleanup(self._clean_log, patch)

    def _clean_log(self, patch):
        patch.stop()
        crashlog.teardown()
        shutil.rmtree(self._logdir, True)

    def _restore_stderr(self):
        sys.stderr = self._stderr

class TestHandler(QuietStderr, unittest.TestCase):
    def test_the_handler_logs_then_shows_the_short_message(self):
        seen = []
        h = crashlog.handler("root", lambda *a: seen.append(a))
        h(*boom())
        self.assertEqual(len(seen), 1)
        root, report, message, path = seen[0]
        self.assertEqual(root, "root")
        self.assertIn("Traceback", report)
        self.assertTrue(message.startswith("ValueError:"))
        self.assertEqual(path, crashlog.log_path())

    def test_a_failing_dialog_does_not_raise(self):
        def broken(*a):
            raise RuntimeError("no display")
        crashlog.handler(None, broken)(*boom())

    def test_interrupts_are_left_to_python(self):
        seen = []
        h = crashlog.handler(None, lambda *a: seen.append(a))
        old = sys.__excepthook__
        calls = []
        sys.__excepthook__ = lambda *a: calls.append(a)
        try:
            h(KeyboardInterrupt, KeyboardInterrupt(), None)
        finally:
            sys.__excepthook__ = old
        self.assertEqual((seen, len(calls)), ([], 1))


class TestInTk(QuietStderr, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import tkinter as tk
            cls.root = tk.Tk()
        except Exception:
            raise unittest.SkipTest("no display")
        cls.root.withdraw()

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()

    def test_a_raising_callback_reaches_the_handler(self):
        seen = []
        root = self.root
        old = root.report_callback_exception
        old_hook = sys.excepthook
        try:
            crashlog.install(root, lambda *a: seen.append(a))
            root.after(0, lambda: 1 / 0)
            root.after(30, root.quit)
            root.mainloop()
        finally:
            root.report_callback_exception = old
            sys.excepthook = old_hook
            crashlog.teardown()
        self.assertEqual(len(seen), 1)
        self.assertIn("ZeroDivisionError", seen[0][2])

    def test_the_dialog_opens_once_and_offers_the_details(self):
        win = crashlog.show_dialog(self.root, "report text", "ValueError: x",
                                   "log.txt")
        self.assertIsNotNone(win)
        try:
            self.assertIsNone(crashlog.show_dialog(self.root, "again", "y"))
            win.copy_details()
            self.assertEqual(win.clipboard_get(), "report text")
        finally:
            win.close_dialog()
        win2 = crashlog.show_dialog(self.root, "r", "m")       # free again
        self.assertIsNotNone(win2)
        win2.close_dialog()


if __name__ == "__main__":
    unittest.main()

"""What happens to an exception nobody caught: it is written to a log file in
the user's home folder and shown in one plain dialog (with the details to copy
into a bug report) instead of vanishing on a console the user never sees.

Everything but ``show_dialog`` is Tk-free, and nothing here may raise: a
failing log must never turn one error into two.
"""

from __future__ import annotations

import logging
import os
import sys
import traceback
from logging.handlers import RotatingFileHandler

import appinfo

LOG_PATH = os.path.join(os.path.expanduser("~"), ".spectradeck.log")
MAX_BYTES = 200_000
BACKUPS = 3
LOGGER = "spectradeck.crash"

_shown = False                  # one dialog at a time: a looping error only logs


def log_path():
    return LOG_PATH


def setup(path=None):
    """Attach the rotating log file (idempotent). True when it will be written."""
    path = path or LOG_PATH
    try:
        log = logging.getLogger(LOGGER)
        log.setLevel(logging.INFO)
        log.propagate = False
        for h in log.handlers:
            if getattr(h, "baseFilename", None) == os.path.abspath(path):
                return True
        handler = RotatingFileHandler(path, maxBytes=MAX_BYTES,
                                      backupCount=BACKUPS, encoding="utf-8")
        handler.setFormatter(logging.Formatter("%(asctime)s %(message)s"))
        log.addHandler(handler)
        return True
    except Exception:           # noqa: BLE001 - an unwritable home is not fatal
        return False


def teardown():
    """Close and detach the log handlers (tests; releases the file on Windows)."""
    log = logging.getLogger(LOGGER)
    for h in list(log.handlers):
        try:
            h.close()
        except Exception:       # noqa: BLE001
            pass
        log.removeHandler(h)


def short_message(exc):
    """One line for the dialog: the exception's type and text, trimmed."""
    text = f"{type(exc).__name__}: {exc}".strip()
    text = " ".join(text.split())
    return text if len(text) <= 300 else text[:297] + "…"


def format_report(exc_type, exc, tb):
    """Plain text for a bug report: the version block, then the traceback."""
    try:
        head = appinfo.version_report()
    except Exception:           # noqa: BLE001
        head = f"{appinfo.NAME} {appinfo.VERSION}"
    try:
        trace = "".join(traceback.format_exception(exc_type, exc, tb))
    except Exception:           # noqa: BLE001
        trace = f"{exc_type.__name__}: {exc}\n"
    return f"{head}\n\n{trace}".rstrip() + "\n"


def record(exc_type, exc, tb):
    """Log an uncaught exception; returns the report text. Never raises."""
    report = format_report(exc_type, exc, tb)
    try:
        logging.getLogger(LOGGER).error("uncaught exception\n%s", report)
    except Exception:           # noqa: BLE001
        pass
    return report


def show_dialog(root, report, message, path=None):
    """The "Something went wrong" window. Returns it (None when one is already
    open or Tk cannot show it)."""
    global _shown
    if _shown:
        return None
    try:
        import tkinter as tk
        from tkinter import ttk
        _shown = True
        win = tk.Toplevel(root)
        win.title(f"{appinfo.NAME}: something went wrong")
        win.transient(root)
        body = ttk.Frame(win, padding=14)
        body.pack(fill="both", expand=True)
        ttk.Label(body, text="Something went wrong.",
                  font=("", 11, "bold")).pack(anchor="w")
        ttk.Label(body, text="You can carry on working; if it keeps "
                  "happening, copy the details and send them with a "
                  "description of what you were doing.", wraplength=460,
                  justify="left").pack(anchor="w", pady=(4, 8))
        ttk.Label(body, text=message, wraplength=460, justify="left"
                  ).pack(anchor="w")
        if path:
            ttk.Label(body, text=f"Saved in {path}", wraplength=460,
                      justify="left").pack(anchor="w", pady=(6, 0))
        row = ttk.Frame(body)
        row.pack(anchor="e", pady=(12, 0))

        def copy():
            try:
                win.clipboard_clear()
                win.clipboard_append(report)
                copy_btn.config(text="Copied")
            except Exception:   # noqa: BLE001
                pass

        def close():
            global _shown
            _shown = False
            win.destroy()

        copy_btn = ttk.Button(row, text="Copy details", command=copy)
        copy_btn.pack(side="left", padx=(0, 8))
        ttk.Button(row, text="Close", command=close).pack(side="left")
        win.protocol("WM_DELETE_WINDOW", close)
        win.bind("<Escape>", lambda e: close())
        win.copy_details, win.close_dialog = copy, close
        return win
    except Exception:           # noqa: BLE001 - the log already has it
        _shown = False
        return None


def handler(root, show=show_dialog):
    """``f(exc_type, exc, tb)``: log, then show the dialog. Interrupts pass
    through to Python's own hook."""
    def handle(exc_type, exc, tb):
        if issubclass(exc_type, (KeyboardInterrupt, SystemExit)):
            sys.__excepthook__(exc_type, exc, tb)
            return
        report = record(exc_type, exc, tb)
        try:
            sys.stderr.write(report)       # a console still gets it
        except Exception:       # noqa: BLE001 - e.g. no console under pythonw
            pass
        try:
            show(root, report, short_message(exc), LOG_PATH)
        except Exception:       # noqa: BLE001
            pass
    return handle


def install(root, show=show_dialog):
    """Route uncaught exceptions (Tk callbacks and everything else) here."""
    setup()
    h = handler(root, show)
    root.report_callback_exception = h
    sys.excepthook = h
    return h

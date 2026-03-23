"""GUI file/directory chooser dialogs (tkinter wrappers).

On macOS, running tkinter in the main process hooks into the system AppKit
event loop and can cause other applications to stall (spinning beach ball).
To avoid this, each dialog is spawned in a short-lived subprocess so the
tkinter event-loop interaction is fully isolated from the main process.
"""

import subprocess
import sys
import textwrap


def _run_dialog(script: str) -> str:
    """Run *script* in a fresh Python subprocess and return stdout (stripped)."""
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def openFile(initial_dir=None):
    """Open a native file-chooser dialog and return the selected path (or '' if cancelled)."""
    script = textwrap.dedent(f"""
        import tkinter as _tk
        import tkinter.filedialog as _fdialog
        root = _tk.Tk()
        root.withdraw()
        path = _fdialog.askopenfilename(initialdir={initial_dir!r})
        root.destroy()
        print(path)
    """)
    return _run_dialog(script)


def openDir(initial_dir=None):
    """Open a native directory-chooser dialog and return the selected path (or '' if cancelled)."""
    script = textwrap.dedent(f"""
        import tkinter as _tk
        import tkinter.filedialog as _fdialog
        root = _tk.Tk()
        root.withdraw()
        path = _fdialog.askdirectory(initialdir={initial_dir!r})
        root.destroy()
        print(path)
    """)
    return _run_dialog(script)

"""GUI file/directory chooser dialogs (tkinter wrappers).

On macOS, running tkinter in the main process hooks into the system AppKit
event loop and can cause other applications to stall (spinning beach ball).
To avoid this, each dialog is spawned in a short-lived subprocess so the
tkinter event-loop interaction is fully isolated from the main process.
"""

import subprocess
import sys
import textwrap

# Sentinel written by the subprocess to mark the dialog result line.
# Using a prefix that cannot appear in a file-system path ensures we can
# reliably extract the result even when libraries (GTK, libpng, etc.) print
# warnings to stdout before the dialog runs.
_SENTINEL = "__DIALOG_RESULT__:"


def _run_dialog(script: str) -> str:
    """Run *script* in a fresh Python subprocess and return the dialog result.

    The subprocess must write its result using the sentinel prefix so that
    any incidental library output (GTK warnings, libpng messages, etc.) that
    ends up on stdout is ignored.
    """
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
    )
    for line in result.stdout.splitlines():
        if line.startswith(_SENTINEL):
            return line[len(_SENTINEL):]
    return ''


def openFile(initial_dir=None):
    """Open a native file-chooser dialog and return the selected path (or '' if cancelled)."""
    script = textwrap.dedent(f"""
        import tkinter as _tk
        import tkinter.filedialog as _fdialog
        import sys
        root = _tk.Tk()
        root.withdraw()
        path = _fdialog.askopenfilename(initialdir={initial_dir!r})
        root.destroy()
        sys.stdout.write({_SENTINEL!r} + (path or '') + '\\n')
        sys.stdout.flush()
    """)
    return _run_dialog(script)


def openDir(initial_dir=None):
    """Open a native directory-chooser dialog and return the selected path (or '' if cancelled)."""
    script = textwrap.dedent(f"""
        import tkinter as _tk
        import tkinter.filedialog as _fdialog
        import sys
        root = _tk.Tk()
        root.withdraw()
        path = _fdialog.askdirectory(initialdir={initial_dir!r})
        root.destroy()
        sys.stdout.write({_SENTINEL!r} + (path or '') + '\\n')
        sys.stdout.flush()
    """)
    return _run_dialog(script)

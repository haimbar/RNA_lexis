#!/usr/bin/env python3
"""Regenerate docs/user_guide.html and docs/user_guide.pdf from docs/user_guide.md.

Run this after editing user_guide.md, before committing, so the checked-in
HTML/PDF snapshots don't go stale. The HTML step reuses the exact same
rendering function the app itself calls from the "Open User Guide" menu
item (rna_lexis.menu.render_user_guide_html), so there is only ever one
copy of the HTML template to keep correct. The PDF has no such built-in
path -- this script is what keeps it from silently drifting out of sync,
which is exactly what happened before (see CHANGELOG / git history).

Usage:
    python docs/build_docs.py

Requires RNA_lexis to be installed (e.g. `pip install -e .`) and, for the
PDF step, the `weasyprint` command-line tool -- not a runtime dependency
of RNA_lexis itself:
    pip install weasyprint
"""

import shutil
import subprocess
import sys
from pathlib import Path

try:
    from rna_lexis.menu import render_user_guide_html
except ImportError:
    print("Could not import rna_lexis. Install the package first, e.g.:\n"
          "    pip install -e .", file=sys.stderr)
    raise SystemExit(1)


def main() -> int:
    docs_dir = Path(__file__).resolve().parent
    guide_md = docs_dir / "user_guide.md"
    guide_html = docs_dir / "user_guide.html"
    guide_pdf = docs_dir / "user_guide.pdf"

    html_path = render_user_guide_html(str(guide_md), str(guide_html))
    if not html_path:
        print(f"Error: {guide_md} not found.", file=sys.stderr)
        return 1
    print(f"Wrote {guide_html}")

    weasyprint = shutil.which("weasyprint")
    if weasyprint is None:
        print("weasyprint not found on PATH -- skipping PDF.\n"
              "Install it with: pip install weasyprint", file=sys.stderr)
        return 0

    # -p / --presentational-hints: without this, legacy HTML attributes like
    # <img width="200"> are ignored and the logo renders at full native size.
    result = subprocess.run(
        [weasyprint, "-p", "--optimize-images", str(guide_html), str(guide_pdf)],
        capture_output=True, text=True,
    )
    if result.stderr:
        print(result.stderr.strip(), file=sys.stderr)
    if result.returncode != 0:
        print("weasyprint failed.", file=sys.stderr)
        return 1
    print(f"Wrote {guide_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

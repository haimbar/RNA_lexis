# CLAUDE.md

Repo-specific operational notes for agent sessions working on RNA_lexis.
For user-facing docs see `README.md` (install/features/troubleshooting),
`docs/user_guide.md` (interactive-menu walkthrough), and `llms.txt` (library
API reference) — this file is about *working on* the repo, not using it.

## Running things

- Tests: `python3 -m pytest -q` (plain `python` may resolve to a Python
  without pytest installed in some environments — use `python3` explicitly).
  59 tests across `test_statistical.py`, `test_algorithms.py`, `test_io.py`,
  `test_alignment.py`, and `test_batch_cli.py` (the `rna_lexis_stat_cli`
  subcommands) — all synthetic-data unit/end-to-end tests, no real
  biological data or network calls. `menu.py` and `plots.py` still have
  **no automated tests** (interactive/visual) — changes there need
  manual/smoke verification (e.g. `python3 -c "from rna_lexis.X import Y;
  ..."`, or rendering a menu with `gen_menu(..., clr=False)` to inspect
  ANSI output directly).
- No linter/formatter/type-checker is configured (no flake8/mypy/ruff/
  pre-commit config, no `[tool.*]` sections in `pyproject.toml`). Match
  existing style by hand.
- Entry points (`pyproject.toml`): `rna_lexis` and `rna_lexis_stat` both map
  to `rna_lexis.menu:main` (the interactive TUI); `rna_lexis_stat_cli` maps
  to `rna_lexis.test_cli:main` (headless batch CLI — this is the one safe to
  drive non-interactively/in CI, since it never calls the auto-open-file
  helpers the TUI does).

## Documentation workflow

- After editing `docs/user_guide.md`, run `python docs/build_docs.py` to
  regenerate `docs/user_guide.html` and `docs/user_guide.pdf` before
  committing. The HTML step reuses `rna_lexis.menu.render_user_guide_html()`
  — the same function the app calls from its own *Open User Guide* menu item
  — so there is exactly one HTML template, not two to keep in sync. The PDF
  step shells out to `weasyprint -p --optimize-images` (the `-p` /
  `--presentational-hints` flag is required — without it, legacy HTML
  attributes like `<img width="200">` are ignored and images render at full
  native size). `weasyprint` is a dev-only tool, not a runtime dependency.
- `PLAN_paper_figures.md` at the repo root is **local planning material
  only** — listed in `.gitignore`, never committed, not part of shipped
  docs. Don't treat it as authoritative or try to "fix" its exclusion.
- Keep README.md / docs/user_guide.md / llms.txt cross-links intact if you
  restructure any of them (README → user_guide.md + llms.txt; user_guide.md
  → llms.txt + README's Troubleshooting; llms.txt → both).

## Module layout (`src/rna_lexis/`)

`algorithms.py` (pure computation, no I/O) → `alignment.py` (Gotoh
global/local) → `statistical.py` (Markov/FDR scoring) → `plots.py`
(matplotlib + Plotly rendering) → `io.py` (file/session I/O, Ensembl fetch,
platform-specific file-opening) → `dialogs.py` (tkinter file/dir pickers,
each spawned in an isolated subprocess) → `menu.py` (interactive TUI,
imports from all of the above) → `test_cli.py` (headless batch CLI) →
`RNAlang.py` (backward-compat shim re-exporting everything). New code should
import from the specific sub-module, not `RNAlang`.

## Known gotchas / established conventions

- **`kaleido` is pinned `>=0.2,<1` deliberately** (see CHANGELOG 0.1.11) —
  `kaleido>=1` needs a separate Chrome install and only works with
  Plotly ≥ 6. If `fig.write_image()` ever raises with a kaleido/plotly
  version-mismatch warning in a dev environment, the fix is
  `pip install "kaleido>=0.2,<1"`, not `pip install -U kaleido` (which
  installs the incompatible version — kaleido's own error message
  recommends the wrong fix).
- **All file `open()` calls should pass `encoding='utf-8'` explicitly.**
  Windows defaults to the locale codepage otherwise, which breaks on any
  non-ASCII content (confirmed concretely: `docs/user_guide.md` contains
  ~23 distinct non-ASCII characters — arrows, Greek letters, math symbols).
- **`open_file_with_default_software()` (io.py) is self-guarded** — wrapped
  in its own try/except, matching `open_pdf()`. Follow that pattern for any
  new "auto-open a saved file" code; don't let `subprocess.run(check=True)`
  or `os.startfile()` propagate uncaught.
- A `UserWarning: Unable to import Axes3D... multiple versions of
  Matplotlib` during `pytest` in this dev environment is benign (system
  `apt` matplotlib + pip matplotlib coexisting) — not a real bug, don't
  chase it.
- `networkx` and `Levenshtein` were removed from `pyproject.toml`
  dependencies — confirmed unused anywhere in the codebase. If either is
  ever actually needed again, re-add deliberately rather than assuming they
  were left in place.
- **There are two separate version strings — bump both together.**
  `pyproject.toml`'s `version` is the packaging source of truth.
  `rna_lexis/__init__.py` has its own hardcoded `__version__ = "..."` that
  does **not** derive from it and must be edited separately. Neither of
  these is what the TUI banner actually displays, though: `menu.py` reads
  `importlib.metadata.version("RNA_lexis")` — the currently *installed*
  package's metadata — which stays stale until the environment is
  reinstalled. Concretely hit this at 0.1.12: both files correctly said
  0.1.12, but this dev environment's editable install still reported 0.1.9
  (via `pip show` and the on-screen banner) until re-running
  `pip install -e . --no-deps` (use `--no-deps` to avoid disturbing an
  already-correct pinned environment, e.g. the `kaleido` pin above).

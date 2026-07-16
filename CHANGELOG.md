# Changelog

## [0.1.14] - 2026-07-16

### Added

- **Bundled NORAD example dataset** — human (`ENST00000565493`, 5401 nt) and
  mouse (`ENSMUST00000192863`, 4945 nt) NORAD transcripts now ship inside
  the pip distribution itself (`src/rna_lexis/data/`, wired via
  `[tool.setuptools.package-data]`), so a fresh `pip install` gives every
  user a ready-to-try real dataset with no fetching, pasting, or account
  needed. Access programmatically via the new
  `rna_lexis.io.example_dataset_path(name)` (`importlib.resources`-based —
  correct whether installed from a wheel, editable install, or a repo
  checkout). Verified end-to-end with a real wheel build installed into a
  throwaway venv, not just the editable dev checkout.
- **"Load example dataset (NORAD)" menu option** — new item in the Choose
  Input menu (position 4; "Python prompt"/"Quit" renumbered to 5), loads
  the bundled sequence instantly then prompts for a save directory,
  mirroring the existing ENST-fetch flow. The FASTA-header-stripping logic
  previously inlined in `choose_file()` is now a shared
  `_load_fasta_or_text()` helper used by both.
- **`tests/test_norad_integration.py`** (8 tests) — end-to-end integration
  tests running the real discovery pipeline (`find_boundary`/`cores`) on
  real biological data, distinct from T.1's synthetic unit tests. Exactly
  reproduces the published RNA-Lexis paper's own validation numbers for
  this transcript (163 xmotifs, 235 cores with default settings; the PRE
  motif `tgtatata` at exactly 10 exact occurrences and 4 approximate
  `tgtaaata` hits under 1-mismatch search) — the deterministic
  literal-search counts are asserted exactly, while discovery-pipeline
  counts (which could shift with future algorithm tuning) are asserted as
  loose sanity bounds rather than pinned exact values, to avoid brittleness.
  Also confirms human NORAD shows more organizational complexity than
  mouse, matching the paper's own comparative finding.
- Documented two related, pre-existing limitations discovered while
  building this: `rna_lexis.io.read_text()` is not FASTA-aware (a `>`
  header line's words get merged into the sequence — use
  `_load_fasta_or_text()` or strip the header yourself), and
  `fetch_enst_cdna()` only accepts human `ENST`-prefixed transcript IDs
  (the mouse dataset was fetched via the lower-level, species-agnostic
  `_fetch_ensembl_json()` instead). See `CLAUDE.md` for details.

## [0.1.13] - 2026-07-16

### Added

- **Unit test coverage expanded from one module to four** — the suite
  previously covered only `rna_lexis.statistical` (6 tests, partial). Added,
  all using small synthetic sequences with known ground truth (no new
  dependencies, no real biological data — see the separate NORAD dataset
  task for that):
  - `tests/test_algorithms.py` (25 tests) — the core motif-discovery engine
    (`find_all_matches`, `find_boundary`/`cores`, `is_bounded`/
    `expand_to_boundary`, `find_with_mutations`, `extend_match_pair`/
    `find_longest_extensions`, `gen_hairpins`, `cover`,
    `compute_default_wd`, `allow_mutation`, `markov_kmer_pvalues`, plus
    smaller helpers), previously **completely untested** despite being the
    most foundational module in the package.
  - `tests/test_io.py` (14 tests) — session save/load round-trips
    (including a non-ASCII regression test for the 0.1.12 encoding fix),
    `is_valid_session`, `read_text`, `load_prefs`/`save_prefs` round-trips
    (including the 0.1.12 config-directory-failure fallback path),
    `_find_valid_sessions`, the summary-CSV hash cache, and a regression
    test locking in 0.1.12's `open_file_with_default_software()` fix (must
    not raise on a bogus path).
  - `tests/test_alignment.py` (8 tests) — `gotoh_global`, `gotoh_local`,
    `make_markers`.
  - `tests/test_batch_cli.py` (6 tests) — end-to-end tests for all four
    `rna_lexis_stat_cli` subcommands (`score-exact`, `rank-cores`,
    `mutation-families`, `gapped-motif`), previously untested despite being
    the most test-friendly entry point in the package (argparse-based,
    file-in/file-out, no interactive prompts).

  Total suite: 59 tests (was 6). `menu.py` and `plots.py` remain untested
  (interactive/visual — out of scope for this pass; see
  `PLAN_paper_figures.md`'s intermediate-task section for the full
  rationale and a separately-scoped follow-on task to bundle a real NORAD
  dataset for integration testing and distribution).

## [0.1.12] - 2026-07-16

### Removed

- **"Extend match pair" menu item** (Sequence operations) — superseded by
  "Motif extensions", which runs the same engine and adds an optional
  approximate/mutation-tolerant seed search that "Extend match pair" never
  had. Running "Motif extensions" with search method 1 (exact) reproduces
  "Extend match pair" exactly. `extend_match_input()` deleted; remaining
  Sequence-ops items renumbered 8–16 → 7–15.
- **Dead code**: `scramble_input()`, `markov_input()`, and the
  `scramble_kmer_pvalues()` algorithm it wrapped — unreachable from any menu,
  not called by the batch CLI or any test, and superseded by
  `markov_kmer_input()` (configurable Markov order, no shuffling required).
- **"Decompose motif"** (`decompose_motif_input()`, `decompose_motif()`) —
  was fully implemented and documented in the README but never wired into
  any menu, making it permanently unreachable. Removed rather than wired in;
  the README section is preserved verbatim in local planning notes in case
  this is revisited later.
- **Unused pip dependencies** `networkx` and `Levenshtein` — confirmed via a
  full-repo search to be imported nowhere in the codebase.

### Fixed

- **Documentation accuracy audit** (`README.md`, `docs/user_guide.md`,
  `llms.txt`) — corrected stale menu-numbering cross-references left over
  from the removals above; fixed several incorrect function
  signatures/defaults in `llms.txt` (`cores()`'s `minclen` default,
  `score_exact_motifs`'s `min_xmotif_type_support` default, missing
  keyword-only params on `rank_core_candidates`/`mutation_family_tests`,
  wrong parameter names on `find_gapped_motif_hits`/`score_gapped_motif`);
  fixed a broken `MarkovBackground` usage example (`.from_sequence()` /
  `.dinuc_prob` don't exist — it's a plain frozen dataclass constructed
  directly, with a `.probability()` method); removed the now-deleted
  `scramble_kmer_pvalues` from `llms.txt`'s import example and API list.
  Added previously-undocumented prompts (`min_occ`, hairpin overlay) to the
  Core neighbors prompt tables, and documented the Kaleido export-timeout /
  HTML-fallback behaviour (added in 0.1.9) for the first time.
- **Windows encoding crash risk** — ~17 `open()` calls across `io.py`,
  `menu.py`, `algorithms.py`, and `plots.py` had no explicit `encoding=`,
  so on Windows without UTF-8 mode opted in they'd use the locale codepage
  instead of UTF-8. Concretely confirmed broken: the "Open User Guide" menu
  item reads/writes `docs/user_guide.md`/`.html`, which contain 23 distinct
  non-ASCII characters (arrows, Greek letters, math symbols) — this would
  raise `UnicodeDecodeError`/`UnicodeEncodeError` on a typical Windows
  machine. All affected calls now pass `encoding='utf-8'` explicitly.
- **`open_file_with_default_software()` (io.py) could crash the caller** —
  unlike its sibling `open_pdf()`, it had no try/except around the
  macOS/Linux `subprocess.run(..., check=True)` calls, so a missing
  `xdg-open`/`open` binary or a headless/SSH/Docker environment would raise
  uncaught. Now self-guarded, matching `open_pdf()`'s pattern; failures
  print a message instead of propagating. This also fixes
  `init_summary()`'s CSV-write retry loop, which only caught
  `PermissionError` and could leak an unhandled exception from the
  now-fixed auto-open call.
- **`load_prefs()` could crash at startup, before any error handling was in
  place** — it ran before `menus()`'s try/except loop began, and called
  `os.makedirs()` (via `_prefs_path()`) with no guard. A restrictive/
  read-only config directory would crash the whole program rather than
  falling back to default preferences. Now guarded; verified the fix cannot
  cause stale preferences to be used, since `os.makedirs(..., exist_ok=True)`
  only fails when the directory doesn't exist and can't be created — meaning
  there was never a pre-existing readable `prefs.json` to lose in that path.
- **`pywin32` added as a Windows-only dependency** (`pywin32; sys_platform
  == 'win32'`) — the optional "auto-close a locked Excel workbook before
  overwriting" convenience (`_close_in_excel()`) previously never activated
  unless a user happened to have `pywin32` installed for an unrelated
  reason. Verified the platform marker is correctly skipped on macOS/Linux
  via a live `pip install --dry-run` dependency resolution.

### Added

- **`docs/build_docs.py`** — regenerates `docs/user_guide.html` and
  `docs/user_guide.pdf` from `docs/user_guide.md`. The PDF previously had no
  automated build path at all (one prior manual "updated the pdf" commit in
  the entire git history) and had gone stale. Uses `weasyprint -p
  --optimize-images`; the `-p`/`--presentational-hints` flag is required —
  without it, legacy HTML attributes like `<img width="200">` are ignored
  and images render at full native size.
- **`rna_lexis.menu.render_user_guide_html()`** — the HTML-rendering logic
  previously inlined in the "Open User Guide" menu handler is now a shared,
  reusable function, used by both that menu item and `build_docs.py`, so
  there is exactly one HTML template instead of two that could drift apart.
- **Troubleshooting sections** in both `README.md` (install/environment
  issues: kaleido/plotly version mismatches, missing `tkinter`/`pdf2svg`/
  `xdg-utils`, the benign multi-matplotlib `Axes3D` warning, legacy Windows
  terminal ANSI rendering) and `docs/user_guide.md` (usage issues: the RNA
  `u`→`t` query-mismatch gotcha, sparse-looking neighbour plots, the generic
  error-recovery message, sessions not auto-loading, minimum-occurrence
  requirements, the xmotif-length auto-expansion loop).
- **"System (non-pip) requirements" section in `README.md`** — documents
  `tkinter`, `pdf2svg`, and `xdg-utils` as external, non-pip-installable
  requirements for specific features, with install commands per platform.
- **Cross-links between `README.md`, `docs/user_guide.md`, and `llms.txt`**
  — each now points to the other two; previously none of them referenced
  each other, so a reader landing on any one had no way to discover the
  others.
- **`CLAUDE.md`** — repo-specific operational notes for future agent
  sessions (test command and coverage gaps, documentation-regeneration
  workflow, module layout, dependency-pin rationale, established
  conventions).

## [0.1.11] - 2026-06-22

### Fixed

- **Kaleido dependency reverted to 0.2.x to eliminate Chrome requirement** —
  version 0.1.9 widened the pin to `kaleido>=1`, which requires a separate
  Chrome installation (`kaleido_get_chrome`) and also requires Plotly ≥ 6 for
  `fig.write_image` to work.  Users with Plotly 5.x got version-mismatch
  warnings and broken static image export.  Reverted to `kaleido>=0.2,<1`,
  which bundles its own Chromium headless and works on Windows, macOS, and
  Linux without any extra browser install.  The 45-second export timeout added
  in 0.1.9 (which was the actual fix for the Windows hang) is retained.

- **Removed spurious `datetime` pip dependency** — `datetime` is a Python
  standard-library module and should not be listed as a pip dependency.
  Removed from `pyproject.toml`.

- **Removed duplicate `matplotlib` entry** — `matplotlib` was listed twice in
  `pyproject.toml`; the duplicate has been removed.

## [0.1.10] - 2026-06-22

### Fixed

- **`min_occ` filter in `plot_seq_nbrs` produced wrong results** — the previous
  implementation counted (s0_position, seq_position) pairs, so a sequence appearing
  once near three `s0` hits scored 3 (false positive) while a sequence appearing
  three times but near only one `s0` hit scored 1 (false negative).  Replaced with
  a simple global occurrence count: `len(find_all_matches(seq, txt))`.

### Added

- **`min_occ` parameter added to `plot_nbrs_condensed`** — the condensed neighbour
  plot now supports the same minimum-occurrence filter as the detailed plot (default 2).

- **`min_occ` exposed in the menu** — both "Core neighbors (detailed)" and "Core
  neighbors (condensed)" now prompt for a minimum-occurrence threshold (default 2,
  enter blank to keep default).

## [0.1.9] - 2026-06-16

### Fixed

- **Export timeout prevents hung Kaleido processes on Windows** — image export now
  runs in a `ThreadPoolExecutor` with a 45-second timeout.  On timeout, an
  interactive HTML fallback is saved automatically and a plain-text error report is
  written next to it.  Widened pins to `plotly>=5, kaleido>=1` (later corrected in
  0.1.11 — the Chrome requirement was an unintended side-effect).

## [0.1.8] - 2026-06-16

### Added

- **`min_occ` parameter for `plot_seq_nbrs`** — neighbours with fewer than
  `min_occ` global occurrences in the source sequence are omitted from the detailed
  neighbours plot.  Mutations are always shown regardless of this threshold.
  Default is 2 (unchanged behaviour).

### Changed

- **Neighbour plot legibility improvements** — annotation colours darkened to
  forest green (overlapping neighbours), dark purple (non-overlapping), dark blue
  (s0), and crimson (mutations); label font increased from 14 to 16 px; legend font
  increased to 13 pt with `entrywidth=260 px` to prevent column overlap.

## [0.1.7] - 2026-06-09

### Fixed

- **Crash on loading sessions saved by older versions** — `workdir = data['dir']`
  now uses `.get('dir')` with a fallback to the directory of the session file,
  so sessions that pre-date the `dir` key load without error.

- **Crash when pasting an empty or non-alphabetic sequence** — `load_from_paste`
  previously raised `ValueError` (uncaught) when the user submitted an empty paste,
  pasted text with no sequence characters, or cancelled the save-directory dialog.
  All three cases now print a message and return `None` (menu continues normally).

- **Any other unexpected exception in a menu handler no longer crashes the program** —
  a broad `except Exception` fallback in the main menu loop prints the error and
  returns to the main menu instead of propagating to the OS.

## [0.1.6] - 2026-06-09

### Fixed

- **"Open Core file" no longer regenerates the CSV on every click** — a SHA-256
  hash of the current inputs (xmotifs, corelist, sequence, and scoring parameters)
  is now stored in a `.chk` sidecar file alongside the CSV.  If the file exists,
  has the required statistical columns, and the inputs have not changed, it is opened
  directly without regeneration.  Regeneration still happens automatically when
  inputs change.

- **"Open Core file" no longer opens a second LibreOffice window on Linux** — before
  calling `xdg-open`, RNA_lexis now checks for a LibreOffice lock file
  (`.~lock.<filename>#`).  If present, `wmctrl` or `xdotool` is tried first to
  raise the existing window; no new `xdg-open` call is made regardless.

- **"Alignment score for two sequences" prompts now support blank-to-cancel** — the
  three numeric input prompts ("Enter start position of first/second sequence",
  "Enter sequence length") now accept a bare Enter as a cancel signal and return
  immediately to the menu.

- **"Alignment score for two sequences" default changed from Local to Global** —
  Smith–Waterman local alignment can return a shorter region than the window the
  user specified, which was surprising.  The default is now Needleman–Wunsch global
  alignment; local alignment is still available as option 2.

- **Local alignment now reports the best-matching region's position** — `AlignmentResult`
  gains four optional fields (`start_a`, `end_a`, `start_b`, `end_b`) populated by
  `gotoh_local`; `print_alignment` shows `"Best local region: seq1[X:Y] (N bp),
  seq2[Z:W] (M bp)"` when these are set.

- **Normalised score and E-value corrected for local alignment** — the self-alignment
  denominator and E-value query length now use the actual aligned region length rather
  than the user-requested window length.

## [0.1.5] - 2026-05-26

### Fixed

- **`cores()` missed short cores when `minxmlen` > `mincorelen`** — the function
  previously iterated only over k values that coincided with actual xmotif lengths,
  skipping intermediate lengths where valid shared cores could exist.  A core of
  length L was silently missed whenever no xmotif had exactly length L (e.g.
  `TGTATATA` at length 8 was found with `minxmlen=7` but lost with `minxmlen=10`,
  even though the same three xmotifs containing it were present in both cases).
  The fix iterates every integer k from the starting length down to `minclen`,
  letting the pool update naturally at each xmotif-length boundary.

## [0.1.3] - 2026-05-21

### Added

- **Mutation support for spacing / periodicity test** — after entering the motif,
  a second prompt asks for a mutation rate (0 = exact only; enter N for 1-per-N-nt,
  maximum N = 6).  When rate > 0, `find_with_mutations` is used and exact +
  approximate positions are merged (deduplicated) before the gap-cluster and
  Rayleigh tests.  The output header shows the match mode and exact/approx counts.

### Changed

- **Sequence operations menu reordered** into three visually distinct color blocks:
  - Items 1–9 (bold cyan) — single-sequence analysis: Find all matches, Search
    with mutations, Motif extensions, Print core, Motif spacing / periodicity test,
    Gapped motif search, Extend match pair, Covered area, Core neighbors (text export).
  - Items 10–13 (bold yellow) — statistical analysis: Rank core motifs (Markov/FDR),
    Mutation-family scoring, Alignment score for two sequences, K-mer Markov analysis.
  - Item 14 (dim white) — other: Export hairpins to CSV.
  A blank line separates each block.

- **`gen_menu` / `show_menu`** accept a new `splits` parameter (list of
  `(boundary_index, color)` pairs) for multi-block menu coloring; existing `split`
  (single int) behavior is unchanged.

## [0.1.2] - 2026-05-21

### Added

- **Motif spacing / periodicity test** (Sequence operations → 6).  Given a motif
  with m ≥ 3 exact occurrences, two statistical tests determine whether the copies
  are spaced more regularly than random placement would produce:

  **Gap cluster test** (primary) — estimates the candidate period T as the median
  consecutive gap, then counts how many of the m−1 gaps fall within δ = max(5, 5%·T)
  of T.  The count is tested against Bin(m−1, 2δ/L) with Bonferroni correction
  for the fact that T is estimated from the same data.  This test is sensitive to
  partial periodicity — a few tightly-clustered gaps among otherwise irregular
  spacing — and handles tandem copies (tiny gaps) and skipped periods (gaps ≈ 2T)
  correctly.

  **Rayleigh test** (confirmatory) — maps all m positions modulo T to angles on a
  circle and tests whether they cluster together (Mardia & Jupp, 2000 approximation).
  The mean resultant length R ∈ [0,1] measures concentration; Z = m·R².

  Either test reaching p < 0.05 gives a positive verdict.  When gap sizes are
  mixed, the output identifies small gaps (< 15%·T) as likely tandem copies and
  large gaps (> 175%·T) as likely skipped periods.

  A warning is shown when m < 6 (low power).  m = 2 reports the single gap only.

### Changed

- **Spacing test output** — gap cluster result is now shown first (it is the
  primary test); Rayleigh follows as confirmation.  The obsolete CV test
  (Monte Carlo coefficient-of-variation) has been removed: it was not sensitive
  to the common biological pattern of a few near-equal gaps among otherwise
  variable spacing.

- **`spacing_periodicity_test` API** — the function no longer accepts `n_sim`
  and no longer returns the CV-related keys (`cv_obs`, `cv_null_median`,
  `cv_null_p5`, `p_cv`, `n_sim`).  The `numpy` import is no longer required in
  `statistical.py`.

## [0.1.1] - 2026-05-21

### Added

- **"Clear workspace" option in the main menu (option 8).**  A new menu entry lets
  the user delete all generated files from the working directory in one step, while
  preserving the loaded input file and any other files with a FASTA extension
  (`.fasta`, `.fa`, `.fna`, `.fas`).

  Before deleting anything the menu shows two lists — files that **will be deleted**
  and files that **will be kept** — and requires the user to type `YES` (exact,
  case-sensitive) to confirm.  Typing anything else, or pressing Enter, cancels the
  operation without touching any files.

  If the user confirms, all non-preserved files are removed and the session is reset
  so that the tool returns to the *Choose Input* screen, ready to reload the original
  sequence fresh.  If no deletable files exist the menu reports this and returns
  without prompting.

  The check for which input file to keep is based on the **actual file path** of the
  loaded sequence rather than its extension alone, so plain-text and other non-FASTA
  source files are protected correctly.

- **"Load new input" is now option 9 and "Quit" is option 10** (shifted by the
  insertion of "Clear workspace" at position 8).

## [0.1.0] - 2026-05-20

### Added

- **Statistical layer in the initialisation summary.**  When a sequence is loaded
  (or the summary is regenerated), `init_summary` now scores every xmotif and
  core against a transcript-specific first-order Markov background and appends
  the following columns to `<session>_test_init.csv`:

  | Column | Description |
  |---|---|
  | `nonoverlap_count` | Non-overlapping occurrence count |
  | `coverage_bp` | Base-pairs covered by non-overlapping occurrences |
  | `area_score` | `length × nonoverlap_count` coverage proxy |
  | `xmotif_type_support` | Number of distinct xmotif families containing this core |
  | `inside_xmotif_count` / `outside_xmotif_count` | Occurrences inside vs. outside xmotifs |
  | `expected_markov` | Expected count under a transcript-specific Markov model |
  | `enrichment_markov` | `observed / expected` enrichment ratio |
  | `p_markov` | Poisson upper-tail p-value for motif enrichment |
  | `q_markov` | Benjamini–Hochberg FDR-adjusted p-value across all tested motifs |
  | `statistically_supported` | Boolean: `q_markov` below threshold **and** enriched |
  | `core_class` | `'supported'`, `'enriched_only'`, or `'not_supported'` |
  | `rank_statistical` | Rank by statistical support (primary) |
  | `rank_coverage` | Rank by coverage (secondary) |

  The CSV is now **sorted by statistical support first** (`statistically_supported`
  descending, `q_markov` ascending), falling back to `p_stable` for non-RNA/DNA
  sessions or sequences where the statistical layer is not applicable.
  Statistical columns are only added when the loaded sequence alphabet is a
  subset of `{a, c, g, t, u}`.

- **Output filename changed: `_init.csv` → `_test_init.csv`.**  The initialisation
  summary is now written to `<session>_test_init.csv`.  Old `_init.csv` files are
  not read or overwritten; they can be safely kept or deleted.  When an existing
  `_test_init.csv` is found and already contains all four statistical columns
  (`expected_markov`, `enrichment_markov`, `p_markov`, `q_markov`), it is opened
  directly without regenerating.  If statistical columns are absent (e.g. a file
  written by an older run), the summary is regenerated automatically.

- **`init_summary` `force` parameter.**  `init_summary(fn, xm, cores, txt,
  force=False)` accepts a new keyword argument.  When `force=True` the summary
  is always regenerated, bypassing the caching check.

- **Three new Sequence operations menu items:**

  - **Rank core motifs (Markov/FDR)** — Interactive wrapper around
    `rank_core_candidates`.  Enumerates all shared substrings of xmotifs within a
    configurable length range, scores each candidate against the Markov background,
    and saves a ranked CSV (`<session>_ranked_cores_markov.csv`).  Prompts: minimum
    and maximum candidate length, minimum xmotif-type support, minimum enrichment,
    FDR threshold, and output file path.

  - **Mutation-family scoring** — Interactive wrapper around
    `mutation_family_tests` / `best_mutation_family_per_motif`.  Tests one or more
    motifs (typed, from the current cores, or from the current xmotifs) at every
    Hamming radius allowed by the mutation cap, scores each family against the
    Markov background, and reports the best-supported radius per motif.  Output:
    `<session>_mutation_family_tests.csv` (all radii) and
    `<session>_mutation_family_tests_best.csv` (one row per motif).

  - **Gapped motif search** — Interactive wrapper around `score_gapped_motif` /
    `find_gapped_motif_hits`.  Searches for all occurrences of a
    `LEFT[gap:min–max]RIGHT` anchor pattern, scores the whole family under the
    Markov background, and saves a CSV with per-hit positions and the family-level
    enrichment statistics (`<session>_gapped_motif.csv`).

- **`rna_lexis.statistical` module** (`statistical.py`).  New pure-Python module
  providing all statistical primitives:
  - `MarkovBackground` — frozen dataclass computing first-order (or higher)
    Markov transition probabilities with pseudocount smoothing.
  - `score_exact_motifs` — scores a list of exact motifs against the Markov
    background; returns enrichment, p-values, FDR q-values, and classification.
  - `rank_core_candidates` — enumerates shared sub-strings of xmotifs, scores
    them, and returns a ranked table.
  - `mutation_family_tests` — tests Hamming-radius families across all allowed
    radii and applies Markov/FDR filtering.
  - `best_mutation_family_per_motif` — selects the strongest accepted radius per
    motif from the full family-test table.
  - `find_gapped_motif_hits` — finds all `LEFT[gap]RIGHT` pattern hits with
    per-hit gap lengths and positions.
  - `score_gapped_motif` — scores a gapped pattern as a single family under the
    Markov background.
  - `write_rows_csv` — convenience helper to write a list of dicts to a CSV.

- **`rna_lexis_stat_cli` console script** (`test_cli.py`).  A new command-line
  interface for batch statistical workflows:

  ```
  rna_lexis_stat_cli score-exact     --fasta FILE --motifs m1 m2 …
  rna_lexis_stat_cli rank-cores      --fasta FILE
  rna_lexis_stat_cli mutation-families --fasta FILE --motifs m1 m2 …
  rna_lexis_stat_cli gapped-motif    --fasta FILE --left LEFT --right RIGHT
  ```

- **`rna_lexis_stat` console script** — alias for `rna_lexis` (same interactive
  menu).

### Changed

- **Sequence operations submenu re-numbered.**  The three new statistical items
  are inserted after *Print core* (item 4) as items 5, 6, and 7.  The existing
  items 5–11 shift to 8–14.

- **`__version__` set to `"0.1.0"`** in `rna_lexis/__init__.py`.

## [0.0.10] - 2026-05-03

### Changed

- **Cores are now maximally extended before being reported.**  After the raw
  core set is identified, each core is grown left and right in the source text
  as long as every occurrence shares the same flanking character.  The extended
  string is strictly more informative than the original: it retains the same
  occurrence count and removes the ambiguity implied by the shorter form.
  Multiple raw cores that collapse to the same maximal string are deduplicated,
  eliminating the sliding-window redundancy that previously inflated core counts
  in repeat-rich sequences (e.g. EBYT lyrics: 87 raw cores → 5 maximal cores).

### Added

- **`extend_core_maximally(txt, core)`** (`algorithms`).  Public helper that
  implements the bidirectional greedy extension described above.  Accepts the
  full source text and a core string; returns the maximally extended string.

## [0.0.9] - 2026-05-03

### Changed

- **K-mer scramble analysis replaced by K-mer Markov analysis.**  The
  shuffle-based permutation test (*K-mer scramble analysis*) has been replaced
  in the Sequence operations menu by an analytical Markov test (*K-mer Markov
  analysis*) that requires no shuffling and runs instantly for any k.

  The key improvement is a new `order` parameter in `markov_kmer_pvalues`
  (default `order=1`, dinucleotide null) that controls how much context the
  null model conditions on.  The previous implementation always used a
  (k−1)-th order model, which is degenerate for long k-mers in short
  sequences (expected ≈ observed → all p-values ≈ 0.5).  With `order=1` the
  dinucleotide frequencies serve as the null, giving meaningful expected counts
  and well-calibrated p-values even for 8-mers or longer in transcripts of a
  few kilobases.

  **New menu prompts:**

  | Prompt | Default |
  |---|:---:|
  | k-mer length | 6 |
  | Markov background order (0 = nucleotide, 1 = dinucleotide, …) | 1 |
  | Output CSV file | `<name>_kmer<k>_order<order>_markov_pvalues.csv` |

  **Output CSV columns** are the same as before except `exceed_count` and
  `below_count` are replaced by `expected_count` (the analytical Markov
  expectation).

  `scramble_kmer_pvalues` and `scramble_input` remain available in the API
  for users who need the permutation-based null.

### Fixed

- `_markov_expected_count_order` added alongside the existing
  `_markov_expected_count`; the new helper supports arbitrary Markov order
  (0 through k−2) via the generalised Prum/Schbath formula.

## [0.0.8] - 2026-05-02

### Added

- **Core neighbors (condensed) plot** (`plot_nbrs_condensed`).  A new *Core
  neighbors (condensed)* option in the Plots menu shows the same neighbourhood
  data as the detailed plot in a compact three-band layout designed for long
  sequences or large neighbour sets.  s0 occurrences are drawn as light-blue
  rectangles at y = 0; all neighbours are collapsed into a single semi-
  transparent **density strip** at y = 1 (green = overlapping with s0, magenta
  = non-overlapping; opacity scales with the number of sequences sharing each
  position); single-nucleotide mutations of s0 are collapsed into a red
  **mutation strip** at y = −1.  Hairpin regions from a companion
  `*_hairpins.csv` file are shown as orange bands below s0.  Accepts the same
  prompts as the detailed plot (sequence, window width, cores/xmotifs,
  title, output file, x-axis range).

- **Core neighbors (text export)** (`export_nbrs_condensed`).  A new *Core
  neighbors (text export)* option in the Sequence operations menu exports the
  sequence regions associated with a query core and its neighbourhood as a CSV
  file.  Each row is one merged region — a contiguous stretch of the full text
  where any combination of s0, its single-nucleotide mutations, or its
  neighbours clusters together.  Regions are built by expanding every
  occurrence of every tracked sequence by ±wd nucleotides and merging
  overlapping intervals, so each row captures the full extent of a cluster
  (including neighbours that bridge gaps between s0 occurrences).  Output
  columns: `start`, `end`, `seq` (= `txt[start:end+1]`).  Results are printed
  to the terminal and the CSV is opened automatically on completion.  Default
  output filename: `<session>_<seq>_regions.csv` in the session directory.

## [0.0.7] - 2026-04-29

### Added

- **Hierarchical motif decomposition (`decompose_motif`).**  A new *Decompose motif*
  option (Sequence operations → 10) tests whether a motif's enrichment is genuine or
  fully explained by shorter sub-sequences.  For every contiguous sub-k-mer of the
  input motif at each length from `len(motif)` down to a user-configurable `min_k`
  (default 4), the expected count is computed analytically from a (k−1)-th order
  Markov model using the Prum/Schbath formula — no shuffling required.  A Poisson
  exact p-value then answers "is this k-mer count surprising given the shorter
  context?"  A single BH-FDR-adjusted p-value (`pvalue_bh`) is reported across all
  sub-k-mers tested, and the shortest level reaching significance is printed to the
  terminal.  Output columns: `level`, `kmer`, `real_count`, `expected_count`,
  `pvalue_over`, `pvalue_under`, `pvalue_bh`, `direction`, `significant`.

- **Analytical Markov k-mer test (`markov_kmer_pvalues`).**  A fast analytical
  counterpart to the shuffle-based scramble test.  For each k-mer the expected count
  under a configurable-order Markov model is computed from the Prum/Schbath formula
  and tested with a Poisson exact p-value, with BH-FDR correction.  Runs instantly
  for any k.  Available in the API (`algorithms.markov_kmer_pvalues`); exposed as
  the *K-mer Markov analysis* menu item as of v0.0.9.

## [0.0.6] - 2026-04-29

### Changed

- **Neighbours plot default sequence set changed to xmotifs.**  In the
  *Visualise neighbourhood* prompt, pressing Enter (or any input other
  than `1`) now selects xmotifs instead of cores.  Enter `1` explicitly
  to use cores.

## [0.0.5] - 2026-04-29

### Added

- **`p_stable` column in `_init.csv`.**  Each xmotif and core is now assigned a
  formal one-sided Binomial p-value measuring sequence stability:

  ```
  p_stable = P(X ≤ nmut)  where  X ~ Bin(L × N_total, mutr)
  ```

  `L` is the sequence length, `N_total = count + n_approx` is the total number
  of occurrences (exact + approximate), `nmut` is the observed total Hamming
  distance across all approximate matches, and `mutr` is the null mutation rate
  (default 1/6 — the same tolerance used during the search).  A small `p_stable`
  means the sequence appears with **fewer** mutations than expected under the null
  rate, i.e., it is **conserved / stable**.  The CSV is now sorted by `p_stable`
  ascending.

- **`n_approx` column in `_init.csv`** — count of approximate (non-exact)
  occurrences, so that `count + n_approx` gives the total number of occurrences
  used in the Binomial test.

### Changed

- **`pmut` in `_init.csv` is now the clean observed mutation rate per site**
  (`nmut / (L × N_total)`).  The previous formula contained a double-counting
  bug (`N + len(mutations)` counted exact matches twice) and used smoothing
  constants that obscured the value.  `pmut` is retained as a diagnostic column
  but `p_stable` is the recommended ranking criterion.

- **K-mer scramble analysis now reports two-sided tests.**  The output CSV
  previously had a single `pvalue` column (upper-tail only, flagging
  over-representation).  It now has two one-sided p-values:
  - `pvalue_over` — fraction of shuffles where the k-mer count was ≥ the
    observed count; a small value indicates the k-mer appears **more** often
    than expected by chance.
  - `pvalue_under` — fraction of shuffles where the k-mer count was ≤ the
    observed count; a small value indicates the k-mer appears **less** often
    than expected by chance.
  - K-mers with no enrichment or depletion will have both p-values near **0.5**.
  - Each direction gets its own E-value (`evalue_over`, `evalue_under`) and
    BH-adjusted p-value (`pvalue_over_bh`, `pvalue_under_bh`).  BH correction
    is applied separately within each family of m tests (where m is the number
    of distinct k-mers observed in the sequence).
  - A `direction` column (`'over'` / `'under'`) records which effect is
    stronger for each k-mer.
  - Rows are sorted by `min(pvalue_over, pvalue_under)` so the most extreme
    k-mers in either direction appear first.
  - The raw counts `exceed_count` and `below_count` are both reported so
    results can be recomputed with a different N without re-running.
  - **Note:** this is a breaking change for any downstream code that parsed
    the old single-column `pvalue` output.

## [0.0.4] - 2026-04-16

### Fixed

- **Cancel in file dialog no longer crashes the program.**  The subprocess-based
  tkinter dialog now writes its result prefixed with a sentinel string
  (`__DIALOG_RESULT__:`), so any incidental library output that ends up on
  stdout (GTK warnings, libpng messages, etc.) is ignored instead of being
  mistaken for a file path.  `choose_file` also guards against any non-empty
  but invalid path returned by `openFile` with an `os.path.isfile` check.
- **Corrupt or missing session JSON no longer crashes the program.**  When the
  auto-detected session file cannot be loaded (`load_session` returns `None`),
  the program now clears the stale reference and falls back to the input-source
  menu instead of crashing with `TypeError`.

## [0.0.3] - 2026-04-16

### Added

- **Persistent preferences** (`io.load_prefs` / `io.save_prefs`): the last used
  directory and a user-defined default data directory are now stored in a
  platform-appropriate config file across sessions:
  - Linux: `~/.config/RNA_lexis/prefs.json`
  - macOS: `~/Library/Application Support/RNA_lexis/prefs.json`
  - Windows: `%APPDATA%\RNA_lexis\prefs.json`
- On startup the tool seeds the initial directory from `last_used_dir`, falling
  back to `default_data_dir`, then the current working directory.
- `last_used_dir` is updated and saved automatically every time a file is opened.
- `default_data_dir` can be set (or cleared) from the *Change setting* menu.  It
  is also shown in *Show settings*.
- New dependency: `platformdirs` (cross-platform config-directory resolution).

## [0.0.2] - 2026-04-16

### Fixed

- **FASTA input:** when loading a `.fa` / `.fasta` file the header line (starting
  with `>`) is now used only to detect the format and extract the gene name; it is
  no longer included in `txt` or `txtb`.  Previously, alphabetic tokens from the
  header (e.g. `RC`, `PastedSequence`, `nt`) were prepended to the sequence,
  corrupting all downstream analysis.  The gene name extracted from the header is
  stored in the session dict under the key `gene_name`.  Multi-record FASTA files
  are handled correctly: all `>`-prefixed lines are stripped before the sequence
  records are concatenated.

## [0.0.1] - 2025-11-01

### Added

- Initial release.

# Changelog

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

# Changelog

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

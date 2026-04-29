# Changelog

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

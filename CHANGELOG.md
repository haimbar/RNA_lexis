# Changelog

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

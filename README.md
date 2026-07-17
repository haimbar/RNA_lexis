# RNA_lexis

**RNA_lexis** is an interactive command-line tool for exploring RNA (or DNA) sequence data.

Starting from a nucleotide sequence, it automatically identifies repeating subsequences
called **xmotifs** and extracts their shorter conserved cores. From there you can:

- Visualise the neighbourhood of any core sequence
- Study k-mer frequency distributions (z-score, rank-frequency, histograms)
- Test k-mer enrichment analytically with Markov-model p-values (no shuffling)
- Generate sequence logos
- Inspect coverage across the full transcript
- Align pairs of subsequences (Gotoh global / Smith–Waterman local)
- Search for exact and approximate (mutation-tolerant) matches
- Extend matching pairs greedily to find the longest shared context
- **Plot self-similarity as a semicircular arc diagram** — pairwise Hamming-bounded extensions among every occurrence of a seed motif, colored by tandem-repeat hop distance
- **Plot a shared-motif diagram against a second sequence** — connects every exact-match motif shared between the loaded sequence and a comparison sequence (paste, file, Ensembl ENST, genomic coordinates, or an ENCODE cCRE accession) as a two-row arc diagram
- **Compare motif-coverage sliding-window curves** — against another transcript (e.g. a cross-species ortholog, with automatic length-mismatch warnings), or against the same sequence's own GC-content or homopolymer/low-complexity-run density
- Fetch sequences directly from the Ensembl REST API by transcript ID
- **Score motifs against a transcript-specific Markov background with FDR correction**
- **Test mutation-tolerant motif families for statistical enrichment**
- **Search for anchor-gap-anchor gapped motif patterns**
- **Test whether motif occurrences are spaced periodically (gap-cluster and Rayleigh tests)**
- Clear workspace (remove all generated files while preserving the input sequence)
- **Try it instantly on a real dataset** — human and mouse NORAD transcripts
  ship with the package itself (`pip install` and go, no fetching or
  pasting required); see *Load example dataset* below.

## Documentation

This README covers installation, features, and troubleshooting. Two other
docs go deeper:

- **[docs/user_guide.md](docs/user_guide.md)** — full menu-by-menu walkthrough
  of the interactive tool, with every prompt, default, and output format
  documented (also available as `docs/user_guide.html`/`.pdf`, or via the
  app's own *Open User Guide* menu item).
- **[llms.txt](llms.txt)** — API reference for using RNA_lexis as a Python
  library (functions, data structures, worked examples), written for both
  developers and AI coding agents.

## Initialisation summary (`_test_init.csv`)

When a sequence is loaded for the first time, an initialisation CSV is written that
characterises every xmotif and core. The file is named `<session>_test_init.csv`.

### Core columns (all sequences)

| column | description |
|---|---|
| `xm` | 1 if the sequence is an xmotif, 0 otherwise |
| `core` | 1 if the sequence is a core (a sequence can be both) |
| `count` | number of exact occurrences in the transcript |
| `n_approx` | number of approximate occurrences (Hamming distance ≤ `maxmut`) |
| `cover` | `len × count` — nucleotides covered by exact matches |
| `numt` | total Hamming distance across all approximate matches |
| `pmut` | observed mutation rate per site: `numt / (len × (count + n_approx))` |
| `p_stable` | Binomial p-value for sequence stability (see below) |
| `maxmut` | maximum mutations allowed by the search (`floor(len × mutr)`, capped at `M`) |

### Statistical columns (RNA/DNA sequences only)

For sequences whose alphabet is a subset of `{a, c, g, t, u}`, the following
columns are also added by scoring each motif against a transcript-specific
first-order Markov background:

| column | description |
|---|---|
| `expected_markov` | Expected count under the Markov null |
| `enrichment_markov` | `observed / expected` enrichment ratio |
| `p_markov` | Poisson upper-tail p-value for motif enrichment |
| `q_markov` | Benjamini–Hochberg FDR-adjusted p-value across all tested motifs |
| `statistically_supported` | `True` if `q_markov` is below threshold and the motif is enriched |
| `core_class` | `'supported'`, `'enriched_only'`, or `'not_supported'` |
| `rank_statistical` | Primary rank (by statistical support) |
| `rank_coverage` | Secondary rank (by coverage) |
| `coverage_bp` | Base-pairs covered by non-overlapping occurrences |
| `xmotif_type_support` | Number of distinct xmotif families containing this core |

When statistical columns are present the CSV is sorted by `statistically_supported`
(descending) then `q_markov` (ascending). Otherwise it is sorted by `p_stable`.

### Stability p-value (`p_stable`)

Each sequence is tested against the null hypothesis that every site mutates independently
at the background rate `mutr` (default 1/6 — the same tolerance used during the search):

```
p_stable = P(X ≤ numt)   where   X ~ Bin(len × N_total, mutr)
```

`N_total = count + n_approx` is the total number of occurrences used as evidence.
A small `p_stable` means the sequence appears **more conserved than expected by
chance**, i.e., it is stable.
A sequence with no strong conservation signal will have `p_stable ≈ 0.5`.

## Maximal core extension

Cores are automatically extended to their **maximal unique form** before being
reported.  Starting from a raw core, the algorithm grows the string one character
at a time in each direction for as long as every occurrence in the transcript
shares the same flanking character.  Multiple raw cores that collapse to the same
maximal string are deduplicated, eliminating the sliding-window redundancy common
in repeat-rich sequences.

## K-mer Markov analysis

The *K-mer Markov analysis* menu option tests every observed k-mer analytically —
no shuffling or random seeds required.  For each k-mer the expected count is
computed from the Prum/Schbath formula under a configurable-order Markov null and
tested with a Poisson exact p-value.  The output CSV reports **two one-sided
p-values** per k-mer:

| column | description |
|---|---|
| `expected_count` | Expected count under the Markov null |
| `pvalue_over` | P(X ≥ obs) — small = **over-represented** |
| `evalue_over` | `pvalue_over × m` (expected false positives among *m* k-mers tested) |
| `pvalue_over_bh` | BH-adjusted p-value for over-representation (FDR) |
| `pvalue_under` | P(X ≤ obs) — small = **under-represented** |
| `evalue_under` | `pvalue_under × m` |
| `pvalue_under_bh` | BH-adjusted p-value for under-representation (FDR) |
| `direction` | `'over'` or `'under'` — which effect is stronger |

The **Markov background order** controls the null: `order=0` conditions only on
nucleotide frequencies (equivalent to the shuffle null); `order=1` (default)
conditions on dinucleotide frequencies and is recommended for transcripts of a
few kilobases.  BH correction is applied separately for each direction across all
*m* k-mers.  Rows are sorted by `min(pvalue_over, pvalue_under)`.

## Statistical motif analysis

### Rank core motifs (Markov/FDR)

The *Rank core motifs (Markov/FDR)* option (Sequence operations → 9) enumerates
all shared substrings of the current xmotifs within a configurable length range,
scores each candidate against the transcript-specific Markov background, and saves
a ranked CSV. Candidates are ranked first by statistical support (`q_markov` below
threshold, enrichment above threshold), then by coverage.

### Mutation-family scoring

The *Mutation-family scoring* option (Sequence operations → 10) tests one or more
motifs at every Hamming radius allowed by the mutation cap. For each motif and
radius, the full neighbourhood of sequences within that distance is counted and
scored against the Markov background. The result shows whether the approximate
family as a whole is enriched beyond chance — going beyond simply allowing
mismatches. The best-supported radius per motif is saved to a `_best.csv` file.

### Gapped motif search

The *Gapped motif search* option (Sequence operations → 6) finds all occurrences
of a pattern of the form `LEFT[gap:min–max]RIGHT`: two exact anchor sequences
separated by a variable-length gap. The whole family is scored under the Markov
background, and individual hits (with exact positions and gap lengths) are saved
to a CSV.

### Motif spacing / periodicity test

The *Motif spacing / periodicity test* option (Sequence operations → 5) takes a
motif, finds all exact occurrences, and tests whether their spacing is more
periodic than expected by chance. Two complementary tests are applied:

**Gap cluster test** — counts how many consecutive gaps fall within a tolerance
window δ = max(5, 5 % × T) of the candidate period T (median gap). Tested against
a Binomial null Bin(n\_gaps, 2δ/L) with Bonferroni correction for estimating T
from the same data.

**Rayleigh test** — maps occurrence positions modulo T to angles and measures the
mean resultant length R (0 = random, 1 = perfectly periodic). A large R indicates
that positions cluster at a consistent phase relative to T.

Either test reaching p < 0.05 is reported as a positive periodicity signal.
Requires at least 3 exact occurrences.

### Batch CLI (`rna_lexis_stat_cli`)

A separate command-line interface supports batch statistical workflows without
entering the interactive menu:

```
rna_lexis_stat_cli score-exact       --fasta FILE --motifs m1 m2 …
rna_lexis_stat_cli rank-cores        --fasta FILE
rna_lexis_stat_cli mutation-families --fasta FILE --motifs m1 m2 …
rna_lexis_stat_cli gapped-motif      --fasta FILE --left LEFT --right RIGHT
```

## Requirements

- Python 3.10 or newer
- Windows 10/11, macOS, or Linux

### System (non-pip) requirements

A few features depend on tools that live outside the Python packaging
ecosystem, so `pip install` can't provide them:

- **`tkinter`** — required for the native file/folder picker dialogs. Bundled
  with the standard Python installer on Windows and macOS, but many minimal
  Linux installs (e.g. server/Docker base images) ship Python without it.
  Install separately if dialogs silently return "no selection":
  `sudo apt install python3-tk` (Debian/Ubuntu) or
  `sudo dnf install python3-tkinter` (Fedora).
- **`pdf2svg`** (optional) — only needed for SVG-format sequence logo output;
  without it, logo plots fall back to PDF automatically. Install via
  `apt install pdf2svg` / `brew install pdf2svg`, or download a Windows build
  manually.
- **`xdg-utils`** (Linux only, optional) — provides `xdg-open`, used to
  auto-open saved plots/CSVs in their default application. Present on
  virtually all Linux desktops by default; may be missing on minimal/headless
  servers or containers, in which case files are still saved normally, just
  not auto-opened.

## Installation

```
pip install .
```

## Usage

After installation, launch the interactive menu with:

```
rna_lexis
```

Or with the statistical-workflow alias:

```
rna_lexis_stat
```

Or, without installing:

```
python -m rna_lexis
```

## Example dataset

Two real datasets ship inside the package itself — no internet connection,
account, or file to prepare — so you can try RNA_lexis immediately after
installing:

- **NORAD (human)** — Ensembl transcript `ENST00000565493`, 5401 nt.
- **NORAD (mouse)** — Ensembl transcript `ENSMUST00000192863`, 4945 nt.

From the interactive menu, choose **Load example dataset (NORAD)** at the
input-source prompt. NORAD is a well-studied lncRNA with a tandem repeat of
the Pumilio Response Element (core motif `tgtatata`), so the default
discovery pipeline reliably surfaces a real, biologically meaningful
repeat structure — a good way to confirm your installation works and see
what the tool finds on real data. Programmatically — note that `read_text()` strips non-alphabetic
characters but doesn't know about FASTA headers, so strip the header line
first:

```python
from rna_lexis.io import example_dataset_path

path = example_dataset_path("NORAD_human")
with open(path, encoding="utf-8") as f:
    txt = "".join(line for line in f if not line.startswith(">"))
txt = "".join(c for c in txt.lower() if c.isalpha())
```

See `src/rna_lexis/data/README.md` in the package source for full
provenance (exact transcript IDs, assembly versions, fetch date).

## Troubleshooting

- **Plots fail to export as PNG/SVG, or you see a kaleido/plotly version-mismatch
  warning.** RNA_lexis pins `kaleido<1` deliberately (see `CHANGELOG.md` 0.1.11) —
  `kaleido>=1` requires a separate Chrome install and only works with Plotly ≥ 6.
  If your environment ends up with an incompatible kaleido anyway (e.g. another
  project in the same environment upgraded it), reinstall with
  `pip install "kaleido>=0.2,<1"`. Do **not** follow kaleido's own suggested fix
  of `pip install -U kaleido` — that installs the incompatible version.
- **A saved plot or CSV doesn't open automatically, or you see "Could not open
  file automatically."** This is expected on headless/SSH/Docker Linux systems
  without `xdg-open`, or if no default application is registered for the file
  type. The file is still saved correctly in both cases — only the auto-open
  convenience step failed.
- **The file/folder picker dialog does nothing — no window appears, the prompt
  just returns immediately.** Usually means `tkinter` isn't installed (common on
  minimal Linux images — see *System requirements* above), or you're in an SSH
  session without X11 forwarding / a display. Paste a file path directly, or
  add `-X` to your SSH command.
- **Logo plots say `Can't find the 'weblogo' executable`.** `pip install` puts
  the `weblogo` script in a directory that may not be on `PATH` — e.g.
  `~/.local/bin` on Linux (user installs) or `Scripts\` inside a Windows venv.
  Activate your venv, or add that directory to `PATH`.
- **SVG logo output comes out as PDF instead.** SVG requires the external
  `pdf2svg` binary (see *System requirements* above); without it, RNA_lexis
  silently falls back to PDF and prints a note explaining why.
- **A `UserWarning: Unable to import Axes3D... multiple versions of Matplotlib`
  appears on every plot.** Harmless — it means both a system package (e.g.
  `apt`'s `python3-matplotlib`) and a pip-installed `matplotlib` are present.
  Plotting still works; to silence it, remove one of the two installations.
- **Menu colors show up as garbled text like `[1;36m`, and the screen doesn't
  clear between menus.** Your terminal isn't interpreting ANSI/VT100 escape
  codes — typical of legacy Windows `cmd.exe`. Use Windows Terminal or
  PowerShell 7+ instead.
- **Can't overwrite the summary CSV on Windows — `PermissionError`.** The file
  is usually open in Excel. RNA_lexis will try to close it automatically via
  `pywin32` (installed automatically on Windows); if that fails, close the
  workbook manually and retry.

## Updating the documentation

After editing `docs/user_guide.md`, regenerate the checked-in HTML/PDF
snapshots before committing:

```
python docs/build_docs.py
```

(The HTML also regenerates itself automatically any time a user opens it
from the app's *Open User Guide* menu item — this script exists mainly to
keep the PDF, which has no built-in refresh path, from going stale.)

## Authors

Haim Bar (haim.bar@uconn.edu) and Assaf Bester (bestera@technion.ac.il)

## License

MIT

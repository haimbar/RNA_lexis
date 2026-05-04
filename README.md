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
- Fetch sequences directly from the Ensembl REST API by transcript ID

## Initialisation summary (`_init.csv`)

When a sequence is loaded for the first time, an initialisation CSV is written that
characterises every xmotif and core. Key columns:

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

### Stability p-value (`p_stable`)

Each sequence is tested against the null hypothesis that every site mutates independently
at the background rate `mutr` (default 1/6 — the same tolerance used during the search):

```
p_stable = P(X ≤ numt)   where   X ~ Bin(len × N_total, mutr)
```

`N_total = count + n_approx` is the total number of occurrences used as evidence.
A small `p_stable` means the sequence appears **more conserved than expected by
chance**, i.e., it is stable.  The CSV is sorted by `p_stable` ascending.
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

## Hierarchical motif decomposition

The *Decompose motif* menu option answers the question: **is a motif's enrichment
genuine, or is it fully explained by shorter sub-sequences?**

For every contiguous sub-k-mer of the input motif at each length from `len(motif)`
down to a configurable minimum (default 4), the expected count is derived analytically
from a (k−1)-th order Markov model using the Prum/Schbath formula:

| k | Expected count formula |
|---|---|
| 2 | `count(a) × count(b) / n` |
| ≥ 3 | `count(w[:-1]) × count(w[1:]) / count(w[1:-1])` |

A Poisson exact p-value then tests whether the observed count is surprising given the
shorter context.  BH-FDR correction is applied across all sub-k-mers tested.  The
output CSV includes `pvalue_over`, `pvalue_under`, `pvalue_bh`, `direction`, and
`significant`; the shortest level reaching significance is printed to the terminal.

## Requirements

- Python 3.10 or newer
- Windows 10/11, macOS, or Linux

## Installation

```
pip install .
```

## Usage

After installation, launch the interactive menu with:

```
rna_lexis
```

Or, without installing:

```
python -m rna_lexis
```

## Authors

Haim Bar (haim.bar@uconn.edu) and Assaf Bester (bestera@technion.ac.il)

## License

MIT

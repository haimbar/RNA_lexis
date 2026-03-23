# RNA_lexis

**RNA_lexis** is an interactive command-line tool for exploring RNA (or DNA) sequence data.

Starting from a nucleotide sequence, it automatically identifies repeating subsequences
called **xmotifs** and extracts their shorter conserved cores. From there you can:

- Visualise the neighbourhood of any core sequence
- Study k-mer frequency distributions (z-score, rank-frequency, histograms)
- Generate sequence logos
- Inspect coverage across the full transcript
- Align pairs of subsequences (Gotoh global / Smith–Waterman local)
- Search for exact and approximate (mutation-tolerant) matches
- Extend matching pairs greedily to find the longest shared context
- Fetch sequences directly from the Ensembl REST API by transcript ID

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

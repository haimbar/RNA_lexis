# Bundled example dataset — NORAD

Two small FASTA files (human and mouse Norad), shipped as part of the
`rna_lexis` package (see `pyproject.toml`'s `[tool.setuptools.package-data]`)
so a fresh `pip install` gives every user a ready-to-try, real dataset —
no fetching, pasting, or account required.

NORAD (**N**oncoding RNA **A**ctivated by **D**NA damage; human gene
`LINC00657`) is a well-studied long noncoding RNA with a distinctive,
easily-rediscoverable feature: roughly a dozen tandem repeats of the
Pumilio Response Element (PRE, consensus `UGUANAUA`), which is what makes
it a good "does the pipeline actually work" smoke-test sequence — the
repeat structure is real biology, not a synthetic pattern, and running
`find_boundary`/`cores` on it should surface PRE-like motifs among the
discovered xmotifs/cores.

## Files

| File | Source | Length |
|---|---|---|
| `NORAD_human.fasta` | Ensembl transcript `ENST00000565493` (NORAD-201), GRCh38, gene `LINC00657` / `ENSG00000260032` | 5401 nt |
| `NORAD_mouse.fasta` | Ensembl transcript `ENSMUST00000192863` (Norad-202), GRCm39, gene `Norad` / `ENSMUSG00000097125` | 4945 nt |

Both fetched 2026-07-16 via the Ensembl REST API (`cdna` sequence type —
note this returns a DNA-alphabet representation, i.e. `t` not `u`, per
standard cDNA convention, even though the gene is noncoding RNA).
`NORAD_human.fasta` was fetched with `rna_lexis.io.fetch_enst_cdna()` — the
same function the app's own "Fetch by Ensembl transcript ID" menu option
uses. `NORAD_mouse.fasta` used the same underlying Ensembl REST endpoint
directly, since `fetch_enst_cdna()` currently only accepts `ENST`-prefixed
(human) IDs, not `ENSMUST` (mouse) — see the note in `io.py`.

## Accessing the bundled files at runtime

Use `rna_lexis.io.example_dataset_path(name)` (`name` is `"NORAD_human"` or
`"NORAD_mouse"`), which resolves the path correctly whether the package is
installed from a wheel, an editable install, or run directly from a repo
checkout. Don't hardcode a path relative to this file — it may not exist
as a plain file on disk once installed (e.g. inside a zipped wheel).

Note `rna_lexis.io.read_text()` strips non-alphabetic characters but is
**not** FASTA-aware — it will merge the `>` header line's words into the
sequence. Either strip the header yourself first, or (from the interactive
menu) use "Load example dataset (NORAD)", which handles this correctly via
the same FASTA-parsing path as "Load from local file".

## Reproducing / updating

```python
from rna_lexis.io import fetch_enst_cdna
tid, label, seq = fetch_enst_cdna("ENST00000565493")
```

For the mouse sequence, use `_fetch_ensembl_json` directly against
`https://rest.ensembl.org/sequence/id/ENSMUST00000192863?type=cdna` (see
above for why `fetch_enst_cdna` itself won't accept this ID).

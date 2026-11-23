# Calibration of the background for whole-unit analysis

This folder holds the simulations behind the defaults in `docs/scoping_codon_units.md`
(section 7): which background model, and which multiple-testing correction, should the
statistical scoring use when the units are codons or letters instead of nucleotides?
Nothing here is shipped in the package.

Everything is reproducible. The tables quoted below are generated into
[results/SUMMARY.md](results/SUMMARY.md) by `summarize.py` from the CSV files in `results/`.

## Result in one paragraph

Use **order 0 with a flat pseudocount of 1 for codons** and **order 1 with a flat pseudocount of 1 for
letter alphabets**, and correct for multiple testing over **every possible word of each length** (K^k
words for an alphabet of K symbols), not over the repeated words only. The shrinkage model proposed in
the scoping doc is not needed; masking a word's own occurrences, a fixed reference codon usage, and
order 2 and above are all anti-conservative and are not recommended as defaults.

## Method

The scoring step of the CLI's `rank-cores` pipeline is re-implemented on one-character-per-unit text in
`harness.py`, with the background model swapped. `harness.selfcheck()` shows that on nucleotide text
this reproduces `MarkovBackground` and `rank_core_candidates` exactly (same expected counts, q-values and
"statistically supported" flags). The test is the package's own: a Poisson upper tail on the observed
count against the expected count, Benjamini-Hochberg, then q <= 0.05 and enrichment >= 10 (the "strict"
rule). A "loose" rule (q <= 0.5, enrichment >= 3) is also used, because the strict rule is so
conservative that it hides bias.

**Null texts** have a known generating model, so the truth is exact:

* *Codons* (`run_codon_background.py`), N = 100, 300, 1500 codons. Null A: independent codons drawn from
  the codon usage of one of 177 real human genes. Null B: amino acids from a first-order Markov chain
  fitted to all 177 genes, codons from that gene's usage given the amino acid.
* *Letters* (`run_letters_background.py`), N = 1,000 to 50,000 letters, from four real texts (see below).
  The null is a Markov chain of order 0, 1 or 2 fitted to the text.

**Backgrounds** are estimated from the text under test, as the package does. For codons: order 0 with
pseudocount 0.1 to 2, a nucleotide-by-codon-position model, shrinkage of codon usage toward that model
(strength s = 4 to 4096), the same with the tested word's own occurrences masked out, order 1 (flat, and
shrunk toward order 0), a fixed reference usage pooled from the other genes, and the true model as an
upper bound. For letters: order 0 to 3 with pseudocount 1.

**Measured**

1. Error against the truth: KL per word, and the share of probability mass that is under-estimated 5x.
2. Bias of the expected count on the words seen at least twice, which are the only ones that can be called
   significant (mean log of expected over true expected; negative means the background expects too little
   and chance repeats look significant).
3. False positives per null text, for three multiple-testing families.
4. Recall of a planted phrase.

Real texts and genes are used only descriptively: they contain real repeats, so a "false positive" is not
defined on them.

## Findings

### 1. The multiple-testing family matters more than the background

`score_exact_motifs` runs Benjamini-Hochberg over the candidate list it is given, and that list is already
selected for being repeated. On codon null text this flags chance repeats in 84-100% of texts for the order-0
backgrounds, **even when the true generating model is the background** (SUMMARY, "false positives per null text, by family").
Counting every word seen at least once does not fix it. Correcting over the possible words does.
Two versions were run. Over the lengths 2 to 6 taken together (`possible`) the false-positive rate is 0.00
for every sound background, but the total is dominated by the longest length (64^6 words) and short words
are almost never called: a 4-codon phrase with 3 copies is recovered in 39% of texts at 100 codons.
**Correcting within each word length k, over that length's K^k possible words (`bylength`),** has a false-positive
rate of at most 2.5% per null text even with the true model as the background (under the nominal 5%), and
recovers the same phrase in 99%. That is the family to use for codons and letters, and the one implemented.

For the existing nucleotide pipeline the discovery step hides most of this when it is strict.
`nucleotide_cli_check.py` runs the CLI's own two calls (`_discover_xmotifs`, `rank_core_candidates`) on iid
nucleotide text with the composition of a human coding sequence, 300 texts per length. With the CLI
defaults (`--min-xmotif-len 7 --min-count 2 --candidate-min-len 5`) at least one core is called supported in
12% of 900-nt texts, 31% of 2,000-nt texts and 49% of 4,500-nt texts. With `--min-xmotif-len 8 --min-count 3
--candidate-min-len 7`, the settings used in the whale and finch analyses, none of the 300 texts at any
length had a supported core, and none had any candidate. The trade-off is that motifs occurring only twice
or shorter than 7 cannot be found. The decision was to leave nucleotide mode as it is and document the CLI
defaults (README, user guide, `llms.txt`, `CHANGELOG.md`).

### 2. Codons: order 0, flat pseudocount 1

* **Error.** The flat pseudocount of 1 is at or near the best at every length (KL 0.54, 0.26, 0.064 at
  100, 300, 1,500 codons for null A). The shrinkage hybrid with s = 64 to 256 is equivalent at best and
  marginally better at 100 codons (0.53). The scoping doc's concern, that 64 pseudocounts distort a short
  transcript, is not borne out.
* **Strong shrinkage and the position-only model are worse.** s >= 1024 and the nucleotide-by-position
  model have KL 0.16 to 0.6 at 1,500 codons against 0.06. Their higher recall on planted repeats comes
  from under-estimating expected counts, not from skill.
* **Order 1 is much worse**, even when the truth has amino-acid-level structure: KL 0.47 to 0.54 at 1,500
  codons against 0.06 to 0.09, and almost no power. 4,096 cells cannot be estimated from a few hundred
  codons.
* **Masking is anti-conservative.** Removing a word's own occurrences from the counts raises recall
  to near the oracle, but on the words seen at least twice it makes the expected count at least 2x too
  small in 48% of cases at 100 codons and 21% at 300 (flat order 0 without masking: 6% and 2%, with a
  slightly positive bias, which is the safe direction).
* **A pooled reference usage is the worst** (61% of repeated words under-estimated 2x at 100 codons for null
  A): a gene's usage differs from the genome's.
* **Recall** at the recommended background: a phrase of 6 codons with 2 copies is recovered in 97%, 94%
  and 80% of texts at 100, 300, 1,500 codons; 4 codons with 4 copies in 97%, 93%, 84%.

### 3. Letters: order 1, flat pseudocount 1

* **Order 1 is calibrated for iid and order-1 text at every length tested.** Strict false-positive rate
  0.00 (worst corpus at most 0.01) for N = 1,000 to 50,000.
* **Orders 0, 2 and 3 are not.** Order 0 and 3 reach 1.00 on order-1 text from 8,000 letters. Order 2
  reaches 0.83 on order-1 text and 0.41 on order-2 text at 50,000 letters: with a flat pseudocount its
  17,576 contexts stay too sparse, and sparse contexts under-estimate the frequent transitions.
* **Order 1 is not enough for real language.** On order-2 text it fails from about 3,000 letters. That
  is expected: words are structure beyond letter pairs, and finding them is the purpose. It means "significant"
  here reads as "repeated beyond what letter-pair frequencies explain".
* **Overshooting the order is safer than undershooting it** for iid text (order 1 on iid text is
  conservative, false positives 0.00) but not once contexts become sparse (order 2 above).
* Recall for planted short phrases at order 1 is low (for example 4 letters with 8 copies: 0.21 at 1,000
  letters), because the correction over 26^k words is strict and a frequent short word is largely
  explained by letter-pair frequencies.

### 4. The rule for the order

Order 1 when N >= K^2 (at least one observation per transition on average), otherwise order 0, and never
above 1 by default. This matches every case tested: codons (K^2 = 4,096, N <= 1,500, so order 0), letters
(K^2 = 576 to 676, N >= 1,000, so order 1) and nucleotides (K^2 = 16, so order 1, as today).
Not tested: codon text with N >= 4,096, where the rule would switch to order 1.

### 5. Discovery does not produce candidates on codon null text

The CLI pipeline (xmotifs, then substrings shared by two distinct xmotifs, then scoring) found candidates
in at most 0.7% of codon null texts for xmotif minimum length 2 to 4, N = 100 to 1,500, so its false-positive
rate there is essentially zero for any background. The discovery parameters were not calibrated beyond
this; that needs planted *families* (a core with variable extensions), which exact repeats do not make.

### 6. Discovery settings per alphabet (phase 5)

`run_pipeline_defaults.py` runs the whole CLI pipeline (`find_boundary`, then `rank_core_candidates` with
`min_xmotif_type_support=2`) on null texts and on planted core families, for codons, Latin letters and amino
acids, over a grid of (shortest xmotif, shortest core candidate, fewest occurrences). Planted families are a core
of 3 to 8 units with 4 or 6 copies, each followed by one of two different 2-unit extensions (the pipeline needs two
distinct xmotifs that share a core). Results (full table in `results/SUMMARY.md`):

* **False positives are zero in every setting, alphabet and null model** (80 texts per cell, per length). They do not
  limit the choice. This is the pipeline's strictness: two distinct repeated words that share a core, scored against
  all `K**L` possible words, essentially never arise by chance in these nulls.
* **Recall is what the settings change**, and it falls in two ways. A minimum core length above the planted core, or a
  minimum xmotif length above the core plus its extension, cannot find the family. And **a minimum of 3 occurrences
  halves recall**: a family with 4 copies split between two extensions has only 2 copies of each xmotif.
* **The defaults** follow from that. The shortest core candidate is the shortest core one wants to find, the shortest
  xmotif is that plus 2, and the minimum count stays at 2: **codon and amino acid 5, 3, 2; Latin letters 6, 4, 2**
  (shortest xmotif, shortest core, fewest occurrences). Settings above these lose recall; settings below them recover
  at most about 1% more here (0.981 against 0.974 for codons) and would mostly add short repeats in real text. Any alphabet without calibrated values gets 5, 3, 3, 2
  if it has more than 8 symbols and the nucleotide values otherwise.
* **Caveat.** The null texts have no real repeats. Real text has many, and the pipeline reports them (that is its
  purpose), so zero false positives here does not mean few supported cores in a real document.

### 7. The order rule, beyond the first calibration

`run_order_rule_checks.py` closes two gaps (table in `results/SUMMARY.md`).

* **Codon text above 4,096 codons**, where the rule switches to order 1. On text with amino-acid structure, order 0
  starts to under-estimate (bias -0.05 to -0.07, and 5% of texts had a supported word at 8,000 codons), while order 1
  is conservative (positive bias, no false positives) and its bias shrinks as the text grows. The rule holds.
* **Amino acids** (K = 21, switch at 441 residues). Both orders are calibrated: order 0 is unbiased on independent
  residues and within 2% on a residue Markov chain, and order 1 is slightly conservative. Following the rule never
  costs safety, only a little power.

## What was not tested

* Real coding sequence as a null (real genes contain real repeats; on the 177 genes, 9% have a word
  called significant at order 0, 3% at order 1).
* Smoother estimators for order 2 and above on letters (for example shrinking each context toward the next
  lower order). The flat-pseudocount order 2 is not recommended; an interpolated one might be.
* Mutation-tolerant (Hamming) families through the full pipeline, and planted families with more than two
  distinct extensions or with mutated copies.
* Languages with spaces or punctuation as symbols, and amino-acid alphabets.

## Data and provenance

* `nucleotide_cli_check.py` also reads `data/human_cds.fasta`, for the nucleotide composition of the random text.
* `data/human_cds.fasta`: coding sequences (stop codon removed) of the canonical transcripts of 177 human
  genes from the Ensembl REST API, fetched by `fetch_cds.py` (the gene list is in the script; only symbols
  are sent). 195,989 codons. Three genes were rejected (MYC and VEGFA do not start at an ATG in their
  canonical CDS; DDX58 returned an HTTP error).
* Letter corpora are **not** in the repository. `run_letters_background.py` reads Genesis and the Constitution from the
  folder in `RNA_LEXIS_TEXTS` (default `calibration/data/texts`), the oration from `RNA_LEXIS_GENEDB` (default
  `/Users/haim/Projects/GeneDB`) and Homer from `RNA_LEXIS_HOMER` (default `/Users/haim/Desktop/Homer`):
  * Mark Antony's oration: `GeneDB/MarkAntony/data/antony_oration.txt`, 4,697 letters, 25 of 26 letters used.
  * Homer, the opening of the Odyssey in polytonic Greek (supplied by the project owner), reduced to its 24
    base letters (diacritics removed, final sigma folded into sigma), 8,123 letters.
  * `genesis.txt`: https://github.com/aparrish/rwet/blob/master/genesis.txt (an excerpt, 3,167 letters).
  * `constitution.txt`: the transcript at https://www.archives.gov/founding-docs/constitution-transcript,
    from "We the People" to the end of Article VII, 20,543 letters.

## Reproducing

```bash
python3 calibration/fetch_cds.py                       # only to refresh data/human_cds.fasta
python3 calibration/harness.py                         # self-check against the package
python3 calibration/run_codon_background.py            # about 45 s
RNA_LEXIS_TEXTS=/path/to/texts python3 calibration/run_letters_background.py --reps 120 --power-reps 60   # about 3 min
python3 calibration/nucleotide_cli_check.py            # about 1 min
python3 calibration/run_pipeline_defaults.py           # about 35 s; needs the texts folder
python3 calibration/run_order_rule_checks.py           # a few seconds
python3 calibration/summarize.py
```

Replicates are seeded, so reruns give the same tables.

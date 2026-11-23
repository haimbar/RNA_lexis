# Scoping: whole-unit (codon-level) analysis in RNA_lexis

Status: scoping only, nothing implemented. Branch: `codon-units-alphabet`.

## 1. Goal

Let the user analyze a sequence at the level of its natural units instead of
single nucleotides, while the engine works unchanged underneath. Every
position, length, mismatch count and plotted motif is then reported in units.

- **Protein-coding sequence:** the unit is the codon (64 of them), so
  synonymous codons stay distinct. Translating to amino acids is optional.
- **Non-biological sequences** (Mark Antony's oration, whale songs, finch
  songs) enter through a codebook that assigns each unit a codon label. This
  is what GeneDB already does with `eng2rna.py` and
  `birdsong_codebook.csv`.
- **lncRNA and other non-coding RNA** stay at the nucleotide level. That is
  the default and does not change.

Out of scope for now: units that are not fixed-width tuples over a base
alphabet, such as whole words with a vocabulary in the thousands. The design
below does not rule them out; see section 4.

## 2. Where the code stands today

The engine is mostly alphabet-agnostic. `find_boundary`, `cores`,
`find_with_mutations`, the neighbor, arc and coverage code, and the alignment
code all work on arbitrary strings. The nucleotide assumptions are in a short
list of places.

| Location | Assumption |
|---|---|
| `statistical.py:30` `ALPHABET = "acgt"`, `:102` `normalize_sequence()` | Silently drops every character outside `acgt`. Scoring or FDR on non-nucleotide text returns nonsense. This is the practical reason GeneDB encodes everything as codons first. |
| `statistical.py:363` `hamming_variants()`, `_variant_count()`, `:531` `_allowed_pattern_chars()` | Variant enumeration and wildcard expansion use `ALPHABET`. With 64 units the family size grows as (63 x length)^radius, so the `max_variants` guard matters much more. The wildcard characters `n . ? *` also collide with any alphabet that contains those symbols. |
| `statistical.py:36` `MarkovBackground` | Already has an `alphabet` field, but `normalize_sequence()` filters to `acgt` before it is used. |
| `algorithms.py:656` `add_mut`, `:674` `allow_mutation`; `plots.py:203, 495, 1089, 1194`; `menu.py:1296, 1420` | Upper case means "this position is a mismatch". Only works for alphabets of at most 26 caseable symbols. |
| `algorithms.py:14` `contains_only_rna`; `menu.py:235, 365, 431` `is_rna` | Gates the nucleotide-only menu items. This is the natural hook for a session-level alphabet. |
| `algorithms.py:805` `extendRNA`, `:845` `gen_hairpins`, `plots.py:640` `plot_logo`, GC-content indicator | Nucleotide-only. |
| `io.py:287` `read_text()` | Keeps `[a-z]` only (spaces optionally). |
| `io.py` `fetch_enst_cdna()` | Returns full cDNA (5'UTR + CDS + 3'UTR). Codon analysis needs the CDS and its reading frame. |
| Segment separators | Bouts and sentences are joined with `_` runs (`"_" * 9` in GeneDB). Works, but is ad hoc, and regex building is not escape-safe for arbitrary symbols. |
| Defaults (`minlen`, `minclen`, mismatch rate 1/6) | Calibrated for a 4-letter alphabet; see section 7. |

## 3. Design

### 3.1 The central idea

Map each unit to a single internal character. An `Alphabet` object holds the
symbol table and does the translation in both directions. The engine runs on
one-character-per-unit strings, so everything it measures is in units without
extra work. Only the display layer (plot labels, CSV motif columns, console
output) decodes back to codons, letters or amino acids.

Why not keep encoding to nucleotides, as GeneDB does now:

- Hamming distance is then counted in nucleotides, so a "mismatch" can be a
  third of a codon and no longer corresponds to a unit-level variant.
- Strings are three times longer, and length thresholds are in nucleotides.
- Results must be decoded afterward, and partial codons at the edges of a
  motif need special handling (`decode_cores.py`, the `lead` frame offset
  in `rnalexis_plots.py`).

### 3.2 The alphabet is declared, not observed

The set of symbols is part of the `Alphabet`, not something learned from the
text. An English document with no `z` is still analyzed over 26 letters; a
transcript that happens never to use one of the 64 codons is still analyzed
over 64. Concretely:

- Every statistic that depends on alphabet size (pseudocounts in
  `MarkovBackground`, the 4^k-style probabilities in
  `markov_pattern_probability`, variant-family sizes, wildcard expansion)
  uses the declared size and the declared symbols.
- Every declared symbol gets an internal character when the `Alphabet` is
  created, whether or not it occurs. The unused ones simply have zero
  counts, and the pseudocount gives them the right small probability.
- A symbol in the text that is not in the alphabet is handled by an
  explicit `on_unknown` option with two settings. The default is `skip`:
  the symbol is removed and a warning reports how many were skipped, which
  ones, and where. The other is `error`, which stops with the first
  offender. Either way it is never silent, unlike `normalize_sequence()`
  today. Skipping makes the two symbols on either side adjacent, which can
  create a motif occurrence that is not in the original text. The skipped
  positions are therefore recorded, so results can be mapped back to the
  original text and occurrences spanning a removal can be flagged.
- The `Alphabet` is saved with a session and written into result files, so a
  later analysis of another text can use the same one.
- **Each alphabet carries its own analysis defaults:** the background
  settings (order and pseudocount) and the motif-length and
  mismatch defaults of section 7. Built-in alphabets ship calibrated
  values. A user-defined alphabet file may specify its own; if it does
  not, the generic rule applies (order 1 when N >= K^2, otherwise order 0,
  with a flat pseudocount of 1). Some defaults are rules that depend on the
  text length, not constants. The command-line and menu options override
  them.

### 3.3 Units as fixed-width tuples, plus extra symbols

A unit is a k-tuple over a base alphabet of size V, so the unit alphabet
has V^k symbols.

- Real codons: V = 4, k = 3, so 64 units.
- Larger base alphabets work the same way (V > 4, k = 3, giving V^3 units).
  Nothing in the engine depends on V or k; they only set the size of the
  symbol table.

The alphabet can also include **extra symbols** beside the tuples: a space,
digits, punctuation, or any other single character. They are first-class
symbols like the codons: they can appear inside motifs, they count in the
background model, and the `Alphabet` knows them. This is needed for text
such as English with word spaces, or numbered items, once it is expressed
in codons. A codon-mode example: 64 codons, plus space, plus the ten
digits, for 75 symbols.

Extra symbols are distinct from **boundary symbols** (section 3.6), which
exist to keep motifs from spanning two segments.

Internal characters come from a range of Unicode code points that have no
case (for example the private use area). `lower()` and `upper()` then leave
them alone, which also removes the case-as-flag problem for alphabets larger
than 26, and avoids any clash between a user's `.` or `*` and the `n . ? *`
pattern wildcards. Raw symbols reappear only on output, escaped as needed.

### 3.4 Statistics

- `normalize_sequence()` and the `ALPHABET` constant take the `Alphabet`.
  The module-level nucleotide constant stays as the default.
- Background model (calibrated; results in `calibration/README.md`). The
  background is a Markov model of the declared units, estimated from the
  text under test with a flat pseudocount of 1 per symbol, as
  `MarkovBackground` does today.
  - **Codons: order 0**, that is the transcript's own codon usage.
  - **Letter alphabets (about 24 to 26 symbols): order 1.**
  - **Rule for the order:** order 1 when N >= K^2 (at least one observation
    per transition on average, K the declared alphabet size), otherwise
    order 0; never above 1 by default. It matches every case tested:
    codons (K^2 = 4,096, N up to 1,500: order 0), letters (K^2 = 576 to
    676, N >= 1,000: order 1) and nucleotides (K^2 = 16: order 1, as
    today). Codon text with N >= 4,096, where the rule would switch to
    order 1, was not tested.
  - **Rejected alternatives,** each measured: shrinking the counts toward a
    nucleotide-by-codon-position model (no better than the flat
    pseudocount, and strong shrinkage is worse); estimating the background
    with the tested word's own occurrences masked out (anti-conservative:
    the expected count is at least 2x too small for 48% of the repeated
    words at 100 codons); order 1 over 64 codons (KL 0.47 to 0.54 against
    0.06 to 0.09 at 1,500 codons); a fixed reference codon usage
    (anti-conservative, since a gene's usage differs from the genome's);
    and order 2 or higher for letters with a flat pseudocount (not
    calibrated even at 50,000 letters).
  - A repeated word raises its own unit counts, which makes the test
    slightly conservative, most at small N. That is the safe direction and
    is left as it is.
  - Nucleotide mode keeps order 1 as the default, so existing results do
    not change.
- **Multiple-testing family: every possible word of each length.** The
  Benjamini-Hochberg correction runs within each tested length k over that
  length's K^k possible words, not over the repeated words only. (Taking
  all lengths together as one family was also measured; it is dominated by
  the longest length and loses nearly all power on short words, so the
  per-length family is the one used.) `score_exact_motifs` corrects over the candidate list it is given,
  and that list is already selected for being repeated; on codon null text
  that flags chance repeats in 84 to 100% of texts even with the true model
  as the background. Counting every word seen at least once does not fix
  it. The per-length correction over all K^k words gives at most 2.5% false
  positives per null text, even with the true model as the background. Nucleotide mode is left as it is (section 9,
  question 9); the new family applies to the new alphabets.
- `hamming_variants` and `_variant_count` take the declared alphabet, and
  the `max_variants` guard needs a unit-level default.

### 3.5 Mutation marker

Replace "upper case means mismatch" in `add_mut`, `allow_mutation` and the
plot highlighters with an explicit marker. Regex building escapes every
symbol. This is a behavior-preserving refactor for nucleotides, and a
prerequisite for anything else.

### 3.6 Segment boundaries

Make segment boundaries a first-class concept: a list of segment spans over
the unit string, plus one reserved internal separator character that is not
part of the declared alphabet. Motifs must not span two segments, and
results should be reportable per segment. This replaces the `"_" * 9`
workaround.

### 3.7 Codon-specific pieces

- **Reading frame:** the user chooses a frame offset (0, 1, 2), or supplies a
  CDS. Trailing 1-2 nucleotides that do not form a codon are dropped, with a
  report of how many.
- **Extra symbols inside a nucleotide stream:** in the codon tokenizer, an
  extra symbol is its own unit and also a unit boundary. Codon reading
  restarts after it. A stretch between two extras whose length is not a
  multiple of 3 is reported with a warning, and its partial codon is
  dropped. A user who wants the frame to continue across a symbol strips
  it before loading. Formatting characters (line breaks, the spaces and
  numbers of a GenBank-style listing) are removed by the reader and never
  reach the tokenizer.
- **Ambiguity codes (`N`, `R`, `Y`, ...) are handled at codon level, not
  nucleotide level.** The whole triple that contains one becomes a single
  unknown unit, which `on_unknown` then skips or rejects. Skipping a lone
  nucleotide would shift the frame for every codon after it, so the
  tokenizer never does that.
- **CDS input:** a new Ensembl fetch for `type=cds`, alongside the existing
  cDNA fetch. `fetch_enst_cdna()` is human-`ENST` only today; the CLAUDE.md
  note about generalizing it applies here too.
- **Alphabet projection:** an optional many-to-one map from codons to amino
  acids (standard genetic code; other codes can be added), so the same
  analysis can be run with synonymous codons collapsed. Running both and
  comparing is the direct way to ask whether synonymous choice matters.
- **Reverse complement:** `extendRNA` stays nucleotide-only. At codon level
  the reverse complement of a codon is another codon, so this could work,
  but it is only meaningful if you also want the opposite-strand frame.
  Disabled in codon mode at first.
- **Stop codons:** treated as ordinary units by default. A flag to split
  segments at stop codons, so that motifs do not span ORFs, is a candidate.
- **Unit distance (later):** Hamming distance between units is 0 or 1 by
  default. A unit-level distance matrix (for example, the nucleotide
  Hamming distance between two codons, 0 to 3) would give partial credit
  for near-synonymous codons. The same hook gives amino-acid substitution
  matrices for the alignment code. Not part of the first release.

### 3.8 Non-biological data

An `Alphabet` can be built from a codebook (the `birdsong_codebook.csv`
format: external symbol to codon label), so the Mark Antony, whale and finch
data enter the same codon-level machinery. The codebook step is the
generalization of `eng2rna.py`. All 26 letters are in the codebook whether or
not a given text uses them (section 3.2), and spaces or digits can be added
as extra symbols.

`eng2rna.py` stays in GeneDB as the nucleotide-level comparison mode; it is
not ported. This package gets its own, more general codec instead:

- **Any external symbol set** (letters, syllable labels, digits, spaces),
  not only a-z, mapped to unit labels of width k over a base alphabet of
  size V (V = 4, k = 3 for codons; larger V or k when the symbol set needs
  more than 64 labels).
- **Codebook in and out:** read and write the codebook file (same format as
  `birdsong_codebook.csv`), so an encoding can be saved with a result and
  reused.
- **Assignment rules:** an explicit mapping, a seeded random permutation,
  or a mapping chosen to balance base frequencies (the idea behind
  `BALANCED_CODEC`). The assignment is recorded in the codebook, so it is
  reproducible without the seed.
- **Round trip:** encode, decode, and decode a motif found at an arbitrary
  offset (the partial-codon handling that GeneDB's `decode_cores.py` does).
- **Unknown symbols** follow the `on_unknown` option from section 3.2.

Since the engine works at the unit level, the codec is only needed when the
user wants the nucleotide-level view (for comparison with `eng2rna` output)
or wants to export a unit sequence as codon text.

### 3.9 Menu and batch CLI

- An alphabet choice at load time replaces the `is_rna` check, and the
  nucleotide-only items are hidden or warn. The menu lists the built-in
  named alphabets (nucleotide, codon, amino acid, Latin letters) plus
  "Load alphabet from file" for user-defined ones, which carry their own
  display name (section 9, question 6).
- `rna_lexis_stat_cli` gets an `--alphabet` option that takes a built-in
  name or an alphabet file, plus `--frame` and `--on-unknown`. Default
  behavior is unchanged.

## 4. What the design does not rule out

Because every unit is one internal character, arbitrary token vocabularies
(words, for example) fit later with only a tokenizer. The costs that would
then show up are the ones in 3.4: the variant-family enumeration and sparse
Markov tables, so that tier would need its own statistical defaults.

## 5. Function triage

| Reported in units once the alphabet is threaded through | Nucleotide-only, gated by alphabet |
|---|---|
| xmotifs, cores, `score_exact_motifs`, FDR, mutation families, gapped motifs, spacing periodicity, self-similarity arcs, neighbor plots, coverage and coverage comparison, shared motifs, k-mer p-values, alignment | hairpins, reverse complement, GC-content indicator, logo plot (for now), `ENST` fetch as it stands |

The logo plot is the one to revisit, and a plain unit-level logo is not the
answer for alphabets of about 20 or more symbols (see section 12, item 4, and
`docs/TODO-unit-level-logo.txt`): the picture needs a different visualization,
and WebLogo does not take a 64-symbol alphabet in any case.

## 6. Packaging and the alias question

Python has no R-style package alias, but there are three ways to get the
same effect.

1. **Second import name in the same distribution.** Ship two top-level
   packages from one wheel: `rna_lexis` and a neutral name such as `lexis`,
   where the second re-exports the first. Setuptools handles this directly
   (`[tool.setuptools.packages.find]`, or an explicit list). Both names
   import from one installation, versions can never drift, and entry points
   can be added under both names.
2. **Second distribution on PyPI.** A tiny package with the neutral name
   whose only dependency is `RNA_lexis`, and which just re-exports it.
   `pip install lexis` pulls in `RNA_lexis`. This keeps the two PyPI pages
   and descriptions separate, at the cost of two things to release.
3. **Invert it later.** Move the code to the neutral name and keep
   `RNA_lexis` as the thin alias distribution that depends on it. This is
   the version of option 2 to choose if the neutral name is going to be the
   main one.

Console-script entry points can be added under any name, independent of
which option is chosen. The `sys.modules` trick (making one name an alias of
another at import time) works but is fragile with tools that look at
`__file__`, so I would avoid it.

Recommendation: do nothing now. Option 1 is the cheapest to add later and
needs no change to the code layout; options 2 and 3 are the choice if
separate PyPI identities matter. Any new name needs a check that it is free
on PyPI, and a look at whether it collides with an existing import name.

## 7. Defaults and calibration

This is the part that needs science rather than plumbing. The current
defaults (xmotif `minlen` 8, core `minclen` 7, mismatch rate 1/6) are tuned
for four symbols. A 7-mer is rare by chance over 4 symbols and trivial over
64. At codon level a motif of 3 to 4 units is already statistically
interesting, but also much more fragile to a single substitution.

Plan: derive `minlen`, `minclen` and the permitted mismatches from the
declared alphabet size and text length (roughly the length at which a chance
occurrence stops being expected), then validate on the existing GeneDB
calibration sweeps (`calibration_sweep.py`, `sim_k_calibration.py`). The
simulations will need a codon-level generator with codon-usage bias, not
just uniform symbols.

**Background calibration (done).** The simulation is in `calibration/`; the
write-up, tables and data provenance are in `calibration/README.md`, and
the tables are generated into `calibration/results/SUMMARY.md`. It used 177
real human coding sequences (195,989 codons) for codon usage and amino-acid
structure, and four real texts for letters (Mark Antony's oration, the
opening of the Odyssey in Greek reduced to 24 base letters, Genesis, and the
US Constitution), with null texts of known generating model and planted
repeats. The result is the recommendation in section 3.4: order 0 for
codons, order 1 for letters, a flat pseudocount of 1, and a correction over
all possible words. The simulation also checks its own scoring against the
package's on nucleotide text.

Still open in this section:

- **Discovery parameters** (`minlen`, `minclen`, the permitted mismatches)
  were not calibrated. On codon null text the discovery step produced
  candidates in at most 0.7% of texts for an xmotif minimum length of 2 to
  4, so false positives from it are negligible, but the power of the full
  pipeline needs planted *families* (a core with variable extensions), not
  exact repeats.
- Codon text longer than 1,500 codons, amino-acid alphabets, and smoother
  estimators for letters at order 2 and above.

## 8. Testing

- The existing 67 tests are the regression guard: nucleotide mode must give
  identical output.
- Round-trip tests for each `Alphabet` and tokenizer, including partial
  trailing codons, the frame offset, and extra symbols.
- Unknown-symbol tests: `skip` removes them, warns with counts and
  positions, and maps results back to the original text; `error` stops at
  the first one.
- Codec tests: codebook write/read round trip, reproducibility from the
  codebook alone, and decoding a motif found at an arbitrary offset.
- Declared-alphabet tests: a text missing one letter (or one codon) gives
  the same background size, pseudocounts and variant families as a text
  that uses it.
- Integration tests in the style of `test_norad_integration.py`:
  - a short CDS, where a planted repeat is found at the right codon
    positions;
  - a small Mark Antony excerpt and one finch recording, checking that the
    unit-level run finds the phrases GeneDB found through codons.
- A synonymous-codon check: two sequences with the same protein but
  different codon choices should give different motifs in codon mode and the
  same motifs after projection to amino acids.
- Menu and plots still have no automated tests, so those changes need the
  usual manual and smoke checks.

## 9. Open questions

1. ~~Default background~~ Calibrated (section 3.4). The default belongs to
   the alphabet (section 3.2): codon order 0, letters order 1, nucleotides
   order 1 unchanged, all with a flat pseudocount of 1, chosen by the rule
   "order 1 when N >= K^2". A user-defined alphabet without its own
   defaults gets that rule. The results are in `calibration/README.md`.
   The shrinkage model originally proposed here was not needed. Related but
   separate: question 2 asks a different question (dependence on
   synonymous choice with the amino acid sequence held fixed), not a better
   background for the same test.
2. ~~A synonymous-shuffle null~~ Decided: wishlist, high priority, not
   in the first release (section 12, item 1). It asks a different question
   from question 1: whether a motif depends on synonymous choice with the
   amino acid sequence held fixed.
3. ~~Where `eng2rna` lives~~ Decided: it stays in GeneDB. The package gets
   a generalized codec instead (section 3.8). Remaining detail: which
   assignment rules (explicit, seeded permutation, balanced) ship first.
4. ~~Extra symbols in a raw nucleotide stream~~ Decided: codon reading
   restarts at each extra symbol, with a warning for a partial codon
   (section 3.7). The question only arises for the built-in codon alphabet
   reading raw nucleotide text. Codebook-defined alphabets work on
   pre-tokenized input, so every symbol is already a unit and there is no
   frame. Ambiguity codes are an unknown-symbol case handled at codon
   level, so the frame is kept (question 5, section 3.7).
5. ~~Unknown symbols~~ Decided: both options exist (`on_unknown='skip'` or
   `'error'`), and the default is skip with a warning (section 3.2).
6. ~~Naming on the CLI and menu~~ Decided: named alphabets. Built-ins are
   named by what they are (nucleotide, codon, amino acid, Latin letters).
   A user-defined alphabet is an alphabet file with its own display name
   (for example "finch syllables"), loaded with `--alphabet <file>` on the
   CLI and "Load alphabet from file" in the menu. The display name appears
   in plot titles and CSV headers. The alphabet file also carries the full
   declared symbol set, extra symbols included (section 3.2).
7. ~~Genetic codes~~ Decided: the standard code only for now. Other codes
   are on the wishlist (section 12).
8. ~~Packaging name~~ Decided: revisit the name when the tool is ready
   (section 6).
9. ~~Multiple-testing family in nucleotide mode~~ Decided: leave nucleotide
   mode as it is, and document the CLI defaults. Codon and letter modes
   correct over all K^k possible words (section 3.4); nucleotide mode keeps
   correcting over the candidate list, so existing results keep their
   meaning. On iid nucleotide text the CLI defaults (`--min-xmotif-len 7
   --min-count 2 --candidate-min-len 5`) call a core supported in 12% of
   900-nt, 31% of 2,000-nt and 49% of 4,500-nt texts, while 8, 3 and 7
   gave none in 300 texts at each length. This is now documented in the
   README, the user guide, `llms.txt` and the changelog. Adding the
   corrected family to nucleotide mode as an option stays possible later.

## 10. Suggested phases

Each phase should be shippable and leave nucleotide behavior unchanged.

1. `Alphabet` (declared symbols, extras, boundary character), tokenizers
   (codon with frame, codebook), segment handling, and the I/O round trip,
   with tests.
2. Statistics and algorithms take the alphabet: `normalize_sequence`,
   `MarkovBackground`, variant enumeration, gapped motifs, mutation marker
   and regex escaping. The background defaults and the multiple-testing
   family are fixed by the calibration (sections 3.4 and 7), so this phase
   can start; the nucleotide-mode family question (section 9, question 9)
   needs a decision first.
3. Display decoding: plots, CSV columns and console output.
4. Menu, batch CLI and documentation (`user_guide.md`, `README.md`,
   `llms.txt` cross-links kept intact).
5. Finish the alphabet-aware defaults and calibration (discovery parameters
   on planted families, codon text above 1,500 codons, amino-acid alphabets,
   smoother estimators for letters at order 2 and above); the amino-acid
   projection; and a unit-level view of the positions around a motif. That
   last item is **not** a plain logo: for alphabets of about 20 or more symbols
   a logo cannot be read, so the visualization needs choosing first (section 12,
   item 4).
6. Wishlist items (section 12), starting with the synonymous-codon test.

### Phase 1 status

Phase 1 is implemented in `alphabet.py`, `tokenizers.py` and `codec.py`, with
tests, and checked against GeneDB: the codec reproduces `eng2rna`'s encoding of
Mark Antony's oration exactly, and the codebook plus finch songs reproduce the
stored encoding of all 299 bouts. Choices made while implementing, where the
text above left room or was refined:

- **Internal characters.** A symbol that is a single lower-case ASCII letter
  is its own character, so nucleotide and Latin text read exactly as before.
  Every other symbol (codons, digits, spaces, syllable names) gets a private
  use character. Section 3.3 said private use for everything.
- **The boundary is always `"_"`** and never an alphabet symbol, so an extra
  symbol `_` in the text cannot collide with it.
- **Aliases** (`u` to `t`) are one character for another in codon mode, applied
  before characters are grouped. In the other modes the target must be a symbol.
- **Codebook assignment.** Explicit mappings and a seeded permutation ship; a
  base-frequency-balanced assignment is not included (section 9, question 3).
- **Segments** are joined by `gap` boundary characters (default 1) and the
  result records the segment spans. The analysis code does not use the spans
  yet, so keeping motifs from spanning a boundary is phase 2 work.

Found for phase 2: `_allowed_pattern_chars` in `statistical.py` treats `n . ? *`
as wildcards, and `n` is also a Latin letter. The wildcard must become a
character that cannot be an alphabet symbol before Latin text can go through
the gapped-motif code.

### Phase 2 status

Phase 2 is implemented in `statistical.py` (every function that depends on the
alphabet takes `alphabet=`), with 25 new tests. Without the argument the output is
unchanged: 36 outputs for the human and mouse NORAD sequences are identical to the
code before the change. Choices made while implementing:

- **Texts and motifs over an alphabet are internal strings**, and
  `normalize_sequence` checks them instead of filtering: a character that is not an
  internal character raises `ValueError`, so text that was not read with `tokenize`
  is caught. The nucleotide default still strips silently, as always.
- **Segment boundaries are respected** by the background (boundaries are not
  counted, no transition crosses one, the slots for a word are counted per
  segment), by literal positions, Hamming windows and gapped hits.
- **Order and family** are the calibrated ones: `default_markov_order` and
  `family="bylength"` (section 3.4). `markov_order=None` means order 1 for
  nucleotides. A list of candidates that leaves out some words treats the missing
  tests as p = 1, which is slightly conservative.
- **Wildcard:** `rna_lexis.alphabet.WILDCARD`, the last private use character,
  replaces `n . ? *` for patterns over an alphabet (the nucleotide code is unchanged).
- **Not calibrated for other alphabets:** the correction in `mutation_family_tests`
  (over the rows listed), and `markov_kmer_pvalues` in `algorithms.py`, which works
  from observed k-mer counts and was left alone.
- The package reproduces the calibration: on simulated codon text 0.0 false
  positives per 100 null texts (simulation: 0.0025), and 99% and 93% recall for the
  two planted cases checked (simulation: 99.5% and 93.5%).
- A bug I introduced while editing `hamming_variants` (the alphabet's characters
  shared a name with a list the loop rewrites) was caught by the before/after
  comparison, not by any existing test; the function now has tests.

Left for phases 3 and 4: decoding results for display (plots, CSV columns,
console), the menu and CLI, the nucleotide-only gating, and `algorithms.py` pieces
that assume nucleotides (`extendRNA`, `gen_hairpins`, GC content).

### Phase 3 status

Phase 3 is implemented: results are decoded where they are shown, with 20 new tests.
Without `alphabet=` every function behaves exactly as before; 20 display outputs and the
36 statistical outputs of phase 2 are identical to the code before the change.

- **Decoding happens at the edge.** Results hold internal strings;
  `Alphabet.show(internal, mask, sep, upper, pad)` turns them into text, and
  `Alphabet.unit_name` / `display_sep` set what axis text says and what goes between
  units. The plot functions, `print_core`, `print_alignment` and `write_rows_csv` take
  `alphabet=`.
- **The "upper case means this unit differs" convention carries over** to codons and
  letters by upper-casing the label of the flagged unit (`atg AAA ttt`); for a
  case-sensitive alphabet, where that would change the symbol, brackets are used.
- **Result rows over an alphabet keep the internal string** instead of the upper-cased
  copy nucleotide rows carry, so they can be fed back into the statistical functions;
  the CSV writer decodes them.
- **Layout follows the label width:** the right margin of the neighbor plot is sized
  from the displayed text, and the detail HTML shows fewer units per line when labels
  are wider (15 codons for a line that holds 60 nucleotides).
- Not done here: `plot_logo`, hairpins, GC content and `extendRNA` stay nucleotide-only
  (phase 4 gates them in the menu), and a unit-level logo is on the wishlist.

Left for phase 4: the menu and CLI (choosing an alphabet, loading a codebook, `--alphabet`,
`--frame`, `--on-unknown`), hiding or warning on the nucleotide-only items, session files
that remember the alphabet, and the user-facing documentation.

### Phase 4 status

Phase 4 is implemented: the menu, the batch CLI and the user-facing documentation, with 41
new tests (menu handlers driven by scripted answers, the CLI end to end). Nucleotide sessions
and CLI calls behave exactly as before (23 menu outputs compared before and after). Choices:

- **Where the choice is made.** When new text is loaded the menu asks how to read it
  (automatic, codons with a frame, Latin letters, an alphabet file); loading a saved
  session restores its alphabet. The CLI takes `--alphabet`, `--frame`, `--on-unknown`.
  Automatic reads non-nucleotide text as Latin letters, which changes the statistical
  items for such text (they previously dropped every character that is not A, C, G, T, U);
  this is in the changelog.
- **Typed queries are read like the text**, so codon motifs are typed as nucleotides; `.`
  stays a one-unit wildcard where the prompt allows it.
- **Per-alphabet defaults are not wired in yet.** The Markov order and the multiple-testing
  family follow the calibration, but the length defaults (minimum xmotif, core and candidate
  lengths) are still the nucleotide ones; the user guide says so and gives codon starting
  points that are not calibrated. Calibrating them, and letting an alphabet file carry its own
  `defaults`, is the remaining phase 5 work.
- **Mixed-up answers in a scripted run are not an error in the menu.** The menu has no
  automated terminal tests; the new tests replace `safe_input` with a script and check the
  handlers directly, including that nothing private-use reaches the console or a file.
- **Not done:** the interactive "Change setting" prompts do not offer per-alphabet defaults,
  the ENST fetch still returns the full transcript (use a CDS file for codon work), and
  `plot_logo`, hairpins and GC content remain nucleotide-only.

## 11. Risks

- The discovery defaults (section 7) are still uncalibrated, and the
  background and family results are from simulation with known truth. Real
  text has structure no Markov model captures, so on language "significant"
  means "repeated beyond what letter-pair frequencies explain". A wrong
  default produces plausible-looking but meaningless motifs.
- The mutation-marker refactor touches `algorithms.py`, `plots.py` and
  `menu.py`, none of which are fully covered by tests (plots and menu have
  none). It should land first, on its own, before any new feature.
- Variant-family size at 64 or more symbols can make the mutation-family
  tests unusable at radius 2 or more without a smarter enumeration.

## 12. Wishlist (deliberately deferred)

Not part of the first release. Listed so the design leaves room for them.

1. **Synonymous-codon test (high priority).** Whether a codon motif
   repeats beyond what its peptide implies, with the amino acid sequence
   held fixed. Not part of the first release, but the design keeps room for
   it: the amino-acid projection (section 3.7) and a replaceable codon table
   are the prerequisites.
   - Analytic, no shuffling. For a codon motif m, let n be the number of
     places where the peptide window equals the translation of m. Under
     the null each such window matches the exact codons with probability
     p, the product along the motif of the conditional codon usage
     q(codon | amino acid). The observed codon count is tested against
     Binomial(n, p), with the existing overlap handling.
   - Splits the Q1 enrichment into peptide-level enrichment (the same
     pipeline on the amino-acid projection) and excess concordance given
     the peptide (this test).
   - Needs at least two peptide-level occurrences, so it applies to
     repeat-rich regions. Single-codon amino acids (Met, Trp) carry no
     information; a report column for the motif's synonymous information
     (roughly the log of the number of synonymous options across its
     positions) tells the user how much the test could have shown.
   - Cannot distinguish selection from a recent duplication, which also
     keeps codons identical. The documentation should say so.
   - q(codon | amino acid) is thin for rare amino acids in one transcript,
     so item 7 below is its natural shrinkage target.
   - A shuffle simulation is useful only as a unit test of the analytic
     p-values.
   - A cheap companion that could ride along with the projection: classify
     each mismatch in an approximate hit as synonymous or non-synonymous.
     Not scheduled.
2. **Other genetic codes:** mitochondrial, and the other NCBI translation
   tables, for the codon-to-amino-acid projection (section 3.7). The first
   release has the standard code only. The projection should be a
   replaceable table, not hardcoded, so this is an addition and not a
   rewrite.
3. **Unit distance matrices:** partial credit for near-synonymous codons,
   and amino-acid substitution matrices (BLOSUM) for the alignment code
   (section 3.7).
4. **Unit-level view of the positions around a motif** (replaces the earlier
   "unit-level logo plot" item; background in `docs/TODO-unit-level-logo.txt`; a first
   version, the log-odds heatmap, is in phase 5 and its status is under
   "Phase 5 status"). A
   logo draws each position as a stack of symbols scaled by information content,
   which is capped at log2(K) bits: about 4.7 for 26 letters, 4.3 for 20 amino
   acids and 6 for 64 codons. In language or protein text the positions around a
   motif are usually close to the background, so the stacks come out as short bars
   of tiny letters that cannot be read. A second problem for skewed alphabets: a
   standard logo measures information against a uniform background, so common
   symbols (e, L) look informative when they are only common. A plain logo is
   reasonable only for a small alphabet. For about 20 symbols or more, directions
   (none chosen):
   - *relative-entropy logo:* measure each symbol against the text's own
     background, the one the statistics already use (sections 3.4 and 7), so only
     enriched symbols stand out;
   - *heatmap:* rows are symbols, columns are positions, colour is log-odds
     against the background; scales to 26 or 64 symbols and shows depleted symbols
     as well as enriched ones;
   - *top-k per position:* the few strongest symbols with the rest grouped as
     "other", plus counts so small samples are visible;
   - *collapsed alphabets:* a logo over symbol classes (vowels and consonants,
     amino-acid chemistry classes); for codons, group by amino acid with the codons
     as a second level, which connects to item 1.
   My leaning, for discussion: the log-odds heatmap against the statistical
   background as the general view, since it works for any alphabet size and uses
   the same null as the tests, with the relative-entropy logo for small alphabets.
5. **Codon-level reverse complement and stop-codon splitting** of segments
   (section 3.7).
6. **Arbitrary token vocabularies,** such as whole words (section 4).
7. **Reference codon usage:** an organism-wide codon usage table. As a
   background in place of the transcript's own it was measured and rejected
   (anti-conservative, since a gene's usage differs from the genome's), so
   the only remaining use is as a prior for the conditional codon usage in
   item 1, where the per-gene counts are thin.
8. **Smoother estimators for letters at order 2 and above:** shrink each
   context toward the next lower order. The flat-pseudocount order 2 is
   anti-conservative even at 50,000 letters; an interpolated one might not
   be (the same idea helped codon order 1 in the calibration, though not
   enough to recommend it).
9. **"Calibrate for my data" command:** for a user-defined alphabet, run
   the section 7 simulation from the user's own text (simulating from its
   fitted background) and report the false-positive rate and a suggested
   background setting. Not scheduled.
10. **Several documents, with and without concatenation.** Analyze a set of
    documents separately and jointly, and compare what each one uses. Texts in the
    same language differ in vocabulary (the US Constitution and *Winnie the Pooh*),
    and texts from different periods differ in spelling ("chuse" in the Constitution,
    "choose" in a 21st-century law), although they share much. Two ways to run a
    corpus, which answer different questions: *separately* (each document keeps its
    own xmotifs, cores and background, and the results are compared), and
    *concatenated* (one text with a segment boundary between documents, so a motif
    cannot span two of them; discovery then finds what recurs across documents, and
    the statistics use one pooled background). Things this needs, none built yet:
    - a corpus object that holds the documents, their alphabet and their segment
      boundaries, and reports for every motif in which documents it occurs and how
      often (the segment machinery of sections 3.5 and 3.6 is the base);
    - a test for whether a motif is enriched in one document relative to another
      or to the rest of the corpus (a two-sample version of the current test, with
      the correction taken over the pooled search space);
    - **variants across documents.** "chuse" and "choose" differ by an insertion, so
      the Hamming-tolerant search of the package will not link them; linking
      spelling variants needs alignment-based (edit-distance) matching, for which
      the Gotoh code is the starting point;
    - uses: classifying a new document (or a gene) by which of a reference
      collection's motifs it contains, and item 11.
11. **Detecting similar documents in a corpus.** For every pair of documents, a
    similarity from the motifs they share (for example the weighted coverage of the
    shared cores, or a Jaccard index of the core sets), then ranking, clustering or an
    ordination of the documents, and a flag for near-duplicates. Questions to settle:
    what null a "similar" pair is compared with (two documents drawn from the same
    background share some motifs by chance, more so the longer they are); how to
    compare documents of very different lengths; and cost, since all pairs is
    quadratic in the number of documents, which for a large corpus calls for an index
    (for example from the k-mers or cores of each document) instead of a pairwise
    search. Depends on item 10 for the representation. The same machinery would find
    genes with similar motif content.
12. **DNA methylation.** A fifth symbol for a methylated cytosine (say `M`), so that
    the alphabet is `a c g t m` and a methylated C is not lost or confused with an
    unmethylated one. First version, as proposed: look for the patterns that would be
    found if the Ms (or some of them) were Cs. Notes for when this is built:
    - the simplest form is a **many-to-one projection** of the alphabet (`m` becomes
      `c`), which is the same mechanism as the codon-to-amino-acid projection of
      phase 5; "some of them" needs a pattern that matches either symbol (a class
      `[cm]`), which is a different feature (matching under an equivalence);
    - the statistics need a decision: a background over five symbols, or over the
      projected four, and which the declared alphabet size K is for the correction;
    - methylation is largely at CpG, and on the opposite strand a methylated C
      appears as the G of a CpG, so reverse complements and strand handling need
      thought (a symbol for the opposite-strand mark, or a rule);
    - where the data come from (bisulfite or enzymatic conversion, nanopore calls) and
      what is a methylated "letter" when only a fraction of reads is methylated: a
      threshold, or a per-position fraction that the symbols cannot hold.

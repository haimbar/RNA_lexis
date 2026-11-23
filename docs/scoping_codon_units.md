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
- Background model (recommended default, pending the calibration in
  section 7; see section 9, question 1). Order 0 over the declared units,
  which for codons is the transcript's own codon usage, with the counts
  shrunk toward a simpler model instead of a flat pseudocount. For a codon
  alphabet the shrinkage target is a nucleotide model by codon position;
  for other alphabets it is uniform. The estimate is
  (count + s * q) / (N + s), where q is the target and s its strength. A
  flat pseudocount of 1 over 64 units would be about 18% of the mass for a
  300-codon transcript, which is why it is not used.
  - Higher orders are available but not the default. Order 1 over 64 units
    has 4,095 free parameters, which overfits a transcript of a few hundred
    codons and absorbs real repeats into the background.
  - When the user does not set the order, it is the highest one for which
    |alphabet|^(order + 1) <= N / 5. This gives order 0 for 64 codons at
    N = 400, and order 1 for 26 letters at N ~ 10,000.
  - Nucleotide mode keeps order 1 as the default, so existing results do
    not change.
  - The background is estimated from the text being tested, so a repeated
    motif inflates its own unit counts. The effect makes the test
    conservative and grows with alphabet size and shorter text. Estimating
    with xmotif occurrences masked is a candidate fix; the calibration
    checks whether it is needed.
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
  extra symbol is its own unit. Whether codon reading restarts after it, or
  continues across it, needs a decision (section 9).
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

The logo plot is the one to revisit: a unit-level logo (64 symbols, or 20
after projection to amino acids) is meaningful and would need its own
renderer, since WebLogo does not take a 64-symbol alphabet.

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

**Background calibration (prerequisite for phase 2).** Before the
background default is fixed, a small simulation compares the candidates on
false-positive rate:

- Text generated from a realistic codon usage table, and from a Markov
  amino-acid model with codon choice, so there is real structure but no
  planted motifs.
- Backgrounds compared: order 0 with a flat pseudocount, the nucleotide
  model by codon position, and the shrinkage hybrid recommended in
  section 3.4, each at a few transcript lengths (about 100, 300, 1,500
  codons).
- Measured: the false-positive rate of `score_exact_motifs` at the nominal
  FDR level, and the power to recover a planted repeat, with and without
  masking the motif when estimating the background.

The hybrid becomes the default only if it is calibrated at the lengths
above; otherwise the section 3.4 recommendation is revised.

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

1. **Default background for codon mode:** recommended, not yet decided.
   Order 0 over the declared units (codon usage), shrunk toward a
   nucleotide model by codon position (section 3.4). Order 1 over 64 units
   is rejected as the default because it overfits and absorbs repeats. The
   decision is confirmed or revised by the calibration in section 7, which
   must run before phase 2. Related but separate: question 2 asks a
   different question (dependence on synonymous choice with the amino acid
   sequence held fixed), not a better background for the same test.
2. **A synonymous-shuffle null:** keep the amino acid sequence, permute the
   codons within each amino acid. It tests directly whether a motif depends
   on synonymous choice. A natural fit for the motivation, but a separate
   feature.
3. ~~Where `eng2rna` lives~~ Decided: it stays in GeneDB. The package gets
   a generalized codec instead (section 3.8). Remaining detail: which
   assignment rules (explicit, seeded permutation, balanced) ship first.
4. **Extra symbols in a raw nucleotide stream:** does codon reading restart
   after each one, or continue across it? With a codebook this does not
   arise, since every unit is already a token.
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

## 10. Suggested phases

Each phase should be shippable and leave nucleotide behavior unchanged.

1. `Alphabet` (declared symbols, extras, boundary character), tokenizers
   (codon with frame, codebook), segment handling, and the I/O round trip,
   with tests.
2. Statistics and algorithms take the alphabet: `normalize_sequence`,
   `MarkovBackground`, variant enumeration, gapped motifs, mutation marker
   and regex escaping. Starts only after the background calibration in
   section 7 has fixed the default.
3. Display decoding: plots, CSV columns and console output.
4. Menu, batch CLI and documentation (`user_guide.md`, `README.md`,
   `llms.txt` cross-links kept intact).
5. Alphabet-aware defaults and calibration; amino-acid projection; unit-level
   logo.
6. Optional: unit distance matrices, synonymous-shuffle null, packaging
   alias.

## 11. Risks

- The statistical defaults (section 7) are the largest unknown. A wrong
  default produces plausible-looking but meaningless motifs.
- The mutation-marker refactor touches `algorithms.py`, `plots.py` and
  `menu.py`, none of which are fully covered by tests (plots and menu have
  none). It should land first, on its own, before any new feature.
- Variant-family size at 64 or more symbols can make the mutation-family
  tests unusable at radius 2 or more without a smarter enumeration.

## 12. Wishlist (deliberately deferred)

Not part of the first release. Listed so the design leaves room for them.

1. **Other genetic codes:** mitochondrial, and the other NCBI translation
   tables, for the codon-to-amino-acid projection (section 3.7). The first
   release has the standard code only. The projection should be a
   replaceable table, not hardcoded, so this is an addition and not a
   rewrite.
2. **Unit distance matrices:** partial credit for near-synonymous codons,
   and amino-acid substitution matrices (BLOSUM) for the alignment code
   (section 3.7).
3. **Synonymous-shuffle null:** permute codons within each amino acid to
   test whether a motif depends on synonymous choice (section 9,
   question 2).
4. **Unit-level logo plot:** its own renderer for 64 or more symbols
   (section 5).
5. **Codon-level reverse complement and stop-codon splitting** of segments
   (section 3.7).
6. **Arbitrary token vocabularies,** such as whole words (section 4).
7. **Reference codon usage as a background:** an organism-wide codon usage
   table as an alternative to the transcript's own, for short transcripts
   where the transcript-specific estimate is poor.

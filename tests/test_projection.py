"""Projections between alphabets: codons to amino acids, and other many-to-one maps."""

import collections
import random
import unittest

from rna_lexis.alphabet import BOUNDARY, Alphabet, Projection, codon, get_alphabet, nucleotide
from rna_lexis.genetic_code import (
    AMINO_ACIDS, STANDARD, STOP, TABLES, amino_acid, codon_to_amino_acid, translate_codons,
)
from rna_lexis.statistical import score_exact_motifs
from rna_lexis.tokenizers import tokenize

A = codon()


class GeneticCodeTests(unittest.TestCase):
    def test_the_standard_code(self):
        self.assertEqual(len(STANDARD), 64)
        self.assertEqual(STANDARD["atg"], "m")
        self.assertEqual(STANDARD["tgg"], "w")
        self.assertEqual({c for c, a in STANDARD.items() if a == STOP}, {"taa", "tag", "tga"})
        sizes = collections.Counter(STANDARD.values())
        self.assertEqual(sizes["l"], 6)
        self.assertEqual(sizes["s"], 6)
        self.assertEqual(sizes["r"], 6)
        self.assertEqual(sizes["m"], 1)
        self.assertEqual(sizes["w"], 1)
        self.assertEqual(sizes[STOP], 3)
        self.assertEqual(sum(sizes.values()), 64)
        self.assertEqual(set(STANDARD.values()), set(AMINO_ACIDS) | {STOP})

    def test_amino_acid_alphabet(self):
        a = amino_acid()
        self.assertEqual(a.size, 21)
        self.assertEqual(a.units_label, "residues")
        self.assertIs(a.char("m"), "m")                       # a letter is its own character
        self.assertNotEqual(a.char(STOP), STOP)               # the stop sign is not a letter
        self.assertEqual(get_alphabet("amino-acid").to_dict(), a.to_dict())

    def test_translation(self):
        t = tokenize("ATGGCCTAAGGGTGA", A)
        p = t.project(codon_to_amino_acid())
        self.assertEqual(p.decode(sep=" "), "m a * g *")
        self.assertEqual(translate_codons(t.text), p.text)
        self.assertEqual(p.alphabet.name, "amino acid")
        with self.assertRaises(ValueError):
            codon_to_amino_acid("martian")

    def test_synonymous_codons_project_to_one_residue(self):
        proj = codon_to_amino_acid()
        t = A.encode(["ctt", "cta", "ctg", "tta", "ttg", "ctc"])
        self.assertEqual(len(set(proj.apply(t))), 1)          # all six leucine codons
        self.assertEqual(len(set(t)), 6)
        self.assertEqual(sorted(proj.groups()["l"]), sorted(["ctt", "cta", "ctg", "tta", "ttg", "ctc"]))

    def test_positions_segments_and_skips_are_unchanged(self):
        t = tokenize(["atgaaa", "tttNNNggg"], A)
        p = t.project(codon_to_amino_acid())
        self.assertEqual(len(p.text), len(t.text))
        self.assertEqual(p.offsets, t.offsets)
        self.assertEqual([(s.start, s.end) for s in p.segments], [(s.start, s.end) for s in t.segments])
        self.assertEqual(len(p.skipped), len(t.skipped))
        self.assertIn(BOUNDARY, p.text)
        self.assertEqual(p.decode(sep=" "), "m k | f g")

    def test_a_tokenized_text_of_another_alphabet_is_refused(self):
        with self.assertRaises(ValueError):
            tokenize("acgt", nucleotide()).project(codon_to_amino_acid())


class ProjectionTests(unittest.TestCase):
    def test_every_symbol_must_be_mapped_to_a_symbol_of_the_target(self):
        src = Alphabet(name="s", symbols=("a", "b", "c"))
        dst = Alphabet(name="d", symbols=("x", "y"))
        with self.assertRaises(ValueError):
            Projection(src, dst, {"a": "x", "b": "y"})            # c is not mapped
        with self.assertRaises(ValueError):
            Projection(src, dst, {"a": "x", "b": "y", "c": "z"})  # z is not in the target
        p = Projection(src, dst, {"a": "x", "b": "x", "c": "y"})
        self.assertEqual(dst.decode(p.apply(src.encode("abcab"))), "xxyxx")

    def test_extras_are_kept_when_the_target_declares_them(self):
        src = Alphabet(name="s", symbols=("a", "b"), extras=(" ",))
        dst = Alphabet(name="d", symbols=("x",), extras=(" ",))
        p = Projection(src, dst, {"a": "x", "b": "x"})
        self.assertEqual(dst.decode(p.apply(src.encode(["a", " ", "b"]))), "x x")
        with self.assertRaises(ValueError):
            Projection(src, Alphabet(name="d2", symbols=("x",)), {"a": "x", "b": "x"})

    def test_foreign_characters_are_refused(self):
        p = codon_to_amino_acid()
        with self.assertRaises(ValueError):
            p.apply("atg")                                        # external text, not an internal string

    def test_methylation_read_as_cytosine(self):
        # a fifth letter m for a methylated c, analyzed "as if" it were c (wishlist item 12)
        five = Alphabet(name="DNA with methylation", symbols=tuple("acgtm"), unit_name="nt")
        four = nucleotide()
        proj = Projection(five, four, {"a": "a", "c": "c", "g": "g", "t": "t", "m": "c"},
                          name="methylated C as C")
        text = tokenize("gtmgtmgtc" * 3, five).text
        self.assertEqual(text.count("gtm"), 6)                    # as written: two kinds of word
        projected = proj.apply(text)
        self.assertEqual(four.decode(projected), "gtcgtcgtc" * 3)
        self.assertEqual(projected.count("gtc"), 9)               # as one pattern
        rows = score_exact_motifs(projected, ["gtc"], alphabet=four, enrichment_threshold=1.0)
        self.assertEqual(rows[0]["exact_count"], 9)

    def test_statistics_over_amino_acids_use_the_declared_size(self):
        a = amino_acid()
        rnd = random.Random(1)
        units = [rnd.choice(AMINO_ACIDS) for _ in range(300)]
        units[10:16] = list("mkvlay")
        units[100:106] = list("mkvlay")
        units[200:206] = list("mkvlay")
        text = a.encode(units)
        row = score_exact_motifs(text, [a.encode("mkvlay")], alphabet=a)[0]
        self.assertEqual(row["exact_count"], 3)
        self.assertAlmostEqual(row["q_markov"], min(1.0, row["p_markov"] * 21 ** 6), places=12)
        self.assertEqual(row["markov_order"], 0)                  # 300 < 21**2 residues


if __name__ == "__main__":
    unittest.main()

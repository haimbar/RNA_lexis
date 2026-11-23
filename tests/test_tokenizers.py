import unittest
import warnings

from rna_lexis.alphabet import BOUNDARY, Alphabet, codon, latin, nucleotide
from rna_lexis.tokenizers import (
    PartialUnitWarning, UnknownSymbolError, UnknownSymbolWarning, tokenize,
)


def quiet(fn, *args, **kwargs):
    """Call fn and return (result, list of warning categories)."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = fn(*args, **kwargs)
    return out, [x.category for x in w]


class CharsModeTests(unittest.TestCase):
    def test_nucleotides_read_as_before(self):
        t, w = quiet(tokenize, "ACGU acgt\nacgt", nucleotide())
        self.assertEqual(t.text, "acgtacgtacgt")
        self.assertEqual(w, [])

    def test_whitespace_is_formatting_unless_declared(self):
        a = latin()
        self.assertEqual(tokenize("to be\nor", a).text, "tobeor")
        b = a.with_extras(" ")
        t = tokenize("to be", b)
        self.assertEqual(len(t.text), 5)
        self.assertEqual(b.decode(t.text), "to be")

    def test_an_unused_letter_is_still_part_of_the_alphabet(self):
        t = tokenize("the quick brown fox", latin())
        self.assertNotIn("z", t.text)
        self.assertEqual(t.alphabet.size, 26)

    def test_offsets_point_into_the_original_text(self):
        t = tokenize("a c g", nucleotide())
        self.assertEqual(t.offsets, [0, 2, 4])

    def test_digits_are_unknown_unless_declared(self):
        t, w = quiet(tokenize, "ab1", latin())
        self.assertEqual(t.text, "ab")
        self.assertEqual(w, [UnknownSymbolWarning])
        t = tokenize("ab1", latin().with_extras("0123456789"))
        self.assertEqual(len(t.text), 3)


class UnknownSymbolTests(unittest.TestCase):
    def test_skip_is_the_default_and_warns_with_details(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            t = tokenize("acgxnacg", nucleotide())
        self.assertEqual(t.text, "acgacg")
        self.assertEqual([s.text for s in t.skipped], ["x", "n"])
        self.assertEqual(len(w), 1)
        msg = str(w[0].message)
        self.assertIn("2 unknown symbol(s)", msg)
        self.assertIn("'x'", msg)
        self.assertIn("offset 3", msg)

    def test_error_stops_at_the_first_unknown_symbol(self):
        with self.assertRaises(UnknownSymbolError) as cm:
            tokenize("acgxnacg", nucleotide(), on_unknown="error")
        self.assertIn("'x'", str(cm.exception))

    def test_invalid_option_is_rejected(self):
        with self.assertRaises(ValueError):
            tokenize("acg", nucleotide(), on_unknown="ignore")

    def test_removal_points_flag_occurrences_across_a_removal(self):
        t, _ = quiet(tokenize, "aaxaa", Alphabet(name="x", symbols=("a",)))
        self.assertEqual(t.text, "aaaa")
        self.assertEqual(t.removal_points, [2])
        self.assertTrue(t.spans_removal(1, 2))     # units 1,2 were not adjacent
        self.assertTrue(t.spans_removal(0, 4))
        self.assertFalse(t.spans_removal(0, 2))
        self.assertFalse(t.spans_removal(2, 2))

    def test_no_warning_when_nothing_is_skipped(self):
        _, w = quiet(tokenize, "acgt", nucleotide())
        self.assertEqual(w, [])


class FixedModeTests(unittest.TestCase):
    def test_reads_codons(self):
        t = tokenize("atgaaattt", codon())
        self.assertEqual(t.decode(sep=" "), "atg aaa ttt")
        self.assertEqual(t.n_units, 3)

    def test_rna_is_read_as_dna_and_case_is_ignored(self):
        t = tokenize("AUGAAAUUU", codon())
        self.assertEqual(t.decode(sep=" "), "atg aaa ttt")

    def test_frame_drops_leading_characters_and_records_them(self):
        t = tokenize("ccatgaaa", codon(), frame=2)
        self.assertEqual(t.decode(sep=" "), "atg aaa")
        self.assertEqual([(d.reason, d.text) for d in t.dropped], [("frame", "cc")])
        self.assertEqual(t.offsets, [2, 5])

    def test_frame_must_be_smaller_than_the_width(self):
        with self.assertRaises(ValueError):
            tokenize("atgaaa", codon(), frame=3)
        with self.assertRaises(ValueError):
            tokenize("acgt", nucleotide(), frame=1)

    def test_trailing_partial_unit_is_dropped_with_a_warning(self):
        t, w = quiet(tokenize, "atgaaat", codon())
        self.assertEqual(t.decode(sep=" "), "atg aaa")
        self.assertEqual(w, [PartialUnitWarning])
        self.assertEqual([(d.reason, d.text) for d in t.dropped], [("partial", "t")])

    def test_ambiguity_code_makes_the_whole_codon_unknown_and_keeps_the_frame(self):
        # skipping the lone 'n' as a character would shift every later codon
        t, w = quiet(tokenize, "atgancaaattt", codon())
        self.assertEqual(t.decode(sep=" "), "atg aaa ttt")
        self.assertEqual([s.text for s in t.skipped], ["anc"])
        self.assertEqual(w, [UnknownSymbolWarning])
        self.assertEqual(t.removal_points, [1])

    def test_ambiguity_code_with_error_option(self):
        with self.assertRaises(UnknownSymbolError):
            tokenize("atgancaaa", codon(), on_unknown="error")

    def test_extras_restart_the_reading(self):
        a = codon().with_extras(" ")
        # "atg" + space + "aaa": the extra is a unit of its own
        t = tokenize("atg aaa", a)
        self.assertEqual(a.decode_list(t.text), ["atg", " ", "aaa"])
        # a run that is not a multiple of 3 loses its tail, and the next run
        # starts fresh at its own first character
        t, w = quiet(tokenize, "atgaa ttt", a)
        self.assertEqual(a.decode_list(t.text), ["atg", " ", "ttt"])
        self.assertEqual(w, [PartialUnitWarning])
        self.assertEqual([d.text for d in t.dropped], ["aa"])

    def test_whitespace_inside_a_run_is_formatting_when_not_declared(self):
        t = tokenize("at gaa\naTTT", codon())
        self.assertEqual(t.decode(sep=" "), "atg aaa ttt")

    def test_frame_applies_to_the_first_run_of_each_segment_only(self):
        a = codon().with_extras(" ")
        t = tokenize(["cATGaaa tttggg"], a, frame=1)
        self.assertEqual(a.decode_list(t.text), ["atg", "aaa", " ", "ttt", "ggg"])

    def test_larger_base_alphabet_works_the_same_way(self):
        # V = 5 base symbols, width 3: units are triples over 5 letters
        import itertools
        units = ["".join(p) for p in itertools.product("abcde", repeat=3)]
        a = Alphabet(name="V5", symbols=tuple(units), mode="fixed", width=3)
        self.assertEqual(a.size, 125)
        t = tokenize("abcdeaedc", a, frame=0)
        self.assertEqual(t.decode(sep=" "), "abc dea edc")


class SegmentTests(unittest.TestCase):
    def test_segments_are_joined_by_the_boundary_character(self):
        t = tokenize(["acg", "tt"], nucleotide())
        self.assertEqual(t.text, "acg" + BOUNDARY + "tt")
        self.assertEqual([(s.start, s.end) for s in t.segments], [(0, 3), (4, 6)])
        self.assertEqual(t.offsets, [0, 1, 2, -1, 0, 1])
        self.assertEqual(t.n_units, 5)

    def test_gap_sets_the_number_of_boundary_characters(self):
        t = tokenize(["a", "c"], nucleotide(), gap=3)
        self.assertEqual(t.text, "a___c")
        with self.assertRaises(ValueError):
            tokenize(["a", "c"], nucleotide(), gap=0)

    def test_each_segment_has_its_own_frame_and_reading(self):
        t, w = quiet(tokenize, ["atgaaat", "gggtt"], codon())
        self.assertEqual(t.decode(sep=" "), "atg aaa | ggg")
        self.assertEqual([d.segment for d in t.dropped], [0, 1])

    def test_removals_are_reported_in_whole_text_positions(self):
        t, _ = quiet(tokenize, ["ac", "agxt"], nucleotide())
        self.assertEqual(t.text, "ac_agt")
        self.assertEqual(t.removal_points, [5])

    def test_empty_input(self):
        t = tokenize("", nucleotide())
        self.assertEqual(t.text, "")
        self.assertEqual(t.n_units, 0)


class WordsModeTests(unittest.TestCase):
    def setUp(self):
        self.a = Alphabet(name="syllables", symbols=("ba", "tu", "ko"), mode="words")

    def test_reads_whitespace_separated_labels(self):
        t = tokenize("ba tu\nko ba", self.a)
        self.assertEqual(t.decode(sep=" "), "ba tu ko ba")
        self.assertEqual(t.offsets, [0, 3, 6, 9])

    def test_unknown_word_is_skipped_whole(self):
        t, w = quiet(tokenize, "ba xx tu", self.a)
        self.assertEqual(t.decode(sep=" "), "ba tu")
        self.assertEqual([s.text for s in t.skipped], ["xx"])
        self.assertEqual(w, [UnknownSymbolWarning])


if __name__ == "__main__":
    unittest.main()

import os
import tempfile
import unittest

from rna_lexis.alphabet import (
    Alphabet, BOUNDARY, codon, get_alphabet, latin, nucleotide,
)


class DeclaredSetTests(unittest.TestCase):
    def test_codon_alphabet_has_all_64_codons(self):
        a = codon()
        self.assertEqual(a.size, 64)
        self.assertEqual(len(set(a.chars)), 64)
        self.assertIn("atg", a)
        self.assertNotIn("atn", a)

    def test_size_counts_declared_symbols_not_observed_ones(self):
        # no 'z' needs to occur for the alphabet to have 26 symbols
        self.assertEqual(latin().size, 26)
        self.assertEqual(nucleotide().size, 4)

    def test_extras_count_in_the_size(self):
        a = codon().with_extras(" 0123456789")
        self.assertEqual(a.size, 64 + 11)
        self.assertEqual(a.extras[0], " ")

    def test_with_extras_does_not_change_the_original_or_duplicate(self):
        a = codon()
        b = a.with_extras("  ")
        self.assertEqual(a.extras, ())
        self.assertEqual(b.extras, (" ",))


class InternalCharacterTests(unittest.TestCase):
    def test_single_lower_case_letters_are_their_own_character(self):
        a = nucleotide()
        self.assertEqual(a.chars, "acgt")
        self.assertEqual(latin().chars, "abcdefghijklmnopqrstuvwxyz")

    def test_other_symbols_use_the_private_use_area_and_have_no_case(self):
        a = codon().with_extras(" 0")
        for label in a.labels:
            c = a.char(label)
            self.assertEqual(len(c), 1)
            self.assertEqual(c.lower(), c)
            self.assertEqual(c.upper(), c)
            self.assertNotEqual(c, BOUNDARY)

    def test_mapping_is_deterministic_and_round_trips(self):
        a, b = codon(), codon()
        self.assertEqual(a.chars, b.chars)
        for label in a.labels:
            self.assertEqual(a.label(a.char(label)), label)

    def test_encode_decode(self):
        a = codon()
        units = ["atg", "aaa", "ttt"]
        s = a.encode(units)
        self.assertEqual(len(s), 3)
        self.assertEqual(a.decode_list(s), units)
        self.assertEqual(a.decode(s, sep=" "), "atg aaa ttt")

    def test_boundary_is_shown_distinctly_when_decoding(self):
        a = latin()
        self.assertEqual(a.decode("ab_cd"), "ab|cd")
        self.assertEqual(a.decode_list("a_b"), ["a", BOUNDARY, "b"])

    def test_unknown_label_raises(self):
        with self.assertRaises(KeyError):
            codon().char("atn")


class ValidationTests(unittest.TestCase):
    def test_duplicate_symbols_are_rejected(self):
        with self.assertRaises(ValueError):
            Alphabet(name="x", symbols=("a", "b", "a"))

    def test_upper_case_labels_need_case_sensitive(self):
        with self.assertRaises(ValueError):
            Alphabet(name="x", symbols=("A", "b"))
        a = Alphabet(name="x", symbols=("A", "a"), case_sensitive=True)
        self.assertNotEqual(a.char("A"), a.char("a"))

    def test_extras_must_be_single_characters(self):
        with self.assertRaises(ValueError):
            Alphabet(name="x", symbols=("a",), extras=("ab",))

    def test_mode_must_match_the_symbol_width(self):
        with self.assertRaises(ValueError):
            Alphabet(name="x", symbols=("aa", "b"), mode="chars")
        with self.assertRaises(ValueError):
            Alphabet(name="x", symbols=("aaa", "bb"), mode="fixed", width=3)

    def test_codes_must_be_distinct_and_equal_in_length(self):
        with self.assertRaises(ValueError):
            Alphabet(name="x", symbols=("a", "b"), codes={"a": "ttt", "b": "ttt"})
        with self.assertRaises(ValueError):
            Alphabet(name="x", symbols=("a", "b"), codes={"a": "ttt", "b": "tt"})
        with self.assertRaises(ValueError):
            Alphabet(name="x", symbols=("a",), codes={"z": "ttt"})

    def test_alias_target_must_be_a_symbol_except_in_fixed_mode(self):
        with self.assertRaises(ValueError):
            Alphabet(name="x", symbols=("a", "b"), aliases={"u": "t"})
        # fixed mode: an alias is one character for another
        Alphabet(name="x", symbols=("aaa",), mode="fixed", width=3, aliases={"u": "t"})

    def test_empty_alphabet_is_rejected(self):
        with self.assertRaises(ValueError):
            Alphabet(name="x", symbols=())


class FileTests(unittest.TestCase):
    def test_save_and_load_round_trip(self):
        a = codon().with_extras(" 01", name="codon + word space")
        a = a._replace(defaults={"order": 0, "min_len": 4})
        with tempfile.TemporaryDirectory() as d:
            path = a.save(os.path.join(d, "a.json"))
            b = Alphabet.load(path)
        self.assertEqual(b.name, a.name)
        self.assertEqual(b.symbols, a.symbols)
        self.assertEqual(b.extras, a.extras)
        self.assertEqual(b.codes, a.codes)
        self.assertEqual(b.mode, "fixed")
        self.assertEqual(b.defaults, {"order": 0, "min_len": 4})
        self.assertEqual(b.chars, a.chars)

    def test_get_alphabet_by_name_and_by_file(self):
        self.assertEqual(get_alphabet("Codon").size, 64)
        self.assertEqual(get_alphabet("nucleotide").size, 4)
        with tempfile.TemporaryDirectory() as d:
            path = latin().save(os.path.join(d, "l.json"))
            self.assertEqual(get_alphabet(path).size, 26)
        with self.assertRaises(ValueError):
            get_alphabet("no such alphabet")


if __name__ == "__main__":
    unittest.main()

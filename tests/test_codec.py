import os
import tempfile
import unittest
import warnings

from rna_lexis.alphabet import BOUNDARY, Alphabet, codon, latin
from rna_lexis.codec import (
    code_width, decode_nt, encode_nt, frame_for_position, make_codebook,
    read_codebook_csv, write_codebook_csv,
)
from rna_lexis.tokenizers import UnknownSymbolError, tokenize


def english_alphabet(seed=1, extras=" "):
    book = make_codebook("abcdefghijklmnopqrstuvwxyz", width=3, seed=seed)
    return Alphabet.from_codebook("English as codons", book, extras=extras)


class MakeCodebookTests(unittest.TestCase):
    def test_codes_are_distinct_and_of_the_right_width(self):
        book = make_codebook("abcdefghijklmnopqrstuvwxyz", width=3, seed=7)
        self.assertEqual(len(book), 26)
        self.assertEqual(len(set(book.values())), 26)
        self.assertTrue(all(len(c) == 3 and set(c) <= set("tcga") for c in book.values()))

    def test_same_seed_gives_the_same_codebook(self):
        self.assertEqual(make_codebook("abc", seed=3), make_codebook("abc", seed=3))
        self.assertNotEqual(make_codebook("abcdefgh", seed=3), make_codebook("abcdefgh", seed=4))

    def test_the_codes_do_not_depend_on_the_order_of_the_symbols(self):
        # the shuffle is of the codes; symbols take them in order, so the
        # first symbol always gets the same code for a given seed
        a = make_codebook("abc", seed=5)
        b = make_codebook("xyz", seed=5)
        self.assertEqual(list(a.values()), list(b.values()))

    def test_too_many_symbols_for_the_width_is_an_error(self):
        with self.assertRaises(ValueError):
            make_codebook([str(i) for i in range(65)], width=3)
        make_codebook([str(i) for i in range(65)], width=4)    # 256 codes

    def test_duplicate_symbols_are_rejected(self):
        with self.assertRaises(ValueError):
            make_codebook(["a", "a"])


class CodebookFileTests(unittest.TestCase):
    def test_write_and_read_round_trip(self):
        book = make_codebook("abc", seed=2)
        with tempfile.TemporaryDirectory() as d:
            path = write_codebook_csv(os.path.join(d, "book.csv"), book)
            self.assertEqual(read_codebook_csv(path), book)

    def test_reads_a_multi_codebook_file_with_where(self):
        # the layout of GeneDB's birdsong_codebook.csv
        rows = ("bird,syllable,codon,count\n"
                "Atama,a,GCT,10\nAtama,b,CCT,5\n"
                "Other,a,AUG,7\n")
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "birds.csv")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(rows)
            book = read_codebook_csv(path, symbol_col="syllable", code_col="codon",
                                     where={"bird": "Atama"})
            self.assertEqual(book, {"a": "gct", "b": "cct"})
            other = read_codebook_csv(path, symbol_col="syllable", code_col="codon",
                                      where={"bird": "Other"})
            self.assertEqual(other, {"a": "atg"})      # upper-case read, U read as T
            with self.assertRaises(ValueError):        # 'a' has two codes without `where`
                read_codebook_csv(path, symbol_col="syllable", code_col="codon")

    def test_missing_column_is_reported(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "x.csv")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("a,b\n1,2\n")
            with self.assertRaises(ValueError) as cm:
                read_codebook_csv(path)
            self.assertIn("symbol", str(cm.exception))


class EncodeDecodeTests(unittest.TestCase):
    def test_round_trip_through_nucleotides(self):
        a = english_alphabet(extras=())
        t = tokenize("markantony", a)
        nt = encode_nt(t, a)
        self.assertEqual(len(nt), 3 * 10)
        back = decode_nt(nt, a)
        self.assertEqual(back.text, t.text)
        self.assertEqual(back.decode(), "markantony")

    def test_segments_round_trip(self):
        a = english_alphabet(extras=())
        t = tokenize(["friends", "romans"], a)
        nt = encode_nt(t, a)
        self.assertIn(BOUNDARY, nt)
        back = decode_nt(nt, a)
        self.assertEqual(back.text, t.text)
        self.assertEqual([(s.start, s.end) for s in back.segments],
                         [(s.start, s.end) for s in t.segments])

    def test_every_symbol_of_an_unused_alphabet_still_decodes(self):
        a = english_alphabet(extras=())
        for letter in "abcdefghijklmnopqrstuvwxyz":
            self.assertEqual(decode_nt(a.codes[letter], a).decode(), letter)

    def test_extras_error_keep_or_drop(self):
        a = english_alphabet(extras=" ")
        t = tokenize("to be", a)
        with self.assertRaises(ValueError):
            encode_nt(t, a)
        kept = encode_nt(t, a, extras="keep")
        self.assertEqual(len(kept), 3 * 4 + 1)
        self.assertEqual(decode_nt(kept, a).text, t.text)
        dropped = encode_nt(t, a, extras="drop")
        self.assertEqual(len(dropped), 3 * 4)
        self.assertEqual(decode_nt(dropped, a).decode(), "tobe")

    def test_decoding_a_motif_found_at_an_offset(self):
        a = english_alphabet(extras=())
        nt = encode_nt(tokenize("brutusisanhonourableman", a), a)
        pos = 3 * 5 + 1                     # a motif found one base into the 6th unit
        motif = nt[pos:pos + 9]
        frame = frame_for_position(pos, code_width(a))
        self.assertEqual(frame, 2)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            units = decode_nt(motif, a, frame=frame)
        self.assertEqual(units.decode(), "is")      # the two whole units after the frame
        # starting on a unit boundary needs no frame
        self.assertEqual(frame_for_position(15, 3), 0)

    def test_unknown_code_in_nucleotide_text(self):
        a = english_alphabet(extras=())
        free = [c for c in ("".join(p) for p in __import__("itertools").product("tcga", repeat=3))
                if c not in a.codes.values()][0]
        with self.assertRaises(UnknownSymbolError):
            decode_nt(free, a, on_unknown="error")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.assertEqual(decode_nt(free + a.codes["a"], a).decode(), "a")

    def test_the_codon_alphabet_is_its_own_code(self):
        c = codon()
        t = tokenize("atgaaattt", c)
        self.assertEqual(encode_nt(t, c), "atgaaattt")
        self.assertEqual(decode_nt("atgaaattt", c).text, t.text)

    def test_alphabet_without_codes_cannot_encode(self):
        with self.assertRaises(ValueError):
            encode_nt("abc", latin())


if __name__ == "__main__":
    unittest.main()


class FindOrfTests(unittest.TestCase):
    def test_longest_atg_to_stop_in_the_three_frames(self):
        from rna_lexis.codec import find_orf
        utr5, utr3 = "ccgcc", "gcgc"
        body = "atg" + "gaa" * 10 + "taa"
        text = utr5 + body + utr3
        got = find_orf(text)
        self.assertEqual(text[got["start"]:got["end"]], body)
        self.assertEqual((got["codons"], got["frame"]), (11, 5 % 3))
        # a longer one in another frame wins; the stop is part of the region
        longer = "ct" + "atg" + "aaa" * 30 + "tga"
        got = find_orf(body + longer)
        self.assertEqual(got["codons"], 31)

    def test_rna_case_and_nothing(self):
        from rna_lexis.codec import find_orf
        self.assertEqual(find_orf("AUGGAAUAA")["codons"], 2)
        self.assertEqual(find_orf("AUG GAA\nUAA")["end"], 9)
        self.assertIsNone(find_orf(""))
        self.assertIsNone(find_orf("atggaaaaa"))             # no stop: not an open reading frame
        self.assertIsNone(find_orf("ccccccccc"))
        self.assertIsNone(find_orf("taaatg"))                 # the stop comes first

    def test_agrees_with_a_plain_scan(self):
        import random
        from rna_lexis.codec import find_orf
        rnd = random.Random(4)
        for _ in range(200):
            t = "".join(rnd.choice("acgt") for _ in range(rnd.randint(0, 300)))
            best = 0
            for i in range(len(t) - 2):
                if t[i:i + 3] != "atg":
                    continue
                for j in range(i + 3, len(t) - 2, 3):
                    if t[j:j + 3] in ("taa", "tag", "tga"):
                        best = max(best, (j - i) // 3)
                        break
            got = find_orf(t)
            self.assertEqual(got["codons"] if got else 0, best)

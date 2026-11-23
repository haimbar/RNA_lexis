"""The faster code paths give the same answers as the plain definitions they replace."""

import random
import unittest

from rna_lexis import alphabet as A
from rna_lexis import statistical as S
from rna_lexis import tokenizers as T
from rna_lexis.algorithms import find_boundary, is_bounded, find_all_matches


def _find_boundary_reference(txt, minlen, maxlen, at_least):
    """find_boundary as it was written: scan the whole text at every length."""
    txt = txt.lower()
    n = len(txt)
    words = []
    for j in range(minlen, maxlen + 1):
        pos = {}
        for i in range(n - j + 1):
            sg = txt[i:i + j]
            if "_" in sg:
                continue
            pos.setdefault(sg, []).append(i)
        qual = {w: p for w, p in pos.items() if len(p) >= at_least}
        if not qual:
            break
        for w, p in qual.items():
            if is_bounded(w, txt, positions=p):
                words.append(w)
    return words


def _hamming_reference(motif, txt, radius):
    L = len(motif)
    hits = []
    if radius < 0 or L == 0 or L > len(txt):
        return hits
    for i in range(len(txt) - L + 1):
        w = txt[i:i + L]
        if "_" in w:
            continue
        d = sum(a != b for a, b in zip(w, motif))
        if d <= radius:
            hits.append((i, w, d))
    return hits


class FindBoundaryTests(unittest.TestCase):
    def test_same_words_in_the_same_order_as_the_plain_scan(self):
        rnd = random.Random(3)
        for trial in range(25):
            base = "".join(rnd.choice("acgt") for _ in range(rnd.randint(40, 400)))
            unit = "".join(rnd.choice("acgt") for _ in range(rnd.randint(6, 30)))
            txt = base
            for _ in range(rnd.randint(0, 5)):
                p = rnd.randrange(len(txt))
                txt = txt[:p] + unit + txt[p:]
            if trial % 3 == 0:                      # segment boundaries in the text
                cut = len(txt) // 2
                txt = txt[:cut] + "_" + txt[cut:]
            for minlen, at_least in ((3, 2), (5, 2), (4, 3), (1, 2), (2, 2)):
                with self.subTest(trial=trial, minlen=minlen, at_least=at_least):
                    self.assertEqual(find_boundary(txt, minlen=minlen, maxlen=40, at_least=at_least),
                                     _find_boundary_reference(txt, minlen, 40, at_least))

    def test_a_codon_text(self):
        cod = A.codon()
        rnd = random.Random(5)
        gene = "".join(rnd.choice("acgt") for _ in range(600))
        txt = T.tokenize(gene + gene[100:400] + gene[:200], cod).text
        self.assertEqual(find_boundary(txt, minlen=3, maxlen=60, at_least=2),
                         _find_boundary_reference(txt, 3, 60, 2))


class HammingScanTests(unittest.TestCase):
    def test_same_hits_as_the_plain_window_scan(self):
        rnd = random.Random(11)
        nuc = A.nucleotide()        # an alphabet keeps the segment boundary "_" that nucleotide mode filters out
        for trial in range(40):
            txt = "".join(rnd.choice("acg_") if rnd.random() < .05 else rnd.choice("acgt")
                          for _ in range(rnd.randint(1, 300)))
            motif = "".join(rnd.choice("acgt") for _ in range(rnd.randint(1, 8)))
            for radius in (-1, 0, 1, 2, 9):
                with self.subTest(trial=trial, radius=radius):
                    self.assertEqual(S.hamming_family_positions(motif, txt, radius, nuc),
                                     _hamming_reference(motif, txt, radius))

    def test_over_a_private_use_alphabet(self):
        cod = A.codon()
        txt = T.tokenize(["atgaaatttggg" * 5, "atgaaaccc" * 3], cod).text
        motif = cod.encode(["atg", "aaa", "ttt"])
        self.assertEqual(S.hamming_family_positions(motif, txt, 2, cod), _hamming_reference(motif, txt, 2))


class NormalizeCacheTests(unittest.TestCase):
    def test_a_long_text_is_checked_once_but_a_changed_text_is_checked_again(self):
        lat = A.latin()
        big = "abcdefghij" * 500
        self.assertEqual(S.normalize_sequence(big, lat), big)
        self.assertEqual(S.normalize_sequence(big, lat), big)            # remembered
        bad = "abcdefghij" * 499 + "abcdefghi"
        with self.assertRaises(ValueError):
            S.normalize_sequence(bad, lat)
        with self.assertRaises(ValueError):                              # not remembered as fine
            S.normalize_sequence(bad, lat)

    def test_the_alphabet_is_part_of_what_is_remembered(self):
        txt = "acgt" * 1000
        self.assertEqual(S.normalize_sequence(txt, A.nucleotide()), txt)
        small = A.Alphabet("ab", ("a", "b"))
        with self.assertRaises(ValueError):
            S.normalize_sequence(txt, small)

    def test_nucleotide_text_is_filtered_as_before(self):
        self.assertEqual(S.normalize_sequence("AUGX-c" * 1000), "atgc" * 1000)
        self.assertEqual(S.normalize_sequence("AUGC"), "atgc")


class FindAllMatchesTests(unittest.TestCase):
    def test_overlapping_matches_and_case(self):
        self.assertEqual(find_all_matches("aa", "AAAA"), [0, 1, 2])
        self.assertEqual(find_all_matches("a.c", "abcaxc"), [0, 3])


if __name__ == "__main__":
    unittest.main()

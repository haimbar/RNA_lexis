"""Statistical functions over alphabets other than nucleotides (phase 2)."""

import itertools
import math
import random
import unittest
import warnings

from scipy.stats import poisson

from rna_lexis.algorithms import find_boundary
from rna_lexis.alphabet import BOUNDARY, WILDCARD, Alphabet, codon, latin, nucleotide
from rna_lexis.statistical import (
    MarkovBackground, bh_adjust, default_markov_order, find_gapped_motif_hits,
    hamming_family_positions, hamming_variants, literal_positions, markov_pattern_probability,
    mutation_family_tests, normalize_sequence, rank_core_candidates, score_exact_motifs,
    score_gapped_motif, _variant_count,
)
from rna_lexis.tokenizers import tokenize

CODON = codon()
CODONS = CODON.symbols


def codon_text(n, seed, planted=()):
    """n random codons as an internal string; planted = [(codons, [positions])]."""
    rnd = random.Random(seed)
    units = [rnd.choice(CODONS) for _ in range(n)]
    for phrase, starts in planted:
        for s in starts:
            units[s:s + len(phrase)] = phrase
    return CODON.encode(units)


class NormalizeTests(unittest.TestCase):
    def test_nucleotide_default_is_unchanged(self):
        self.assertEqual(normalize_sequence("AUGC_xx ng"), "atgcg")

    def test_alphabet_keeps_the_boundary_and_checks_symbols(self):
        self.assertEqual(normalize_sequence("acg_tt", nucleotide()), "acg_tt")
        with self.assertRaises(ValueError):
            normalize_sequence("acgx", nucleotide())
        with self.assertRaises(ValueError):
            normalize_sequence("atg", CODON)            # external text, not an internal string

    def test_nucleotide_alphabet_reads_u_as_t(self):
        self.assertEqual(normalize_sequence("acgu", nucleotide()), "acgt")


class DefaultOrderTests(unittest.TestCase):
    def test_rule(self):
        self.assertEqual(default_markov_order(64, 300), 0)       # codons
        self.assertEqual(default_markov_order(64, 4095), 0)
        self.assertEqual(default_markov_order(64, 4096), 1)
        self.assertEqual(default_markov_order(26, 1000), 1)      # letters
        self.assertEqual(default_markov_order(26, 675), 0)
        self.assertEqual(default_markov_order(4, 16), 1)         # nucleotides

    def test_scoring_uses_it_when_no_order_is_given(self):
        txt = codon_text(300, 1)
        rows = score_exact_motifs(txt, [txt[:3]], alphabet=CODON)
        self.assertEqual(rows[0]["markov_order"], 0)
        rnd = random.Random(2)
        ltxt = "".join(rnd.choice("abcdefghijklmnopqrstuvwxyz") for _ in range(1500))
        rows = score_exact_motifs(ltxt, [ltxt[:4]], alphabet=latin())
        self.assertEqual(rows[0]["markov_order"], 1)
        rows = score_exact_motifs(ltxt, [ltxt[:4]], alphabet=latin(), markov_order=0)
        self.assertEqual(rows[0]["markov_order"], 0)


class MarkovBackgroundTests(unittest.TestCase):
    def test_declared_size_is_used_whether_or_not_symbols_occur(self):
        txt = CODON.encode(["atg", "aaa", "atg"])
        bg = MarkovBackground(txt, order=0, alphabet=CODON)
        self.assertEqual(len(bg.alphabet), 64)
        unseen = CODON.char("ttt")
        self.assertAlmostEqual(bg.probability(unseen), 1 / (3 + 64))
        self.assertAlmostEqual(bg.probability(CODON.char("atg")), 3 / (3 + 64))

    def test_nucleotide_alphabet_agrees_with_the_legacy_background(self):
        rnd = random.Random(3)
        txt = "".join(rnd.choice("acgt") for _ in range(300))
        for order in (0, 1, 2):
            old = MarkovBackground(txt, order=order)
            new = MarkovBackground(txt, order=order, alphabet=nucleotide())
            for w in ("acgt", "tt", "gattaca"):
                self.assertEqual(old.probability(w), new.probability(w))
                self.assertEqual(old.expected_count(w), new.expected_count(w))

    def test_boundaries_are_not_symbols_and_no_transition_crosses_one(self):
        a = nucleotide()
        one = MarkovBackground("aaaa_cccc", order=1, alphabet=a)
        # 8 letters, none of them the boundary; the a->c step across "_" is not counted
        self.assertEqual(one._n, 8)
        self.assertEqual(one._transition_counts.get(("a", "c"), 0), 0)
        self.assertEqual(one._transition_counts[("a", "a")], 3)
        self.assertEqual(one._transition_counts[("c", "c")], 3)
        # slots: a word of length 4 fits once in each segment
        self.assertEqual(one.slots(4), 2)
        self.assertEqual(one.slots(5), 0)
        self.assertEqual(one.probability("a_c"), 0.0)

    def test_without_a_boundary_slots_are_the_usual_ones(self):
        bg = MarkovBackground("acgtacgt", order=1, alphabet=nucleotide())
        self.assertEqual(bg.slots(3), 6)


class VariantTests(unittest.TestCase):
    def test_nucleotide_variants_are_all_distinct_and_within_radius(self):
        out = list(hamming_variants("acgta", 2))
        self.assertEqual(len(out), 1 + 5 * 3 + 10 * 9)
        self.assertEqual(len(set(out)), len(out))
        self.assertEqual(out[0], "acgta")
        for v in out:
            self.assertLessEqual(sum(x != y for x, y in zip(v, "acgta")), 2)
        self.assertEqual(out[1:5], ["ccgta", "gcgta", "tcgta", "aagta"])

    def test_variant_count_follows_the_alphabet_size(self):
        self.assertEqual(_variant_count(6, 1), 1 + 6 * 3)
        self.assertEqual(_variant_count(6, 1, 64), 1 + 6 * 63)
        self.assertEqual(_variant_count(6, 2, 64), 1 + 6 * 63 + 15 * 63 * 63)

    def test_codon_variants(self):
        motif = CODON.encode(["atg", "aaa", "ttt"])
        out = list(hamming_variants(motif, 1, alphabet=CODON))
        self.assertEqual(len(out), 1 + 3 * 63)
        self.assertEqual(len(set(out)), len(out))
        with self.assertRaises(ValueError):
            list(hamming_variants(CODON.encode(["atg"] * 6), 3, alphabet=CODON))   # too many variants

    def test_a_boundary_in_a_motif_is_rejected(self):
        with self.assertRaises(ValueError):
            list(hamming_variants("ac_gt", 1, alphabet=nucleotide()))


class MultipleTestingTests(unittest.TestCase):
    def test_bh_with_the_number_of_tests_given(self):
        p = [0.04, 0.01, 0.03]
        self.assertEqual(bh_adjust(p, m=3), bh_adjust(p))
        # two p-values out of m = 100 tests, the rest taken as p = 1
        q = bh_adjust([0.001, 0.01], m=100)
        self.assertAlmostEqual(q[0], 0.1)
        self.assertAlmostEqual(q[1], 0.5)
        with self.assertRaises(ValueError):
            bh_adjust([0.1, 0.2], m=1)

    def test_q_value_is_corrected_over_all_possible_words_of_that_length(self):
        txt = codon_text(300, 4, planted=[(CODONS[:4], [10, 100, 200])])
        motif = CODON.encode(CODONS[:4])
        row = score_exact_motifs(txt, [motif], alphabet=CODON)[0]
        self.assertEqual(row["exact_count"], 3)
        self.assertAlmostEqual(row["q_markov"], min(1.0, row["p_markov"] * 64 ** 4), places=12)
        # corrected over the candidate list only, as nucleotide mode does, it is far smaller
        cand = score_exact_motifs(txt, [motif], alphabet=CODON, family="candidates")[0]
        self.assertAlmostEqual(cand["q_markov"], row["p_markov"])
        self.assertTrue(cand["statistically_supported"])

    def test_a_chance_repeat_is_not_called_supported_by_length(self):
        # a codon pair that happens to occur twice is not significant against 64**2 possible pairs
        txt = codon_text(300, 5)
        pairs = {}
        for i in range(len(txt) - 1):
            pairs.setdefault(txt[i:i + 2], []).append(i)
        repeated = [w for w, pos in pairs.items() if len(pos) >= 2]
        self.assertTrue(repeated)
        rows = score_exact_motifs(txt, repeated, alphabet=CODON)
        self.assertFalse(any(r["statistically_supported"] for r in rows))

    def test_default_family_depends_on_the_alphabet(self):
        txt = codon_text(300, 4, planted=[(CODONS[:4], [10, 100, 200])])
        motif = CODON.encode(CODONS[:4])
        default = score_exact_motifs(txt, [motif], alphabet=CODON)[0]["q_markov"]
        bylength = score_exact_motifs(txt, [motif], alphabet=CODON, family="bylength")[0]["q_markov"]
        self.assertEqual(default, bylength)
        with self.assertRaises(ValueError):
            score_exact_motifs(txt, [motif], alphabet=CODON, family="nope")


class ScoringMatchesTheCalibrationTests(unittest.TestCase):
    def test_expected_count_is_the_order0_flat_pseudocount_model(self):
        txt = codon_text(500, 6)
        motif = txt[40:43]
        row = score_exact_motifs(txt, [motif], alphabet=CODON)[0]
        n = len(txt)
        counts = {c: txt.count(c) for c in set(txt)}
        p = 1.0
        for c in motif:
            p *= (counts.get(c, 0) + 1) / (n + 64)
        self.assertAlmostEqual(row["expected_markov"], (n - 3 + 1) * p, places=12)
        self.assertAlmostEqual(row["p_markov"], float(poisson.sf(row["exact_count"] - 1, row["expected_markov"])))


class SegmentTests(unittest.TestCase):
    def test_motifs_do_not_span_a_boundary(self):
        tok = tokenize(["acgt", "acgt"], nucleotide())
        self.assertEqual(tok.text, "acgt_acgt")
        self.assertEqual(literal_positions("tacg", tok.text, nucleotide()), [])
        self.assertEqual(literal_positions("acgt", tok.text, nucleotide()), [0, 5])
        with self.assertRaises(ValueError):
            literal_positions("t_a", tok.text, nucleotide())

    def test_hamming_windows_and_gapped_hits_skip_boundaries(self):
        a = nucleotide()
        hits = hamming_family_positions("acgt", "acgt_acgt", 1, a)
        self.assertEqual([h[0] for h in hits], [0, 5])
        gapped = find_gapped_motif_hits("aa_cc", "a", "c", min_gap=0, max_gap=3, alphabet=a)
        self.assertEqual(gapped, [])
        gapped = find_gapped_motif_hits("aacc", "a", "c", min_gap=0, max_gap=3, alphabet=a)
        self.assertTrue(gapped)


class GappedMotifTests(unittest.TestCase):
    def test_wildcard_probability_equals_the_sum_over_every_fill(self):
        abc = Alphabet(name="abc", symbols=("a", "b", "c"))
        rnd = random.Random(7)
        txt = "".join(rnd.choice("abc") for _ in range(200))
        for order in (0, 1, 2):
            bg = MarkovBackground(txt, order=order, alphabet=abc)
            pattern = "a" + WILDCARD * 2 + "c"
            brute = sum(bg.probability("a" + x + y + "c") for x in "abc" for y in "abc")
            self.assertAlmostEqual(markov_pattern_probability(pattern, bg), brute, places=12)

    def test_letters_n_and_dot_are_ordinary_symbols(self):
        # in the nucleotide code "n" is a wildcard; over Latin letters it is a letter
        lat = latin()
        txt = "banana" * 20
        bg = MarkovBackground(txt, order=0, alphabet=lat)
        self.assertAlmostEqual(markov_pattern_probability("nan", bg), bg.probability("nan"))
        hits = find_gapped_motif_hits("banana", "b", "n", min_gap=1, max_gap=1, alphabet=lat)
        self.assertEqual([h["gap_length"] for h in hits], [1])
        with self.assertRaises(ValueError):
            markov_pattern_probability("a?c", bg)                # "?" is neither a letter nor the wildcard

    def test_score_gapped_motif_over_codons(self):
        left, right = CODON.encode(CODONS[:2]), CODON.encode(CODONS[5:7])
        gap = CODON.encode(CODONS[20:23])
        txt = codon_text(300, 8)
        txt = txt[:50] + left + gap + right + txt[50 + 2 + 3 + 2:]
        res = score_gapped_motif(txt, left, right, min_gap=0, max_gap=4, alphabet=CODON)
        self.assertGreaterEqual(res["observed_count"], 1)
        self.assertEqual(res["markov_order"], 0)
        self.assertGreater(res["expected_markov"], 0)


class MutationFamilyTests(unittest.TestCase):
    def test_family_size_uses_the_declared_alphabet(self):
        motif = CODON.encode(CODONS[:6])
        txt = codon_text(400, 9, planted=[(CODONS[:6], [10, 120])])
        txt = txt[:250] + motif[:3] + CODON.char(CODONS[40]) + motif[4:] + txt[256:]   # one-codon variant
        rows = mutation_family_tests(txt, [motif], alphabet=CODON, mutr=1 / 3)
        row = next(r for r in rows if r["radius"] == 1)
        self.assertEqual(row["hamming_variant_count"], 1 + 6 * 63)
        self.assertEqual(row["exact_count"], 2)
        self.assertEqual(row["approx_count"], 1)


class PipelineTests(unittest.TestCase):
    def test_a_planted_core_family_is_found_end_to_end(self):
        core = CODONS[:5]
        ext1, ext2 = CODONS[30:32], CODONS[40:42]
        rnd = random.Random(11)
        units = [rnd.choice(CODONS) for _ in range(300)]
        spots = [(10, ext1), (70, ext1), (140, ext2), (210, ext2)]
        for s, ext in spots:
            units[s:s + len(core) + 2] = list(core) + list(ext)
        txt = CODON.encode(units)
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()):
            xm = sorted(set(find_boundary(txt, minlen=4, maxlen=60, at_least=2)), key=lambda x: (-len(x), x))
        rows = rank_core_candidates(txt, xm, candidate_min_len=3, alphabet=CODON)
        core_s = CODON.encode(core)
        supported = [r["motif"].lower() for r in rows if r["statistically_supported"]]
        self.assertTrue(any(core_s in m or m in core_s for m in supported), supported)
        # no chance repeats: everything called supported overlaps the planted material
        planted = set(core) | set(ext1) | set(ext2)
        for m in supported:
            self.assertTrue(set(CODON.decode_list(m)) & planted)


if __name__ == "__main__":
    unittest.main()

import unittest

import numpy as np

from rna_lexis.algorithms import (
    allow_mutation,
    compute_default_wd,
    contains_only_rna,
    cores,
    count_kgrams,
    cover,
    expand_to_boundary,
    extend_match_pair,
    extendRNA,
    find_all_matches,
    find_boundary,
    find_longest_extensions,
    find_with_mutations,
    gen_hairpins,
    hop_distance,
    is_bounded,
    markov_kmer_pvalues,
    zscore,
)


class ContainsOnlyRnaTests(unittest.TestCase):
    def test_pure_dna_accepted(self):
        self.assertTrue(contains_only_rna("acgtacgt"))

    def test_pure_rna_accepted(self):
        self.assertTrue(contains_only_rna("acguacgu"))

    def test_rejects_other_letters(self):
        self.assertFalse(contains_only_rna("acgtx"))


class CountKgramsTests(unittest.TestCase):
    def test_counts_overlapping_kmers(self):
        counts = count_kgrams("aaaa", 2, at_least=1)
        self.assertEqual(counts["aa"], 3)


class FindAllMatchesTests(unittest.TestCase):
    def test_finds_all_positions(self):
        self.assertEqual(find_all_matches("ata", "atatata"), [0, 2, 4])

    def test_wildcard_dot(self):
        matches = find_all_matches("a.a", "axaaya", ret="str")
        self.assertEqual(matches, ["axa", "aya"])


class BoundaryAndCoresTests(unittest.TestCase):
    def test_find_boundary_finds_repeated_flanked_motif(self):
        # "tgtatata" occurs 3x, each time with a different left/right flank
        # (c/a/g on the left, a/g/t on the right), satisfying is_bounded().
        txt = "ctgtatataatgtatatagtgtatatat"
        found = find_boundary(txt, minlen=8, maxlen=8, at_least=3)
        self.assertIn("tgtatata", found)

    def test_cores_finds_substring_shared_by_noncontaining_xmotifs(self):
        txt = "aaatgtatataccctgtatataaaatgtatatagggg"
        xmotifs = ["aaatgtatataccc", "gggtgtatataaaa", "aaatgtatatagggg"]
        found = cores(txt, xmotifs, minclen=6)
        self.assertTrue(any("tgtatata" in c for c in found))

    def test_cores_empty_without_shared_substring(self):
        self.assertEqual(cores("", [], minclen=6), [])

    def test_is_bounded_true_for_maximal_word(self):
        txt = "ccc" + "gattaca" + "nnn" + "gattaca" + "zzz"
        self.assertTrue(is_bounded("gattaca", txt))

    def test_expand_to_boundary_grows_until_maximal(self):
        txt = "ccc" + "gattaca" + "nnn" + "gattaca" + "zzz"
        self.assertEqual(expand_to_boundary("attac", txt, dir="b"), "gattaca")


class MutationSearchTests(unittest.TestCase):
    def test_find_with_mutations_separates_exact_and_approximate(self):
        txt = "aaatgcatgcgggtgcaggcttt"
        exact, approx, maxmut = find_with_mutations("tgcatgc", txt, mutr=1 / 6)
        self.assertEqual(exact, [3])
        self.assertEqual(len(approx), 1)
        self.assertEqual(approx[0][2], 1)  # hamming distance
        self.assertEqual(maxmut, 1)

    def test_allow_mutation_returns_unchanged_for_short_strings(self):
        # len <= 6 is returned as-is (no mutation slots available)
        self.assertEqual(allow_mutation("abcdef", n=1), "abcdef")

    def test_allow_mutation_builds_one_variant_per_position(self):
        pattern = allow_mutation("gattacagattaca", n=1)
        self.assertIn("[^g]attacagattaca", pattern)
        self.assertEqual(pattern.count("|"), len("gattacagattaca") - 1)


class ExtensionTests(unittest.TestCase):
    def setUp(self):
        self.txt = "ccc" + "gattaca" + "nnn" + "gattaca" + "zzz"

    def test_find_longest_extensions_extends_past_seed(self):
        results = find_longest_extensions("gattaca", self.txt, mutr=1 / 6)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual((r["pos1"], r["pos2"]), (3, 13))
        self.assertEqual(r["total_len"], 8)
        self.assertEqual(r["hamming"], 1)

    def test_extend_match_pair_agrees_with_find_longest_extensions(self):
        left, right, ext1, ext2, ham = extend_match_pair(
            "gattaca", self.txt, 3, 13, mutr=1 / 6
        )
        self.assertEqual((left, right), (0, 1))
        self.assertEqual(ham, 1)

    def test_find_longest_extensions_empty_below_two_occurrences(self):
        self.assertEqual(find_longest_extensions("gattaca", "gattaca"), [])


class CoverageAndWidthTests(unittest.TestCase):
    def test_cover_positive_for_nonoverlapping_occurrences(self):
        txt = "ccc" + "gattaca" + "nnn" + "gattaca" + "zzz"
        self.assertGreater(cover("gattaca", txt), 0)

    def test_cover_zero_for_single_occurrence(self):
        self.assertEqual(cover("gattaca", "xxxgattacayyy"), 0)

    def test_compute_default_wd_within_configured_bounds(self):
        txt = "ccc" + "gattaca" + "nnn" + "gattaca" + "zzz"
        wd = compute_default_wd(txt, ["gattaca"])
        self.assertGreaterEqual(wd, 20)
        self.assertLessEqual(wd, 80)

    def test_compute_default_wd_falls_back_to_max_when_no_cores(self):
        self.assertEqual(compute_default_wd("acgtacgt", []), 80)


class HairpinTests(unittest.TestCase):
    def test_gen_hairpins_finds_self_complementary_stem_loop(self):
        stem1 = "gcgcgcgcg"
        stem2 = extendRNA(stem1, xtype=3)  # reverse complement
        rna = stem1 + "aaaa" + stem2
        hairpins = gen_hairpins(rna, minSL=7, minLL=3, maxLL=40)
        self.assertEqual(len(hairpins), 1)
        hp = hairpins[0]
        self.assertEqual(hp.stem_len, 9)
        self.assertEqual(hp.loop_len, 4)
        self.assertEqual(hp.stem_seq, stem1)
        self.assertEqual(hp.loop_seq, "aaaa")

    def test_gen_hairpins_empty_when_no_stem(self):
        self.assertEqual(gen_hairpins("a" * 30, minSL=8), [])


class ZscoreTests(unittest.TestCase):
    def test_zscore_centers_around_zero(self):
        z = zscore(np.array([1, 2, 3, 4, 5]), robust=False)
        self.assertAlmostEqual(float(np.mean(z)), 0.0, places=6)


class MarkovKmerPvaluesTests(unittest.TestCase):
    def test_enriched_kmer_gets_low_over_pvalue(self):
        txt = "gattaca" * 10 + "ttggccaaatgcgcgcgcatgcatatatgcgc" * 5
        rows = markov_kmer_pvalues(txt, k=7, order=1)
        by_kmer = {r["kmer"]: r for r in rows}
        self.assertLess(by_kmer["gattaca"]["pvalue_over"], 0.01)
        self.assertEqual(by_kmer["gattaca"]["direction"], "over")


class HopDistanceTests(unittest.TestCase):
    def test_exact_multiples_of_unit(self):
        self.assertEqual(hop_distance(290, 290), 1)
        self.assertEqual(hop_distance(580, 290), 2)
        self.assertEqual(hop_distance(870, 290), 3)

    def test_rounds_to_nearest_unit(self):
        # 592 / 290 = 2.04 -> rounds to 2, not 3
        self.assertEqual(hop_distance(592, 290), 2)
        # 875 / 290 = 3.02 -> rounds to 3
        self.assertEqual(hop_distance(875, 290), 3)

    def test_never_returns_zero(self):
        # A spacing much smaller than the unit still counts as >= 1 hop --
        # two distinct positions are never "0 units apart".
        self.assertEqual(hop_distance(10, 290), 1)
        self.assertEqual(hop_distance(0, 290), 1)

    def test_unit_is_not_hardcoded(self):
        # Same physical spacing, different repeat-unit length -> different
        # hop bucket. Guards against ever reintroducing a fixed default.
        self.assertEqual(hop_distance(600, 300), 2)
        self.assertEqual(hop_distance(600, 600), 1)


if __name__ == "__main__":
    unittest.main()

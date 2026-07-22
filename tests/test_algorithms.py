import unittest

import numpy as np

from rna_lexis.algorithms import (
    allow_mutation,
    binary_coverage,
    compute_default_wd,
    contains_only_rna,
    coverage_comparison_warnings,
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
    gc_content_indicator,
    gen_hairpins,
    homopolymer_indicator,
    hop_distance,
    is_bounded,
    markov_kmer_pvalues,
    sliding_coverage_fraction,
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


class BinaryCoverageTests(unittest.TestCase):
    def test_marks_full_span_of_each_occurrence(self):
        # "gattaca" occurs at 3-10 and 13-20 in this text.
        txt = "ccc" + "gattaca" + "nnn" + "gattaca" + "zzz"
        cov = binary_coverage(txt, ["gattaca"])
        self.assertEqual(list(cov[:3]), [0, 0, 0])
        self.assertEqual(list(cov[3:10]), [1] * 7)
        self.assertEqual(list(cov[10:13]), [0, 0, 0])
        self.assertEqual(list(cov[13:20]), [1] * 7)

    def test_overlapping_motifs_union_without_double_counting(self):
        cov = binary_coverage("gattacagattaca", ["gattaca"])
        self.assertEqual(list(cov), [1] * 14)
        self.assertTrue(set(cov.tolist()) <= {0, 1})

    def test_no_motifs_gives_all_zero(self):
        cov = binary_coverage("acgtacgt", [])
        self.assertEqual(list(cov), [0] * 8)


class GcContentIndicatorTests(unittest.TestCase):
    def test_marks_g_and_c_only(self):
        self.assertEqual(list(gc_content_indicator("ggccaaat")), [1, 1, 1, 1, 0, 0, 0, 0])

    def test_case_insensitive(self):
        self.assertEqual(list(gc_content_indicator("GgCcAaTt")), [1, 1, 1, 1, 0, 0, 0, 0])


class HomopolymerIndicatorTests(unittest.TestCase):
    def test_marks_runs_at_or_above_min_run(self):
        # aaaa(4)=yes, tt(2)=no, gggg(4)=yes, c(1)=no
        cov = homopolymer_indicator("aaaattggggc", min_run=4)
        self.assertEqual(list(cov), [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0])

    def test_min_run_is_configurable(self):
        cov = homopolymer_indicator("aattggc", min_run=2)
        self.assertEqual(list(cov), [1, 1, 1, 1, 1, 1, 0])

    def test_no_run_meets_threshold(self):
        cov = homopolymer_indicator("acgtacgt", min_run=4)
        self.assertEqual(list(cov), [0] * 8)


class SlidingCoverageFractionTests(unittest.TestCase):
    def test_uniform_coverage_stays_uniform(self):
        cov = np.ones(100, dtype=np.int8)
        curve = sliding_coverage_fraction(cov, window=10, n_points=20)
        self.assertEqual(len(curve), 20)
        self.assertTrue(np.allclose(curve, 1.0))

    def test_no_coverage_stays_zero(self):
        cov = np.zeros(100, dtype=np.int8)
        curve = sliding_coverage_fraction(cov, window=10, n_points=20)
        self.assertTrue(np.allclose(curve, 0.0))

    def test_window_larger_than_input_does_not_crash(self):
        cov = np.array([1, 0, 1, 0], dtype=np.int8)
        curve = sliding_coverage_fraction(cov, window=1000, n_points=10)
        self.assertEqual(len(curve), 10)
        self.assertTrue(np.all(curve >= 0) and np.all(curve <= 1))


class CoverageComparisonWarningsTests(unittest.TestCase):
    def test_no_warnings_for_similar_lengths(self):
        # Real reference case: bundled human vs mouse NORAD (5401 vs 4945 nt).
        self.assertEqual(coverage_comparison_warnings(5401, 4945, window=200), [])

    def test_warns_on_large_length_ratio(self):
        warnings = coverage_comparison_warnings(5401, 2000, window=200)
        self.assertEqual(len(warnings), 1)
        self.assertIn("differ substantially", warnings[0])

    def test_warns_on_short_sequences_regardless_of_ratio(self):
        # Equal lengths (ratio 1.0) but both short relative to the window.
        warnings = coverage_comparison_warnings(300, 300, window=200)
        self.assertEqual(len(warnings), 1)
        self.assertIn("short relative to the window", warnings[0])

    def test_can_return_both_warnings_at_once(self):
        warnings = coverage_comparison_warnings(1000, 100, window=200)
        self.assertEqual(len(warnings), 2)

    def test_ratio_threshold_is_configurable(self):
        # 1000 vs 900 -> ratio ~1.111: passes a looser threshold, fails a tighter one.
        self.assertEqual(coverage_comparison_warnings(1000, 900, window=1, ratio_threshold=1.2), [])
        warnings = coverage_comparison_warnings(1000, 900, window=1, ratio_threshold=1.05)
        self.assertTrue(any("differ substantially" in w for w in warnings))


if __name__ == "__main__":
    unittest.main()

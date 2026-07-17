import unittest

from rna_lexis.statistical import (
    find_gapped_motif_hits,
    hamming_family_positions,
    score_gapped_motif,
    mutation_family_tests,
    normalize_sequence,
    rank_core_candidates,
    score_exact_motifs,
    shared_exact_motifs,
)


class StatisticalWorkflowTests(unittest.TestCase):
    def test_normalize_sequence_accepts_rna_and_dna(self):
        self.assertEqual(normalize_sequence("AAUUxxCCgg"), "aattccgg")

    def test_hamming_family_positions_counts_exact_and_approximate_windows(self):
        hits = hamming_family_positions("TGTATATA", "TGTATATAAACGTATATA", 1)
        distances = sorted(d for _, _, d in hits)
        self.assertEqual(distances, [0, 1])

    def test_exact_scoring_returns_markov_support_columns(self):
        txt = "TGTATATA" * 5 + "ACCGGTTACGACCGT"
        rows = score_exact_motifs(txt, ["TGTATATA"], enrichment_threshold=1.0)
        self.assertEqual(rows[0]["motif"], "TGTATATA")
        self.assertIn("expected_markov", rows[0])
        self.assertIn("q_markov", rows[0])

    def test_core_ranking_scores_shared_xmotif_substrings(self):
        txt = "AAATGTATATACCCGGGTGTATATAAAATGTATATAGGG"
        xmotifs = ["AAATGTATATACCC", "GGGTGTATATAAAA", "AAATGTATATAGGG"]
        rows = rank_core_candidates(
            txt,
            xmotifs,
            candidate_min_len=5,
            candidate_max_len=10,
            enrichment_threshold=1.0,
        )
        motifs = {row["motif"] for row in rows}
        self.assertIn("TGTATATA", motifs)

    def test_mutation_family_tests_include_decision_column(self):
        txt = "TGTATATACCCCTGTATATAGGGGTGCATATAAAAATGTATATA"
        rows = mutation_family_tests(txt, ["TGTATATA"], enrichment_threshold=1.0)
        self.assertTrue(rows)
        self.assertIn("decision", rows[0])
        self.assertTrue(any(row["radius"] == 1 for row in rows))

    def test_gapped_motif_search_and_score(self):
        txt = "AAATGTCCCTATAAAATGTGGGGTATA"
        hits = find_gapped_motif_hits(txt, "TGT", "TATA", min_gap=3, max_gap=4)
        self.assertEqual(len(hits), 2)
        score = score_gapped_motif(txt, "TGT", "TATA", min_gap=3, max_gap=4)
        self.assertEqual(score["observed_count"], 2)
        self.assertIn("expected_markov", score)


class SharedExactMotifsTests(unittest.TestCase):
    def setUp(self):
        # Two differently-flanked 11-nt variants of the core "gattaca",
        # each repeated with distinct outer flanks -- this is what makes
        # cores() derive "gattaca" as a shared substring of two
        # non-containing xmotifs (see algorithms.cores()'s docstring), the
        # same mechanism the joint concat-discovery step relies on.
        core = "gattaca"
        v1 = "xx" + core + "yy"
        v2 = "pp" + core + "qq"
        self.txt_a = "a" + v1 + "b" + "c" + v1 + "d" + "e" + v2 + "f" + "g" + v2 + "h"

    def test_finds_planted_shared_motif_via_joint_discovery(self):
        txt_b = "zzz" + "gattaca" + "www"
        shared = shared_exact_motifs(self.txt_a, txt_b)
        self.assertIn("gattaca", shared)

    def test_excludes_motif_absent_from_txt_b(self):
        txt_b = "no matching core here"
        shared = shared_exact_motifs(self.txt_a, txt_b)
        self.assertNotIn("gattaca", shared)

    def test_cores_a_supplement_adds_motif_missed_by_joint_scan(self):
        # "uniquemotif" occurs only once in each sequence -- below the
        # combined-repeat threshold find_boundary() needs to find it on its
        # own -- but is picked up via the cores_a supplement path.
        txt_a = "z" * 20 + "uniquemotif" + "z" * 20
        txt_b = "w" * 10 + "uniquemotif" + "w" * 10
        shared = shared_exact_motifs(txt_a, txt_b, cores_a=["uniquemotif"])
        self.assertIn("uniquemotif", shared)

    def test_no_supplement_effect_without_cores_a(self):
        txt_a = "z" * 20 + "uniquemotif" + "z" * 20
        txt_b = "w" * 10 + "uniquemotif" + "w" * 10
        shared = shared_exact_motifs(txt_a, txt_b)
        self.assertNotIn("uniquemotif", shared)

    def test_min_len_filters_short_candidates(self):
        shared = shared_exact_motifs("xxxxxxxx", "yyyyyyyy", cores_a=["ab"], min_len=6)
        self.assertEqual(shared, [])

    def test_n_top_caps_result_count(self):
        cores_a = ["motifone", "motiftwo", "motifthre"]
        txt_a = "".join(cores_a)
        txt_b = "".join(cores_a)
        shared = shared_exact_motifs(txt_a, txt_b, cores_a=cores_a, n_top=2)
        self.assertEqual(len(shared), 2)


if __name__ == "__main__":
    unittest.main()

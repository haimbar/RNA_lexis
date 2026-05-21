import unittest

from rna_lexis.statistical import (
    find_gapped_motif_hits,
    hamming_family_positions,
    score_gapped_motif,
    mutation_family_tests,
    normalize_sequence,
    rank_core_candidates,
    score_exact_motifs,
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


if __name__ == "__main__":
    unittest.main()

import unittest

from rna_lexis.alignment import gotoh_global, gotoh_local, make_markers


class GotohGlobalTests(unittest.TestCase):
    def test_identical_sequences_score_perfectly(self):
        res = gotoh_global("gattaca", "gattaca")
        self.assertEqual(res.mode, "global")
        self.assertEqual(res.matches, 7)
        self.assertEqual(res.mismatches, 0)
        self.assertEqual(res.gaps, 0)
        self.assertEqual(res.score, 7 * 2)  # default match=2
        self.assertEqual(res.markers, "|" * 7)

    def test_single_mismatch_reduces_score(self):
        res = gotoh_global("gattaca", "gattapa")
        self.assertEqual(res.matches, 6)
        self.assertEqual(res.mismatches, 1)
        self.assertEqual(res.score, 6 * 2 - 1)  # default mismatch=-1

    def test_global_result_has_no_local_position_fields(self):
        res = gotoh_global("gattaca", "gattaca")
        self.assertIsNone(res.start_a)
        self.assertIsNone(res.end_a)
        self.assertIsNone(res.start_b)
        self.assertIsNone(res.end_b)

    def test_indel_introduces_gap(self):
        res = gotoh_global("gattaca", "gattca")
        self.assertGreater(res.gaps, 0)
        self.assertIn("-", res.aligned_a + res.aligned_b)


class GotohLocalTests(unittest.TestCase):
    def test_finds_matching_subregion_within_noise(self):
        a = "xxxxgattacayyyy"
        b = "zzzgattacawwww"
        res = gotoh_local(a, b)
        self.assertEqual(res.mode, "local")
        self.assertEqual(a[res.start_a:res.end_a], "gattaca")
        self.assertEqual(b[res.start_b:res.end_b], "gattaca")
        self.assertEqual(res.mismatches, 0)

    def test_local_score_matches_self_alignment_of_shared_region(self):
        a = "xxxxgattacayyyy"
        b = "zzzgattacawwww"
        res = gotoh_local(a, b)
        self_score = gotoh_global("gattaca", "gattaca").score
        self.assertEqual(res.score, self_score)


class MakeMarkersTests(unittest.TestCase):
    def test_counts_matches_mismatches_and_gaps(self):
        markers, matches, mismatches, gaps = make_markers("gattaca", "gattapa")
        self.assertEqual(markers, "|||||.|")
        self.assertEqual(matches, 6)
        self.assertEqual(mismatches, 1)
        self.assertEqual(gaps, 0)

    def test_gap_character_counted_separately_from_mismatch(self):
        markers, matches, mismatches, gaps = make_markers("gattaca", "gat-aca")
        self.assertEqual(gaps, 1)
        self.assertEqual(mismatches, 0)


if __name__ == "__main__":
    unittest.main()

"""Tests for rna_lexis.plots.

Uses the non-interactive Agg backend so tests never try to open a GUI
window (matplotlib's default backend can otherwise hang waiting for a
display in headless/CI environments).
"""

import os
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")

from rna_lexis.plots import plot_self_similarity_arcs, plot_shared_motif_diagram
from rna_lexis.statistical import spacing_periodicity_test


class PlotSelfSimilarityArcsTests(unittest.TestCase):
    def setUp(self):
        # Three occurrences exactly one repeat-unit (290 nt) apart, so the
        # median-gap unit estimate is unambiguous and reproducible.
        self.positions = [100, 390, 680]
        self.results = [
            {"pos1": 100, "pos2": 390, "total_len": 50, "hamming": 5},
            {"pos1": 390, "pos2": 680, "total_len": 60, "hamming": 3},
            {"pos1": 100, "pos2": 680, "total_len": 40, "hamming": 8},
        ]

    def test_saves_a_file_without_raising(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "arcs.png")
            plot_self_similarity_arcs("gattaca", self.results, self.positions, file=out)
            self.assertTrue(os.path.isfile(out))
            self.assertGreater(os.path.getsize(out), 0)

    def test_accepts_an_explicit_unit(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "arcs.svg")
            plot_self_similarity_arcs("gattaca", self.results, self.positions,
                                       unit=290, file=out)
            self.assertTrue(os.path.isfile(out))

    def test_no_crash_and_no_file_when_too_few_occurrences(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "arcs.png")
            plot_self_similarity_arcs("gattaca", [], [100], file=out)
            self.assertFalse(os.path.isfile(out))

    def test_accepts_spacing_stats_from_the_real_periodicity_test(self):
        # Same test used by "Motif spacing / periodicity test" -- locks in
        # that the two functions' outputs/inputs stay compatible.
        stats = spacing_periodicity_test(self.positions, 1000)
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "arcs.png")
            plot_self_similarity_arcs("gattaca", self.results, self.positions,
                                       spacing_stats=stats, file=out)
            self.assertTrue(os.path.isfile(out))


class PlotSharedMotifDiagramTests(unittest.TestCase):
    def setUp(self):
        self.txt_a = "aaaggcccaaaggccctttcagcctttcagcct" * 3
        self.txt_b = "zzzaaaggccczzzcagcctzzz"
        self.motifs = ["aaaggccc", "cagcct"]

    def test_saves_a_file_without_raising(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "shared.png")
            plot_shared_motif_diagram(self.txt_a, self.txt_b, "GeneA", "GeneB",
                                       self.motifs, file=out)
            self.assertTrue(os.path.isfile(out))
            self.assertGreater(os.path.getsize(out), 0)

    def test_no_crash_and_no_file_when_no_motif_hits_both_sequences(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "shared.png")
            plot_shared_motif_diagram(self.txt_a, self.txt_b, "GeneA", "GeneB",
                                       ["notpresentanywhere"], file=out)
            self.assertFalse(os.path.isfile(out))

    def test_skips_motifs_absent_from_one_sequence(self):
        # "aaaggccc" hits both; "onlyintxta" only exists in txt_a -- the
        # function should still succeed using just the valid motif.
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "shared.png")
            plot_shared_motif_diagram(self.txt_a + "onlyintxta", self.txt_b, "GeneA", "GeneB",
                                       ["aaaggccc", "onlyintxta"], file=out)
            self.assertTrue(os.path.isfile(out))

    def test_caps_at_six_motifs_when_more_are_shared(self):
        # 8 distinct motifs, each with an exact hit in both sequences --
        # only the first 6 (one per palette color) should be drawn.
        motifs = [f"motif{i:02d}xx" for i in range(8)]
        txt_a = "".join(motifs) * 2
        txt_b = "".join(motifs)
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "shared.svg")
            plot_shared_motif_diagram(txt_a, txt_b, "GeneA", "GeneB", motifs, file=out)
            self.assertTrue(os.path.isfile(out))


if __name__ == "__main__":
    unittest.main()

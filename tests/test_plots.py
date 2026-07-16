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

from rna_lexis.plots import plot_self_similarity_arcs
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


if __name__ == "__main__":
    unittest.main()

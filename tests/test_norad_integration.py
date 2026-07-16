"""Integration tests using the bundled NORAD example dataset (T.2).

Runs the real parsedata() discovery pipeline (find_boundary + cores) end to
end on real biological sequences, and asserts against NORAD's published,
peer-reviewed biology -- the tandem Pumilio Response Element (PRE) repeats
-- rather than exact golden-file counts, so the test stays valid across
unrelated future algorithm tuning. See src/rna_lexis/data/README.md for
dataset provenance.
"""

import os
import unittest

from rna_lexis.algorithms import find_all_matches, find_with_mutations
from rna_lexis.io import example_dataset_path
from rna_lexis.menu import _load_fasta_or_text, parsedata


def _load(name):
    return _load_fasta_or_text(example_dataset_path(name))["txt"]


class ExampleDatasetPathTests(unittest.TestCase):
    def test_human_and_mouse_paths_resolve_to_real_files(self):
        for name in ("NORAD_human", "NORAD_mouse"):
            self.assertTrue(os.path.isfile(example_dataset_path(name)))

    def test_unknown_dataset_name_raises(self):
        with self.assertRaises(ValueError):
            example_dataset_path("not_a_real_dataset")


class NoradHumanPreRepeatTests(unittest.TestCase):
    """NORAD's defining feature is its tandem Pumilio Response Element (PRE)
    repeats (consensus UGUANAUA / TGTANATA in cDNA form). These counts come
    from a literal substring/Hamming search, not a tunable discovery
    heuristic, so -- unlike the discovery-pipeline tests below -- they're
    asserted exactly.
    """

    @classmethod
    def setUpClass(cls):
        cls.txt = _load("NORAD_human")

    def test_sequence_length_matches_ensembl_transcript(self):
        self.assertEqual(len(self.txt), 5401)

    def test_core_pre_motif_occurs_exactly_ten_times(self):
        self.assertEqual(len(find_all_matches("tgtatata", self.txt)), 10)

    def test_pre_variant_found_under_one_mismatch_search(self):
        _exact, approx, _maxmut = find_with_mutations("tgtatata", self.txt, mutr=1 / 6)
        variant_hits = [matched for _pos, matched, _dist in approx if matched == "tgtaaata"]
        self.assertEqual(len(variant_hits), 4)


class NoradDiscoveryPipelineTests(unittest.TestCase):
    """Runs the real parsedata() pipeline end to end and checks that the
    known PRE motif is (a) discovered and (b) discovery finds a
    biologically-plausible number of xmotifs/cores. Bounds are loose by
    design -- exact counts (163 xmotifs / 235 cores for human, with default
    settings, matching the published RNA-Lexis paper's validation numbers
    for this exact transcript) would make this test brittle against
    unrelated future tuning of find_boundary()/cores().
    """

    @staticmethod
    def _discover(name):
        defvals = {"minxmlen": 7, "maxxmlen": 60, "mincorelen": 6, "mincount": 2}
        return parsedata(_load(name), defvals)

    def test_human_norad_discovers_the_pre_motif_among_cores(self):
        parsed = self._discover("NORAD_human")
        self.assertIn("tgtatata", parsed["corelist"])
        self.assertGreater(len(parsed["xmotifs"]), 50)
        self.assertGreater(len(parsed["corelist"]), 50)

    def test_mouse_norad_discovers_the_pre_motif_among_cores(self):
        parsed = self._discover("NORAD_mouse")
        self.assertIn("tgtatata", parsed["corelist"])
        self.assertGreater(len(parsed["xmotifs"]), 20)
        self.assertGreater(len(parsed["corelist"]), 20)

    def test_human_shows_more_organizational_complexity_than_mouse(self):
        # Matches the paper's finding that human NORAD is more densely
        # organized into recurrent motifs than its mouse ortholog.
        human = self._discover("NORAD_human")
        mouse = self._discover("NORAD_mouse")
        self.assertGreater(len(human["xmotifs"]), len(mouse["xmotifs"]))
        self.assertGreater(len(human["corelist"]), len(mouse["corelist"]))


if __name__ == "__main__":
    unittest.main()

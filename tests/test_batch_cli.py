"""End-to-end tests for rna_lexis.test_cli (the rna_lexis_stat_cli batch CLI).

Unlike the interactive menu, this entry point is argparse-based and purely
file-in/file-out with no interactive prompts, making it the most
straightforward part of the package to exercise end-to-end.
"""

import csv
import os
import tempfile
import unittest

from rna_lexis.test_cli import main


class BatchCliTestCase(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)

    def _path(self, name):
        return os.path.join(self.tmpdir.name, name)

    def _write(self, name, content):
        path = self._path(name)
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
        return path

    def _read_csv(self, path):
        with open(path, newline="", encoding="utf-8") as f:
            return list(csv.DictReader(f))


class ScoreExactTests(BatchCliTestCase):
    def test_scores_a_planted_repeated_motif(self):
        fasta = self._write("seq.fasta", ">t\ngattacagattacagattaca\n")
        motifs = self._write("motifs.txt", "gattaca\n")
        out = self._path("out.csv")

        main([
            "score-exact", "--input", fasta, "--motifs", motifs,
            "--output", out, "--enrichment-threshold", "1.0",
        ])

        rows = self._read_csv(out)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["motif"], "GATTACA")
        self.assertEqual(rows[0]["exact_count"], "3")
        self.assertIn("q_markov", rows[0])


class RankCoresTests(BatchCliTestCase):
    def test_ranks_substring_shared_by_explicit_xmotifs(self):
        # Auto-discovery (find_boundary) only surfaces xmotifs that
        # themselves repeat; to get a *family* of related xmotifs sharing
        # an inner substring (what rank-cores is actually designed to
        # score), pass --xmotifs explicitly, mirroring the fixture already
        # verified for algorithms.cores() in test_algorithms.py.
        fasta = self._write(
            "seq.fasta", ">t\naaatgtatataccctgtatataaaatgtatatagggg\n"
        )
        xmotifs = self._write(
            "xmotifs.txt",
            "aaatgtatataccc\ngggtgtatataaaa\naaatgtatatagggg\n",
        )
        out = self._path("out.csv")

        main([
            "rank-cores", "--input", fasta, "--xmotifs", xmotifs,
            "--output", out, "--candidate-min-len", "6",
            "--enrichment-threshold", "1.0",
        ])

        rows = self._read_csv(out)
        motifs = {r["motif"] for r in rows}
        self.assertIn("TGTATATA", motifs)

    def test_writes_legacy_cores_output_when_requested(self):
        fasta = self._write(
            "seq.fasta", ">t\naaatgtatataccctgtatataaaatgtatatagggg\n"
        )
        xmotifs = self._write(
            "xmotifs.txt",
            "aaatgtatataccc\ngggtgtatataaaa\naaatgtatatagggg\n",
        )
        out = self._path("out.csv")
        legacy_out = self._path("legacy.csv")

        main([
            "rank-cores", "--input", fasta, "--xmotifs", xmotifs,
            "--output", out, "--legacy-cores-output", legacy_out,
            "--candidate-min-len", "6", "--enrichment-threshold", "1.0",
        ])

        self.assertTrue(os.path.isfile(legacy_out))
        legacy_rows = self._read_csv(legacy_out)
        self.assertTrue(any("TGTATATA" in r["core"] for r in legacy_rows))


class MutationFamiliesTests(BatchCliTestCase):
    def test_reports_a_decision_per_radius(self):
        fasta = self._write("seq.fasta", ">t\ngattacagattacagattaca\n")
        motifs = self._write("motifs.txt", "gattaca\n")
        out = self._path("out.csv")

        main([
            "mutation-families", "--input", fasta, "--motifs", motifs,
            "--output", out, "--enrichment-threshold", "1.0",
        ])

        rows = self._read_csv(out)
        self.assertTrue(rows)
        self.assertIn("decision", rows[0])
        self.assertTrue(any(r["radius"] == "0" for r in rows))

    def test_writes_best_radius_output_when_requested(self):
        fasta = self._write("seq.fasta", ">t\ngattacagattacagattaca\n")
        motifs = self._write("motifs.txt", "gattaca\n")
        out = self._path("out.csv")
        best_out = self._path("best.csv")

        main([
            "mutation-families", "--input", fasta, "--motifs", motifs,
            "--output", out, "--best-output", best_out,
            "--enrichment-threshold", "1.0",
        ])

        best_rows = self._read_csv(best_out)
        self.assertEqual(len(best_rows), 1)  # one row per motif
        self.assertEqual(best_rows[0]["motif"], "GATTACA")


class GappedMotifTests(BatchCliTestCase):
    def test_finds_anchor_pairs_within_gap_range(self):
        fasta = self._write("seq.fasta", ">t\ngattacagattacagattaca\n")
        out = self._path("out.csv")

        main([
            "gapped-motif", "--input", fasta, "--left", "gat",
            "--right", "aca", "--min-gap", "0", "--max-gap", "3",
            "--output", out,
        ])

        rows = self._read_csv(out)
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(r["gap_length"] == "1" for r in rows))
        self.assertEqual(rows[0]["matched_sequence"], "GATTACA")


if __name__ == "__main__":
    unittest.main()

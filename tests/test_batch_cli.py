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


# ---------------------------------------------------------------------------
# Reading the text over an alphabet (--alphabet, --frame, --on-unknown)
# ---------------------------------------------------------------------------
import contextlib
import io
import random

from rna_lexis.alphabet import Alphabet, codon

CODONS = codon().symbols
CORE = list(CODONS[10:17])


def _codon_fasta(n=400, seed=1, planted=(20, 100, 180, 260)):
    rnd = random.Random(seed)
    units = [rnd.choice(CODONS) for _ in range(n)]
    for s in planted:
        units[s:s + 7] = CORE
    return "".join(units)


def _private_use(text):
    return [c for c in text if 0xE000 <= ord(c) <= 0xF8FF]


class AlphabetCliTests(BatchCliTestCase):
    def _run(self, argv):
        buf, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(err):
            main(argv)
        return buf.getvalue(), err.getvalue()

    def test_score_exact_over_codons_shows_codon_labels(self):
        fasta = self._write("cds.fasta", ">cds\n" + _codon_fasta() + "\n")
        motifs = self._write("motifs.txt", "".join(CORE) + "\n")
        out = self._path("out.csv")
        stdout, _ = self._run(["score-exact", "--input", fasta, "--motifs", motifs, "--output", out,
                               "--alphabet", "codon", "--enrichment-threshold", "1.0"])
        self.assertIn("Read 400 codons over 'codon'", stdout)
        rows = self._read_csv(out)
        self.assertEqual(rows[0]["motif"], " ".join(c.upper() for c in CORE))
        self.assertEqual(rows[0]["exact_count"], "4")
        self.assertEqual(rows[0]["markov_order"], "0")                  # the codon default
        self.assertEqual(rows[0]["statistically_supported"], "True")
        self.assertEqual(_private_use(open(out, encoding="utf-8").read()), [])

    def test_q_value_is_corrected_over_all_possible_words(self):
        fasta = self._write("cds.fasta", ">cds\n" + _codon_fasta() + "\n")
        motifs = self._write("motifs.txt", "".join(CORE) + "\n")
        out = self._path("out.csv")
        self._run(["score-exact", "--input", fasta, "--motifs", motifs, "--output", out, "--alphabet", "codon"])
        row = self._read_csv(out)[0]
        self.assertAlmostEqual(float(row["q_markov"]), min(1.0, float(row["p_markov"]) * 64 ** 7), delta=1e-9 + 1e-6 * float(row["q_markov"]))

    def test_default_without_an_alphabet_is_unchanged(self):
        fasta = self._write("seq.fasta", ">t\ngattacagattacagattaca\n")
        motifs = self._write("motifs.txt", "gattaca\n")
        a, b = self._path("a.csv"), self._path("b.csv")
        self._run(["score-exact", "--input", fasta, "--motifs", motifs, "--output", a, "--enrichment-threshold", "1.0"])
        self._run(["score-exact", "--input", fasta, "--motifs", motifs, "--output", b, "--enrichment-threshold", "1.0",
                   "--markov-order", "1"])
        self.assertEqual(open(a, encoding="utf-8").read(), open(b, encoding="utf-8").read())

    def test_frame_and_unknown_symbols(self):
        text = "cc" + "".join(CORE) * 2 + "nnn" + "".join(CODONS[:3])
        fasta = self._write("cds.fasta", ">cds\n" + text + "\n")
        motifs = self._write("motifs.txt", "".join(CORE) + "\n")
        out = self._path("out.csv")
        stdout, err = self._run(["score-exact", "--input", fasta, "--motifs", motifs, "--output", out,
                                 "--alphabet", "codon", "--frame", "2", "--enrichment-threshold", "1.0"])
        self.assertIn("1 skipped", stdout)
        self.assertIn("unknown symbol(s) skipped", err)
        self.assertEqual(self._read_csv(out)[0]["exact_count"], "2")
        with self.assertRaises(SystemExit) as cm:
            self._run(["score-exact", "--input", fasta, "--motifs", motifs, "--output", out,
                       "--alphabet", "codon", "--frame", "2", "--on-unknown", "error"])
        self.assertIn("not a symbol", str(cm.exception))

    def test_fasta_records_are_segments_that_motifs_cannot_span(self):
        half = "".join(CORE)
        fasta = self._write("two.fasta", f">a\n{half[:12]}\n>b\n{half[12:]}\n")
        motifs = self._write("motifs.txt", half + "\n")
        out = self._path("out.csv")
        self._run(["score-exact", "--input", fasta, "--motifs", motifs, "--output", out,
                   "--alphabet", "codon", "--enrichment-threshold", "1.0"])
        self.assertEqual(self._read_csv(out)[0]["exact_count"], "0")        # split across two records

    def test_rank_cores_over_codons(self):
        # a core family: the core with two different extensions, as the pipeline needs two distinct xmotifs
        rnd = random.Random(3)
        units = [rnd.choice(CODONS) for _ in range(400)]
        for k, s in enumerate((20, 100, 180, 260, 340)):
            units[s:s + 9] = CORE + (list(CODONS[30:32]) if k % 2 == 0 else list(CODONS[40:42]))
        fasta = self._write("cds.fasta", ">cds\n" + "".join(units) + "\n")
        out = self._path("ranked.csv")
        self._run(["rank-cores", "--input", fasta, "--output", out, "--alphabet", "codon",
                   "--min-xmotif-len", "4", "--candidate-min-len", "3"])
        rows = self._read_csv(out)
        supported = [r["motif"] for r in rows if r["statistically_supported"] == "True"]
        core = " ".join(c.upper() for c in CORE)
        self.assertTrue(any(m in core for m in supported), supported)
        self.assertEqual(_private_use(open(out, encoding="utf-8").read()), [])

    def test_mutation_families_count_variants_over_64_symbols(self):
        fasta = self._write("cds.fasta", ">cds\n" + _codon_fasta() + "\n")
        motifs = self._write("motifs.txt", "".join(CORE) + "\n")
        out = self._path("fam.csv")
        self._run(["mutation-families", "--input", fasta, "--motifs", motifs, "--output", out,
                   "--alphabet", "codon", "--mutr", "0.34"])
        rows = {r["radius"]: r for r in self._read_csv(out)}
        self.assertEqual(rows["1"]["hamming_variant_count"], str(1 + 7 * 63))

    def test_gapped_motif_over_codons(self):
        fasta = self._write("cds.fasta", ">cds\n" + _codon_fasta() + "\n")
        out = self._path("gap.csv")
        left, right = "".join(CORE[:2]), "".join(CORE[4:6])
        self._run(["gapped-motif", "--input", fasta, "--left", left, "--right", right,
                   "--min-gap", "0", "--max-gap", "3", "--output", out, "--alphabet", "codon"])
        rows = self._read_csv(out)
        self.assertEqual(len(rows), 4)
        self.assertEqual(rows[0]["gap_length"], "2")
        self.assertEqual(rows[0]["gap_sequence"], " ".join(c.upper() for c in CORE[2:4]))
        self.assertTrue(rows[0]["pattern"].startswith(" ".join(c.upper() for c in CORE[:2]) + "[gap:0-3]"))

    def test_an_alphabet_file_with_non_latin_letters(self):
        greek = Alphabet(name="Greek", symbols=tuple("αβγδεζ"))
        path = greek.save(self._path("greek.json"))
        text = ("αβγδ" * 5 + "ζε" * 3 + "\n")
        fasta = self._write("g.txt", text)
        motifs = self._write("m.txt", "αβγδ\n")
        out = self._path("out.csv")
        self._run(["score-exact", "--input", fasta, "--motifs", motifs, "--output", out,
                   "--alphabet", path, "--enrichment-threshold", "1.0"])
        row = self._read_csv(out)[0]
        self.assertEqual(row["motif"], "ΑΒΓΔ")
        self.assertEqual(row["exact_count"], "5")

    def test_saved_session_motifs_are_used_as_they_are(self):
        import json
        a = codon()
        core = a.encode(CORE)
        session = {"corelist": [core], "xmotifs": [core], "alphabet": a.to_dict()}
        path = self._write("sess.json", json.dumps(session))
        fasta = self._write("cds.fasta", ">cds\n" + _codon_fasta() + "\n")
        out = self._path("out.csv")
        self._run(["score-exact", "--input", fasta, "--motifs", path, "--output", out, "--alphabet", "codon",
                   "--enrichment-threshold", "1.0"])
        self.assertEqual(self._read_csv(out)[0]["exact_count"], "4")
        other = self._write("sess2.json", json.dumps(dict(session, alphabet=Alphabet(name="x", symbols=("a",)).to_dict())))
        with self.assertRaises(SystemExit):
            self._run(["score-exact", "--input", fasta, "--motifs", other, "--output", out, "--alphabet", "codon"])

    def test_unknown_alphabet_is_an_error(self):
        fasta = self._write("cds.fasta", ">cds\n" + _codon_fasta() + "\n")
        motifs = self._write("motifs.txt", "".join(CORE) + "\n")
        with self.assertRaises(SystemExit):
            self._run(["score-exact", "--input", fasta, "--motifs", motifs, "--output", self._path("o.csv"),
                       "--alphabet", "no-such-alphabet"])

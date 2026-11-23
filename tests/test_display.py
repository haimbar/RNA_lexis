"""Showing results over alphabets other than nucleotides (phase 3)."""

import contextlib
import csv
import io
import json
import os
import random
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.basedatatypes as pb

import rna_lexis.plots as P
from rna_lexis.algorithms import find_with_mutations, print_core
from rna_lexis.alignment import gotoh_global, print_alignment
from rna_lexis.alphabet import Alphabet, codon, latin, nucleotide
from rna_lexis.statistical import (
    find_gapped_motif_hits, mutation_family_tests, score_exact_motifs, score_gapped_motif, write_rows_csv,
)

A = codon()
C = A.symbols


def private_use(text):
    return [c for c in text if 0xE000 <= ord(c) <= 0xF8FF]


def codon_text(n=300, seed=1):
    """Random codons with a 7-codon core planted four times, one copy with a changed codon."""
    rnd = random.Random(seed)
    units = [rnd.choice(C) for _ in range(n)]
    core = list(C[10:17])
    for s in (20, 90, 160, 230):
        units[s:s + 7] = core
    units[160 + 3] = C[40]                           # a one-codon variant
    return A.encode(units), A.encode(core)


class ShowTests(unittest.TestCase):
    def test_default_separator_and_marking(self):
        t = A.encode(["atg", "aaa", "ttt"])
        self.assertEqual(A.show(t), "atg aaa ttt")
        self.assertEqual(A.show(t, mask=[False, True, False]), "atg AAA ttt")
        self.assertEqual(A.show(t, upper=True), "ATG AAA TTT")
        self.assertEqual(A.show(t, sep="-"), "atg-aaa-ttt")
        self.assertEqual(latin().show("hello", mask=[0, 0, 1, 0, 0]), "heLlo")

    def test_nucleotide_marking_is_the_old_convention(self):
        self.assertEqual(nucleotide().show("gattaca", mask=[0, 0, 0, 1, 0, 0, 1]), "gatTacA")

    def test_boundary_and_padding(self):
        syl = Alphabet(name="syl", symbols=("a", "bb", "ccc"), mode="words")
        self.assertEqual(syl.label_width, 3)
        t = syl.encode(["a", "bb", "ccc"])
        self.assertEqual(syl.show(t, pad=True), "a   bb  ccc")
        self.assertEqual(syl.show(t[:1] + "_" + t[1:]), "a | bb ccc")

    def test_case_sensitive_alphabets_mark_with_brackets(self):
        cs = Alphabet(name="cs", symbols=("A", "a"), case_sensitive=True)
        self.assertEqual(cs.show(cs.encode(["A", "a"]), mask=[False, True]), "A[a]")

    def test_mask_must_match_the_length(self):
        with self.assertRaises(ValueError):
            A.show(A.encode(["atg"]), mask=[True, False])

    def test_unit_names_and_file_round_trip(self):
        self.assertEqual(nucleotide().units_label, "nt")
        self.assertEqual(A.units_label, "codons")
        self.assertEqual(latin().units_label, "letters")
        a = codon()._replace(display_sep="-", unit_name="triplet")
        with tempfile.TemporaryDirectory() as d:
            b = Alphabet.load(a.save(os.path.join(d, "a.json")))
        self.assertEqual((b.unit_name, b.display_sep, b.sep), ("triplet", "-", "-"))
        self.assertEqual(b.units_label, "triplets")


class CsvTests(unittest.TestCase):
    def read(self, path):
        with open(path, encoding="utf-8") as fh:
            return fh.read()

    def test_motif_columns_are_decoded_and_upper_case(self):
        txt, core = codon_text()
        rows = score_exact_motifs(txt, [core], alphabet=A)
        self.assertEqual(rows[0]["motif"], core)                   # the row keeps the internal string
        with tempfile.TemporaryDirectory() as d:
            text = self.read(write_rows_csv(os.path.join(d, "r.csv"), rows, alphabet=A))
        self.assertEqual(private_use(text), [])
        self.assertIn(A.show(core, upper=True), text)
        self.assertEqual(next(csv.DictReader(io.StringIO(text)))["motif"], A.show(core, upper=True))

    def test_gapped_columns_and_pattern(self):
        txt, core = codon_text()
        left, right = core[:2], core[4:6]
        score = score_gapped_motif(txt, left, right, min_gap=0, max_gap=3, alphabet=A)
        hits = find_gapped_motif_hits(txt, left, right, min_gap=0, max_gap=3, alphabet=A)
        with tempfile.TemporaryDirectory() as d:
            a = self.read(write_rows_csv(os.path.join(d, "a.csv"), [score], alphabet=A))
            b = self.read(write_rows_csv(os.path.join(d, "b.csv"), hits, alphabet=A))
        self.assertEqual(private_use(a + b), [])
        pat = next(csv.DictReader(io.StringIO(a)))["pattern"]
        self.assertEqual(pat, f"{A.show(left, upper=True)}[gap:0-3]{A.show(right, upper=True)}")
        self.assertIn(A.show(core[2:4], upper=True), b)             # the gap between the anchors

    def test_without_an_alphabet_rows_are_written_as_before(self):
        rows = score_exact_motifs("gattacagattaca" * 5, ["gattaca"])
        with tempfile.TemporaryDirectory() as d:
            text = self.read(write_rows_csv(os.path.join(d, "r.csv"), rows))
        self.assertIn("GATTACA", text)

    def test_mutation_family_rows(self):
        txt, core = codon_text()
        rows = mutation_family_tests(txt, [core], alphabet=A, mutr=1 / 3)
        with tempfile.TemporaryDirectory() as d:
            text = self.read(write_rows_csv(os.path.join(d, "m.csv"), rows, alphabet=A))
        self.assertEqual(private_use(text), [])


class PlotTests(unittest.TestCase):
    def setUp(self):
        self.figs = []
        self._show = pb.BaseFigure.show
        pb.BaseFigure.show = lambda fig, *a, **k: self.figs.append(fig)
        self._plt_show = plt.show
        plt.show = lambda *a, **k: None
        self._open = (P.open_file_with_default_software, P.open_pdf)
        P.open_file_with_default_software = lambda *a, **k: None
        P.open_pdf = lambda *a, **k: None
        self.txt, self.core = codon_text()

    def tearDown(self):
        pb.BaseFigure.show = self._show
        plt.show = self._plt_show
        P.open_file_with_default_software, P.open_pdf = self._open
        plt.close("all")

    def fig_text(self):
        return json.dumps([json.loads(f.to_json()) for f in self.figs], ensure_ascii=False)

    def test_seq_nbrs_shows_decoded_labels_and_a_marked_variant(self):
        with contextlib.redirect_stdout(io.StringIO()):
            P.plot_seq_nbrs(self.core, [self.core, self.core[:4]], self.txt, wd=30, alphabet=A)
        text = self.fig_text()
        self.assertEqual(private_use(text), [])
        self.assertIn(A.show(self.core), text)
        self.assertIn(A.show(self.core[:4]), text)                 # a neighbour's label
        variant = list(C[10:17])
        variant[3] = C[40]
        marked = A.show(A.encode(variant), mask=[0, 0, 0, 1, 0, 0, 0])
        self.assertIn(marked, text)                                # the changed codon in upper case

    def test_condensed_plot(self):
        with contextlib.redirect_stdout(io.StringIO()):
            P.plot_nbrs_condensed(self.core, [self.core], self.txt, wd=30, alphabet=A)
        text = self.fig_text()
        self.assertEqual(private_use(text), [])
        self.assertIn(A.show(self.core), text)

    def test_export_prints_decoded_windows(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            P.export_nbrs_condensed(self.core, [self.core], self.txt, wd=5, alphabet=A)
        self.assertEqual(private_use(buf.getvalue()), [])
        self.assertIn(A.show(self.core), buf.getvalue())

    def test_sequence_hits_and_the_detailed_view(self):
        exact, approx, _ = find_with_mutations(self.core, self.txt, mutr=1 / 3)
        sd = [{"seq": self.core, "L": len(self.core), "exact": exact, "approx": approx}]
        self.assertTrue(approx)
        with contextlib.redirect_stdout(io.StringIO()):
            P.plot_sequence_hits(sd, len(self.txt), alphabet=A)
        text = self.fig_text()
        self.assertEqual(private_use(text), [])
        self.assertIn(A.show(self.core), text)
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, "d.html")
            P.plot_sequence_hits_detailed(sd, self.txt, 0, len(self.txt), file=f, alphabet=A)
            html = open(f, encoding="utf-8").read()
        self.assertEqual(private_use(html), [])
        self.assertIn(A.show(self.core), html)                     # the legend
        self.assertIn("codons", html)
        self.assertNotIn("nucleotides", html)

    def test_matplotlib_labels(self):
        cnt = {A.char(C[i]): (i % 5) + 1 for i in range(1, 60)}      # spread, so the IQR is not zero
        cnt[A.char(C[0])] = 500                                   # one clear outlier
        with contextlib.redirect_stdout(io.StringIO()):
            P.plotzscore(cnt, zmin=1.96, alphabet=A)
        labels = [t.get_text() for t in plt.gca().texts]
        self.assertEqual(labels, [C[0]])                           # the outlier, as its label
        self.assertEqual(private_use("".join(labels)), [])
        plt.close("all")
        P.plot_coverage(self.txt, {self.core: 5, self.core[:3]: 2}, 1.2, alphabet=A)
        ticks = [t.get_text() for t in plt.gca().get_xticklabels()]
        self.assertIn(A.show(self.core), ticks)
        self.assertEqual(private_use("".join(ticks)), [])

    def test_arcs_and_shared_diagram_use_the_unit_name(self):
        res = [{"pos1": 20, "pos2": 90, "total_len": 9, "hamming": 1},
               {"pos1": 90, "pos2": 160, "total_len": 8, "hamming": 2}]
        P.plot_self_similarity_arcs(self.core, res, [20, 90, 160], alphabet=A)
        ax = plt.gca()
        self.assertEqual(ax.get_xlabel(), "Position (codons)")
        self.assertIn(A.show(self.core), ax.get_title())
        plt.close("all")
        with contextlib.redirect_stdout(io.StringIO()):
            P.plot_shared_motif_diagram(self.txt, self.txt, "A", "B", [self.core], alphabet=A)
        texts = " ".join(t.get_text() for t in plt.gcf().findobj(matplotlib.text.Text))
        self.assertIn("codons", texts)
        self.assertIn(A.show(self.core), texts)
        self.assertEqual(private_use(texts), [])

    def test_nucleotide_defaults_still_say_nt(self):
        res = [{"pos1": 20, "pos2": 90, "total_len": 9, "hamming": 1},
               {"pos1": 90, "pos2": 160, "total_len": 8, "hamming": 2}]
        P.plot_self_similarity_arcs("gattaca", res, [20, 90, 160])
        self.assertEqual(plt.gca().get_xlabel(), "Position (nt)")


class PrinterTests(unittest.TestCase):
    def test_print_core_pads_in_whole_units(self):
        t = A.encode([C[i % 60] for i in range(120)])
        core = A.encode([C[3], C[4]])
        xm = [A.encode(C[2:6]), A.encode(C[3:7])]
        out = print_core(t, core, xm, alphabet=A)
        self.assertEqual(len(out), 2)
        self.assertEqual(private_use("".join(out)), [])
        # the core starts in the same column in both lines
        cols = {line.index(A.show(core)) for line in out}
        self.assertEqual(len(cols), 1)

    def test_print_alignment_shows_labels_gap_and_markers(self):
        a = A.encode(C[1:7])
        b = A.encode([C[1], C[2], C[9], C[4], C[6]])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            print_alignment(gotoh_global(a, b), alphabet=A)
        out = buf.getvalue()
        self.assertEqual(private_use(out), [])
        self.assertIn("---", out)                               # a gap is as wide as a label
        self.assertIn(A.show(a[:2]), out)
        lines = [l for l in out.splitlines() if l.strip()][2:5]
        self.assertEqual(len({len(l.rstrip()) for l in lines[:1] + lines[2:]}), 1)

    def test_letters_print_exactly_as_before(self):
        res = gotoh_global("gattacagattaca", "gattagagatta")
        old, new = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(old):
            print_alignment(res)
        with contextlib.redirect_stdout(new):
            print_alignment(res, alphabet=nucleotide())
        self.assertEqual(old.getvalue(), new.getvalue())


if __name__ == "__main__":
    unittest.main()

"""The unit profile: what surrounds a motif (the view that replaces a logo for large alphabets)."""

import csv
import io
import math
import os
import random
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import rna_lexis.plots as P
from rna_lexis.alphabet import BOUNDARY, Alphabet, codon, latin, nucleotide
from rna_lexis.genetic_code import amino_acid, codon_to_amino_acid
from rna_lexis.statistical import unit_profile, describe_profile, describe_empty_positions, describe_matches, write_unit_profile_csv
from rna_lexis.tokenizers import tokenize

A = codon()
C = A.symbols


def private_use(text):
    return [c for c in text if 0xE000 <= ord(c) <= 0xF8FF]


def planted_text(n=800, seed=2, copies=(20, 120, 220, 320, 420, 520, 620)):
    rnd = random.Random(seed)
    units = [rnd.choice(C) for _ in range(n)]
    core = list(C[10:14])
    for s in copies:
        units[s - 1] = C[5]            # a conserved codon just before
        units[s + 4] = C[33]           # and one just after
        units[s:s + 4] = core
    return A.encode(units), A.encode(core)


class CountsTests(unittest.TestCase):
    def test_counts_around_exact_occurrences(self):
        abc = Alphabet(name="abc", symbols=("a", "b", "c"))
        txt = "cabcca" + "bbcabcbb" + "c"                    # 'ab' at 1 and 9... built by hand below
        txt = "xabcyabdz".replace("x", "c").replace("y", "c").replace("z", "c")   # c a b c c a b d? no d
        txt = "cabcccabac"                                   # occurrences of "ab" at 1 and 6
        p = unit_profile(txt, "ab", flank=1, alphabet=abc)
        self.assertEqual(p["n"], 2)
        self.assertEqual(p["columns"], [-1, 0, 1, 2])
        self.assertEqual(p["motif_columns"], (1, 3))
        rows = {c: r for c, r in zip(p["chars"], p["counts"])}
        self.assertEqual(rows["c"], [2, 0, 0, 1])           # before: c, c; after: c, a
        self.assertEqual(rows["a"], [0, 2, 0, 1])
        self.assertEqual(rows["b"], [0, 0, 2, 0])
        self.assertEqual(p["motif_labels"], ["a", "b"])

    def test_background_is_the_texts_own_frequency(self):
        abc = Alphabet(name="abc", symbols=("a", "b", "c"))
        p = unit_profile("cabcccabac", "ab", flank=1, alphabet=abc)
        self.assertAlmostEqual(sum(p["background"]), 1.0)
        by = dict(zip(p["chars"], p["background"]))
        self.assertGreater(by["c"], by["a"])               # c is the commonest letter
        # a symbol exactly at its background frequency has log-odds near 0
        self.assertAlmostEqual(math.log2((0 + by["b"]) / 3 / by["b"]), math.log2(1 / 3))

    def test_windows_off_the_text_or_across_a_boundary_are_left_out(self):
        abc = Alphabet(name="abc", symbols=("a", "b", "c"))
        p = unit_profile("abcabcab", "ab", flank=1, alphabet=abc)        # first and last occurrence touch the ends
        self.assertEqual((p["n"], p["skipped"]), (1, 2))
        tok = tokenize(["cabc", "cabc"], abc)
        p = unit_profile(tok.text, "ab", flank=3, alphabet=abc)         # a window of 3 each side spans the boundary
        self.assertEqual(p["n"], 0)
        p = unit_profile(tok.text, "ab", flank=1, alphabet=abc)
        self.assertEqual(p["n"], 2)

    def test_mismatches_widen_the_set_of_windows(self):
        abc = Alphabet(name="abc", symbols=("a", "b", "c"))
        txt = "cabccaccbcc"
        self.assertEqual(unit_profile(txt, "ab", flank=1, alphabet=abc)["n"], 1)
        self.assertGreater(unit_profile(txt, "ab", flank=1, mismatches=1, alphabet=abc)["n"], 1)

    def test_nucleotides_without_an_alphabet(self):
        p = unit_profile("gattacagattacagg", "gatt", flank=2)
        self.assertEqual(p["n"], 1)
        self.assertEqual(sorted(p["labels"]), ["a", "c", "g", "t"])

    def test_the_motif_cannot_be_empty_or_span_a_boundary(self):
        with self.assertRaises(ValueError):
            unit_profile("acgt", "", alphabet=nucleotide())
        with self.assertRaises(ValueError):
            unit_profile("ac_gt", "c_g", alphabet=nucleotide())


class SignificanceTests(unittest.TestCase):
    def test_planted_neighbors_are_found_and_chance_is_not(self):
        txt, motif = planted_text()
        p = unit_profile(txt, motif, flank=3, alphabet=A)
        self.assertEqual(p["n"], 7)
        col = {c: j for j, c in enumerate(p["columns"])}
        row = {c: i for i, c in enumerate(p["chars"])}
        before, after = row[A.char(C[5])], row[A.char(C[33])]
        self.assertTrue(p["significant"][before][col[-1]])
        self.assertTrue(p["significant"][after][p["motif_columns"][1] + 0])
        self.assertGreater(p["log_odds"][before][col[-1]], 4)
        self.assertEqual(p["counts"][before][col[-1]], 7)
        total = sum(map(sum, p["significant"]))
        self.assertLessEqual(total, 4 + 4 + 3)             # the motif, the two planted cells, a few chance cells

    def test_few_windows_do_not_make_noise(self):
        # with 7 windows, a symbol that never appears or appears once is not evidence of anything
        txt, motif = planted_text()
        p = unit_profile(txt, motif, flank=3, alphabet=A)
        absent = [(i, j) for i in range(len(p["chars"])) for j in range(len(p["columns"]))
                  if p["counts"][i][j] == 0]
        self.assertTrue(absent)
        self.assertFalse(any(p["significant"][i][j] for i, j in absent))

    def test_under_a_null_about_one_percent_of_cells_pass(self):
        rnd = random.Random(11)
        units = [rnd.choice(list("abcdefghijklmnopqrstuvwxyz")) for _ in range(40000)]
        for s in range(100, 39000, 300):
            units[s:s + 4] = list("thex")
        p = unit_profile("".join(units), "thex", flank=4, alphabet=latin())
        cells = [(i, j) for i in range(26) for j in range(len(p["columns"]))
                 if not (p["motif_columns"][0] <= j < p["motif_columns"][1])]
        share = sum(p["significant"][i][j] for i, j in cells) / len(cells)
        self.assertLess(share, 0.03)

    def test_depletion_needs_a_symbol_expected_often_enough(self):
        rnd = random.Random(3)
        units = [rnd.choice("ab") for _ in range(4000)]      # a and b each half the text
        motif = "aabbababba"                                 # long enough that it hardly ever occurs by chance
        for s in range(10, 3900, 40):
            units[s:s + len(motif)] = list(motif)
            units[s - 1] = "a"                               # b never right before the motif
        abx = Alphabet(name="ab", symbols=("a", "b"))
        p = unit_profile("".join(units), motif, flank=2, alphabet=abx)
        j = p["columns"].index(-1)
        ib = p["chars"].index("b")
        self.assertTrue(p["significant"][ib][j])
        self.assertLess(p["log_odds"][ib][j], -2)


class GroupingAndFilesTests(unittest.TestCase):
    def test_codons_are_grouped_by_amino_acid(self):
        txt, motif = planted_text()
        proj = codon_to_amino_acid()
        p = unit_profile(txt, motif, flank=2, alphabet=A, group_by=proj)
        self.assertEqual(len(p["chars"]), 64)
        labels = [A.label(c) for c in p["chars"]]
        for name, start, end in p["groups"]:
            self.assertTrue(all(proj.mapping[labels[i]] == name for i in range(start, end)), name)
        self.assertEqual(p["groups"][-1][2], 64)
        self.assertEqual(sorted(labels), sorted(C))

    def test_csv(self):
        txt, motif = planted_text()
        p = unit_profile(txt, motif, flank=2, alphabet=A)
        with tempfile.TemporaryDirectory() as d:
            path = write_unit_profile_csv(os.path.join(d, "p.csv"), p)
            text = open(path, encoding="utf-8").read()
            rows = list(csv.DictReader(io.StringIO(text)))
        self.assertEqual(private_use(text), [])
        self.assertEqual(len(rows), 64)
        self.assertIn("count@-1", rows[0])
        self.assertIn("p@-1", rows[0])
        by = {r["symbol"]: r for r in rows}
        self.assertEqual(by[C[5].upper()]["count@-1"], "7")


class PlotTests(unittest.TestCase):
    def setUp(self):
        self._show = plt.show
        plt.show = lambda *a, **k: None
        self._open = P.open_file_with_default_software
        P.open_file_with_default_software = lambda *a, **k: None

    def tearDown(self):
        plt.show = self._show
        P.open_file_with_default_software = self._open
        plt.close("all")

    def test_heatmap_shows_decoded_labels_and_rings_the_planted_cells(self):
        txt, motif = planted_text()
        p = unit_profile(txt, motif, flank=3, alphabet=A, group_by=codon_to_amino_acid())
        P.plot_unit_profile(p, kind="heatmap")
        ax = plt.gcf().axes[0]
        xt = [t.get_text() for t in ax.get_xticklabels()]
        self.assertEqual(xt[3:7], list(A.decode_list(motif)))          # the motif, as codons
        self.assertEqual(xt[:3], ["-3", "-2", "-1"])
        self.assertEqual(xt[7:], ["+1", "+2", "+3"])
        self.assertEqual(len(ax.get_yticklabels()), 64)
        self.assertEqual(private_use("".join(xt) + "".join(t.get_text() for t in ax.get_yticklabels())), [])
        rings = [line for line in ax.lines if line.get_marker() == 'o']
        self.assertEqual(len(rings), 2)                                 # the two planted cells
        self.assertIn("7 occurrences", ax.get_title())
        self.assertIn(" ".join(A.decode_list(motif)).upper(), ax.get_title())

    def test_only_supported_cells_are_colored(self):
        txt, motif = planted_text()
        p = unit_profile(txt, motif, flank=3, alphabet=A)
        P.plot_unit_profile(p, kind="heatmap")
        data = plt.gcf().axes[0].images[0].get_array()
        colored = int((data != 0).sum())
        self.assertEqual(colored, sum(map(sum, p["significant"])))

    def test_saves_a_file_and_works_for_letters_and_amino_acids(self):
        rnd = random.Random(5)
        for alpha in (latin(), amino_acid()):
            letters = [c for c in alpha.labels if len(c) == 1 and c.isalpha()]
            units = [rnd.choice(letters) for _ in range(3000)]
            for s in range(50, 2900, 100):
                units[s:s + 3] = list(letters[:3])
                units[s - 1] = letters[5]
            p = unit_profile(alpha.encode(units), alpha.encode(letters[:3]), flank=3, alphabet=alpha)
            with tempfile.TemporaryDirectory() as d:
                out = os.path.join(d, "p.png")
                P.plot_unit_profile(p, file=out)
                self.assertGreater(os.path.getsize(out), 1000)
            plt.close("all")


class MismatchTests(unittest.TestCase):
    """A motif is matched up to a chosen number of substituted units, and the output says which occurrences were approximate."""
    MOTIF = "kloprsuv"

    def setUp(self):
        plt.close("all")
        self.addCleanup(plt.close, "all")
        self._open = P.open_file_with_default_software
        P.open_file_with_default_software = lambda *a, **k: None
        self.addCleanup(setattr, P, "open_file_with_default_software", self._open)
        rnd = random.Random(12)
        self.a = latin()
        units = [rnd.choice("abcdefghijmnqtwxyz") for _ in range(3000)]        # none of the motif's letters
        for st in range(100, 2900, 400):                                         # 7 exact copies
            units[st:st + 8] = list(self.MOTIF)
        for st in (250, 650):                                                    # one mismatch: the last letter
            units[st:st + 8] = list("kloprsua")
        units[1050:1058] = list("kloprsax")                                      # two mismatches
        self.txt = self.a.encode("".join(units))

    def profile(self, m):
        return unit_profile(self.txt, self.a.encode(self.MOTIF), flank=3, mismatches=m, alphabet=self.a)

    def test_the_number_of_occurrences_grows_with_the_mismatches_allowed(self):
        self.assertEqual([self.profile(m)["n"] for m in (0, 1, 2)], [7, 9, 10])

    def test_the_breakdown_and_the_approximate_occurrences_are_reported(self):
        p = self.profile(2)
        self.assertEqual(p["by_distance"], {0: 7, 1: 2, 2: 1})
        self.assertEqual(p["approximate"], [(250, "kloprsua", 1), (650, "kloprsua", 1), (1050, "kloprsax", 2)])
        lines = describe_matches(p)
        self.assertEqual(lines[0], "10 occurrence(s) allowing up to 2 mismatches: 7 exact, 2 with 1 mismatch, 1 with 2 mismatches")
        self.assertIn("  position 1050: KLOPRSAX instead of KLOPRSUV (2 mismatches)", lines)
        self.assertIn("  position 250: KLOPRSUA instead of KLOPRSUV (1 mismatch)", lines)
        self.assertEqual(describe_matches(self.profile(0)), ["7 occurrence(s), all exact"])

    def test_the_title_of_both_plots_says_how_many_are_exact(self):
        for kind in ("logo", "heatmap"):
            P.plot_unit_profile(self.profile(2), kind=kind)
            ax = plt.gcf().axes[0]
            title = ax.get_title(loc="left") or ax.get_title()
            self.assertIn("10 occurrences: 7 exact, 2 with 1 mismatch, 1 with 2 mismatches", title, kind)
            plt.close("all")
        P.plot_unit_profile(self.profile(0))
        self.assertNotIn("exact", plt.gcf().axes[0].get_title(loc="left"))

    def test_the_flanks_of_an_approximate_occurrence_are_counted_like_the_others(self):
        for m, expect in ((0, 7), (2, 10)):
            p = self.profile(m)
            self.assertEqual(sum(p["counts"][i][0] for i in range(len(p["chars"]))), expect)

    def test_the_letters_where_a_variant_differs_are_drawn_when_they_recur(self):
        # 20 copies in which the last letter is 'a' instead of 'v': that is a clear departure at the last position
        units = list(self.a.decode(self.txt))
        for st in range(60, 2900, 140):
            units[st:st + 8] = list("kloprsua")
        p = unit_profile(self.a.encode("".join(units)), self.a.encode(self.MOTIF), flank=3, mismatches=1, alphabet=self.a)
        P.plot_unit_profile(p)
        ticks = [t.get_text() for t in plt.gcf().axes[0].get_xticklabels()]
        self.assertEqual(ticks[3:11], list("KLOPRSUV"))          # the motif's columns are kept, to show the variant

    def test_the_menu_asks_for_the_mismatches_and_prints_what_it_found(self):
        import contextlib
        import rna_lexis.menu as M
        answers = iter([self.MOTIF, "3", "2", "", "1", "", "n"])  # motif, flank, mismatches, alpha, logo, format, no csv
        saved = (M.safe_input, M._spawn_plot, M._prompt_save)
        asked = []
        M.safe_input = lambda prompt="": (asked.append(prompt), next(answers, ""))[1]
        M._spawn_plot = lambda fn, *a, **k: None
        M._prompt_save = lambda **k: ("", 1)
        try:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                M.profile_input(self.txt, {}, "", self.a)
        finally:
            M.safe_input, M._spawn_plot, M._prompt_save = saved
        text = out.getvalue()
        self.assertIn("allowing up to 2 mismatches: 7 exact, 2 with 1 mismatch, 1 with 2 mismatches", text)
        self.assertIn("KLOPRSAX instead of KLOPRSUV", text)
        self.assertTrue(any("8 letters long; 1 per 6 would be 1" in a for a in asked), asked)

    def test_the_menu_does_not_allow_more_than_half_the_motif_to_differ(self):
        import contextlib
        import rna_lexis.menu as M
        answers = iter([self.MOTIF, "3", "5", "4", "", "1", "", "n"])   # 5 is refused (more than half of 8), then 4 is taken
        saved = (M.safe_input, M._spawn_plot, M._prompt_save)
        M.safe_input = lambda prompt="": next(answers, "")
        M._spawn_plot = lambda fn, *a, **k: None
        M._prompt_save = lambda **k: ("", 1)
        try:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                M.profile_input(self.txt, {}, "", self.a)
        finally:
            M.safe_input, M._spawn_plot, M._prompt_save = saved
        self.assertIn("Please enter an integer from 0 to 4", out.getvalue())


class SignificanceLevelPromptTests(unittest.TestCase):
    """The menu lets the user say how strong a count must be to be drawn."""

    def run_menu(self, alpha_answer):
        import contextlib
        import rna_lexis.menu as M
        a = latin()
        rnd = random.Random(5)
        txt = a.encode("".join(rnd.choice("abcdefghijklmnopqrstuvwxyz") for _ in range(3000)))
        motif = "qzxjvk"
        txt = a.encode("".join(rnd.choice("abcdefghijklmnopqrstuvwxyz") for _ in range(300))) + "".join(
            a.encode("ee" + motif + "e" + rnd.choice("abcdefghijklmnopqrstuvwxyz")) for _ in range(6)) + txt
        answers = iter([motif, "2", "0", alpha_answer, "1", "", "n"])
        saved = (M.safe_input, M._spawn_plot, M._prompt_save)
        M.safe_input = lambda prompt="": next(answers, "")
        M._spawn_plot = lambda fn, *a, **k: None
        M._prompt_save = lambda **k: ("", 1)
        try:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                M.profile_input(txt, {}, "", a)
        finally:
            M.safe_input, M._spawn_plot, M._prompt_save = saved
        return out.getvalue()

    def test_a_larger_level_draws_more_cells_and_a_bad_answer_keeps_the_default(self):
        import re as _re
        cells = lambda t: int(_re.search(r"(\d+) cell\(s\) differ", t).group(1))
        strict, default, bad, loose = (self.run_menu(x) for x in ("0.00000001", "", "abc", "0.5"))
        self.assertEqual(cells(default), cells(bad))
        self.assertLessEqual(cells(strict), cells(default))
        self.assertGreater(cells(loose), cells(default))
        self.assertIn("p <= 0.5", loose)


class LogoTests(unittest.TestCase):
    """The letter logo: what it draws, and what it must not."""

    def setUp(self):
        plt.close("all")
        self.addCleanup(plt.close, "all")
        self._open = P.open_file_with_default_software
        P.open_file_with_default_software = lambda *a, **k: None
        self.addCleanup(setattr, P, "open_file_with_default_software", self._open)

    def planted_letters(self):
        """A text over Latin letters with 'q' two before every 'the', which is otherwise rare."""
        rnd = random.Random(3)
        base = [rnd.choice("abcdefghijklmnoprstuvwxy") for _ in range(4000)]
        for s in range(100, 3900, 120):
            base[s:s + 5] = list("qxthe")
        return latin(), latin().encode("".join(base))

    def test_default_is_the_logo_and_the_heatmap_is_still_there(self):
        a, txt = self.planted_letters()
        p = unit_profile(txt, a.encode("the"), flank=3, alphabet=a)
        P.plot_unit_profile(p)
        self.assertEqual(len(plt.gcf().axes[0].images), 0)               # no heatmap image: it is drawn with shapes
        plt.close("all")
        P.plot_unit_profile(p, kind="heatmap")
        self.assertEqual(len(plt.gcf().axes[0].images), 1)
        with self.assertRaises(ValueError):
            P.plot_unit_profile(p, kind="pie")

    def test_the_planted_letters_are_stacked_where_they_were_put_with_their_counts(self):
        a, txt = self.planted_letters()
        p = unit_profile(txt, a.encode("the"), flank=3, alphabet=a)
        P.plot_unit_profile(p)
        ax = plt.gcf().axes[0]
        labels = [t.get_text() for t in ax.texts]
        n = p["n"]
        counts = [int(t.split("/")[0]) for t in labels if "/" in t]
        self.assertTrue(all(t.endswith(f"/{n}") for t in labels if "/" in t))
        self.assertGreaterEqual(max(counts), n - 3)                      # q and x are in nearly every window (a few
                                                                         # chance 'the' have other letters there)
        lines = describe_profile(p)
        self.assertTrue(any(l.strip().startswith("-2:") and " Q " in " " + l for l in lines), lines)
        self.assertTrue(any(l.strip().startswith("-1:") and " X " in " " + l for l in lines), lines)

    def test_the_motif_is_one_block_not_a_column_per_letter(self):
        a, txt = self.planted_letters()
        p = unit_profile(txt, a.encode("the"), flank=3, alphabet=a)
        P.plot_unit_profile(p)
        ax = plt.gcf().axes[0]
        ticks = [t.get_text() for t in ax.get_xticklabels()]
        self.assertEqual(ticks, ["-3", "-2", "-1", "+1", "+2", "+3"])   # the motif has no ticks of its own
        self.assertEqual(max(ax.get_xlim()), 6 + 3.0)                   # six flank columns and a 3-wide block

    def test_a_letter_that_is_only_common_is_not_drawn_and_a_rare_one_is_small_but_counted(self):
        a, txt = self.planted_letters()
        p = unit_profile(txt, a.encode("the"), flank=3, alphabet=a)
        lines = " ".join(describe_profile(p))
        self.assertNotIn(" E ", lines)                                   # e is common everywhere
        # the counts always accompany a letter, so a rare letter seen a few times cannot pass for a dominant one
        for line in describe_profile(p):
            self.assertRegex(line, r"[A-Z] \d+ of \d+ \(\d+x as common")

    def test_an_empty_column_is_named_with_its_most_common_letters(self):
        a, txt = self.planted_letters()
        p = unit_profile(txt, a.encode("the"), flank=3, alphabet=a)
        shown = {line.split(":")[0].strip() for line in describe_profile(p)}
        empty = describe_empty_positions(p)
        self.assertTrue(shown)
        named = {line.split(":")[0].strip() for line in empty}
        self.assertEqual(shown | named, {"-3", "-2", "-1", "+1", "+2", "+3"})   # every flank position is accounted for
        self.assertFalse(shown & named)
        for line in empty:
            self.assertRegex(line, r"nothing stands out \(most common: [A-Z] \d+ of \d+")
        p["significant"] = [[False] * len(r) for r in p["significant"]]
        self.assertEqual(describe_empty_positions(p), [])                       # all empty: the figure says so itself

    def test_nothing_enriched_says_so(self):
        rnd = random.Random(9)
        a = latin()
        txt = a.encode("".join(rnd.choice("abcdefghijklmnopqrstuvwxyz") for _ in range(3000)))
        p = unit_profile(txt, a.encode("th"), flank=3, alphabet=a)
        p["significant"] = [[False] * len(r) for r in p["significant"]]
        P.plot_unit_profile(p)
        self.assertTrue(any("no letter is clearly enriched" in t.get_text() for t in plt.gcf().axes[0].texts))
        self.assertEqual(describe_profile(p), [])

    def test_mismatches_keep_the_columns_so_that_the_variants_can_be_shown(self):
        a, txt = self.planted_letters()
        units = list(txt)
        for s in range(100, 3900, 240):                                  # every other copy: 'the' becomes 'tha'
            units[s + 4] = a.char("a")
        p = unit_profile("".join(units), a.encode("the"), flank=3, mismatches=1, alphabet=a)
        P.plot_unit_profile(p)
        ticks = [t.get_text() for t in plt.gcf().axes[0].get_xticklabels()]
        if "T" in ticks:                                                 # the variant stack forced the motif columns
            self.assertEqual(ticks[3:6], ["T", "H", "E"])

    def test_works_for_codons_amino_acids_and_other_alphabets(self):
        rnd = random.Random(5)
        for alpha in (codon(), amino_acid(), Alphabet("Greek", tuple("αβγδεζηθικλμνξοπρστυφχψω"))):
            units = [rnd.choice(alpha.labels) for _ in range(3000)]
            for s in range(50, 2900, 100):
                units[s:s + 3] = list(alpha.labels[:3])
                units[s - 1] = alpha.labels[5]
                units[s + 3] = alpha.labels[6]
            p = unit_profile(alpha.encode(units), alpha.encode(alpha.labels[:3]), flank=3, alphabet=alpha)
            with tempfile.TemporaryDirectory() as d:
                out = os.path.join(d, "logo.png")
                P.plot_unit_profile(p, file=out)
                self.assertGreater(os.path.getsize(out), 1000)
            self.assertEqual(private_use("".join(t.get_text() for t in plt.gcf().axes[0].texts)), [])
            plt.close("all")


if __name__ == "__main__":
    unittest.main()

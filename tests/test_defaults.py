"""Discovery settings that belong to an alphabet (phase 5)."""

import argparse
import contextlib
import csv
import io
import json
import os
import random
import tempfile
import unittest

import rna_lexis.menu as M
from rna_lexis.alphabet import (
    NUCLEOTIDE_SETTINGS, Alphabet, codon, default_settings, get_alphabet, latin, nucleotide,
)
from rna_lexis.genetic_code import amino_acid
from rna_lexis.test_cli import _resolve_settings, build_parser, main


class DefaultSettingsTests(unittest.TestCase):
    def test_nucleotide_text_keeps_its_settings(self):
        self.assertEqual(default_settings(None), {"minxmlen": 7, "mincorelen": 6,
                                                  "candidate_min_len": 5, "mincount": 2})
        self.assertEqual(default_settings(nucleotide()), default_settings(None))

    def test_built_in_alphabets_carry_calibrated_settings(self):
        self.assertEqual(default_settings(codon()), {"minxmlen": 5, "mincorelen": 3, "candidate_min_len": 3, "mincount": 2})
        self.assertEqual(default_settings(amino_acid()), default_settings(codon()))
        self.assertEqual(default_settings(latin()), {"minxmlen": 6, "mincorelen": 4, "candidate_min_len": 4, "mincount": 2})
        # the shortest xmotif is the shortest core plus the 2-unit extension the calibration planted
        for a in (codon(), latin(), amino_acid()):
            s = default_settings(a)
            self.assertGreaterEqual(s["minxmlen"], s["candidate_min_len"] + 1)

    def test_other_alphabets_get_a_rule_or_their_own_values(self):
        greek = Alphabet(name="Greek", symbols=tuple("αβγδεζηθικλμν"))
        self.assertEqual(default_settings(greek)["candidate_min_len"], 3)
        small = Alphabet(name="five", symbols=tuple("acgtm"))
        self.assertEqual(default_settings(small), default_settings(None))          # 8 symbols or fewer
        own = Alphabet(name="Greek", symbols=tuple("αβγδεζηθικλμν"), defaults={"minxmlen": 8, "mincount": 3, "junk": 1})
        s = default_settings(own)
        self.assertEqual((s["minxmlen"], s["mincount"], s["candidate_min_len"]), (8, 3, 3))
        self.assertNotIn("junk", s)

    def test_defaults_are_saved_with_an_alphabet_file(self):
        a = Alphabet(name="x", symbols=tuple("abcdefghij"), defaults={"minxmlen": 9})
        with tempfile.TemporaryDirectory() as d:
            b = get_alphabet(a.save(os.path.join(d, "x.json")))
        self.assertEqual(default_settings(b)["minxmlen"], 9)

    def test_a_session_saved_before_the_calibration_gets_the_built_in_values(self):
        old = latin().to_dict()
        old["defaults"], old["unit_name"] = {}, "unit"
        a = Alphabet.from_dict(old)
        self.assertEqual(default_settings(a), default_settings(latin()))
        self.assertTrue(a.same_symbols(latin()))
        # an alphabet with the same name but other symbols does not borrow them
        odd = Alphabet.from_dict({"name": "Latin letters", "symbols": list("abc")})
        self.assertEqual(odd.defaults, {})
        self.assertFalse(odd.same_symbols(latin()))

    def test_same_symbols_ignores_defaults_and_display(self):
        a = codon()
        b = a._replace(defaults={}, unit_name="triplet", display_sep="-")
        self.assertTrue(a.same_symbols(b))
        self.assertFalse(a.same_symbols(a.with_extras(" ")))


class MenuDefaultsTests(unittest.TestCase):
    def fresh(self):
        return {'minxmlen': 7, 'maxxmlen': 60, 'mincorelen': 6, 'mincount': 2}

    def test_switching_alphabets_changes_untouched_settings(self):
        d = self.fresh()
        M._apply_alphabet_defaults(d, codon())
        self.assertEqual((d['minxmlen'], d['mincorelen'], d['mincount'], d['candidate_min_len']), (5, 3, 2, 3))
        M._apply_alphabet_defaults(d, latin())
        self.assertEqual((d['minxmlen'], d['mincorelen'], d['candidate_min_len']), (6, 4, 4))
        M._apply_alphabet_defaults(d, None)                         # back to nucleotides
        self.assertEqual((d['minxmlen'], d['mincorelen'], d['mincount'], d['candidate_min_len']), (7, 6, 2, 5))

    def test_a_setting_the_user_changed_is_kept(self):
        d = self.fresh()
        M._apply_alphabet_defaults(d, codon())
        d['minxmlen'] = 9                                           # the user's choice in "Change setting"
        M._apply_alphabet_defaults(d, latin())
        self.assertEqual(d['minxmlen'], 9)
        self.assertEqual(d['mincorelen'], 4)                        # untouched, so it follows the alphabet

    def test_nucleotide_sessions_are_unchanged(self):
        d = self.fresh()
        M._apply_alphabet_defaults(d, None)
        self.assertEqual({k: d[k] for k in ('minxmlen', 'maxxmlen', 'mincorelen', 'mincount')},
                         {'minxmlen': 7, 'maxxmlen': 60, 'mincorelen': 6, 'mincount': 2})

    def test_rank_core_motifs_prompt_uses_the_alphabet_default(self):
        rnd = random.Random(3)
        A = codon()
        units = [rnd.choice(A.symbols) for _ in range(300)]
        core = list(A.symbols[10:15])
        for k, s in enumerate((20, 100, 180, 260)):
            units[s:s + 7] = core + list(A.symbols[30:32] if k % 2 == 0 else A.symbols[40:42])
        txt = A.encode(units)
        defvals = {'minxmlen': 5, 'maxxmlen': 60, 'mincorelen': 3, 'mincount': 2}
        with contextlib.redirect_stdout(io.StringIO()):
            strs = M.parsedata(txt, defvals, A)
        with tempfile.TemporaryDirectory() as d:
            fp = os.path.join(d, "s.fa")
            M.open_file_with_default_software = lambda *a, **k: None
            answers = iter([""] * 8)
            M.safe_input = lambda prompt="": next(answers, "")
            with contextlib.redirect_stdout(io.StringIO()):
                M.statistical_core_input(fp, txt, strs, A)
            rows = list(csv.DictReader(open(os.path.join(d, "s_ranked_cores_markov.csv"), encoding="utf-8")))
        lengths = {int(r["length"]) for r in rows}
        self.assertEqual(min(lengths), 3)                           # codon default, not the nucleotide 5


class CliDefaultsTests(unittest.TestCase):
    def parse(self, *argv):
        args = build_parser().parse_args(argv)
        _resolve_settings(args)
        return args

    def test_nucleotide_defaults_are_unchanged(self):
        a = self.parse("rank-cores", "--input", "x", "--output", "y")
        self.assertEqual((a.min_xmotif_len, a.candidate_min_len, a.min_count), (7, 5, 2))

    def test_alphabet_defaults(self):
        for name, expect in (("codon", (5, 3, 2)), ("latin", (6, 4, 2)), ("amino-acid", (5, 3, 2))):
            a = self.parse("rank-cores", "--input", "x", "--output", "y", "--alphabet", name)
            self.assertEqual((a.min_xmotif_len, a.candidate_min_len, a.min_count), expect, name)

    def test_explicit_options_win(self):
        a = self.parse("rank-cores", "--input", "x", "--output", "y", "--alphabet", "codon",
                       "--min-xmotif-len", "8", "--candidate-min-len", "6", "--min-count", "3")
        self.assertEqual((a.min_xmotif_len, a.candidate_min_len, a.min_count), (8, 6, 3))

    def test_score_exact_gets_the_count_default_too(self):
        a = self.parse("score-exact", "--input", "x", "--motifs", "m", "--output", "y")
        self.assertEqual(a.min_count, 2)

    def test_rank_cores_finds_a_planted_core_with_no_length_options(self):
        A = codon()
        rnd = random.Random(7)
        units = [rnd.choice(A.symbols) for _ in range(400)]
        core = list(A.symbols[10:16])
        for k, s in enumerate((20, 100, 180, 260, 340)):
            units[s:s + 8] = core + list(A.symbols[30:32] if k % 2 == 0 else A.symbols[40:42])
        with tempfile.TemporaryDirectory() as d:
            fasta = os.path.join(d, "cds.fasta")
            open(fasta, "w").write(">c\n" + "".join(units) + "\n")
            out = os.path.join(d, "o.csv")
            with contextlib.redirect_stdout(io.StringIO()):
                main(["rank-cores", "--input", fasta, "--output", out, "--alphabet", "codon"])
            rows = list(csv.DictReader(open(out, encoding="utf-8")))
        supported = [r["motif"] for r in rows if r["statistically_supported"] == "True"]
        core_s = " ".join(c.upper() for c in core)
        self.assertTrue(any(m in core_s for m in supported), supported)


if __name__ == "__main__":
    unittest.main()

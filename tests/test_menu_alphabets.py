"""The interactive menu over alphabets other than nucleotides (phase 4).

The menu has no terminal in the tests: ``safe_input`` is replaced by a script of
answers, plots are drawn in-process, and the system file opener is stubbed.
"""

import contextlib
import io
import json
import os
import random
import tempfile
import unittest
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.basedatatypes as pb

import rna_lexis.io as IO
import rna_lexis.menu as M
import rna_lexis.plots as P
from rna_lexis.alphabet import Alphabet, codon, latin

A = codon()
C = A.symbols
CORE = list(C[10:17])                                  # 7 codons
CORE_NT = "".join(CORE)                                # what a user types for it


def private_use(text):
    return [c for c in text if 0xE000 <= ord(c) <= 0xF8FF]


class Script:
    """Answers for safe_input, in order; blank when they run out."""

    def __init__(self, answers):
        self.answers = list(answers)
        self.calls = 0

    def __call__(self, prompt=""):
        self.calls += 1
        if self.calls > 80:
            raise RuntimeError("too many prompts (a loop?)")
        return self.answers.pop(0) if self.answers else ""


def codon_text(n=400, seed=1):
    rnd = random.Random(seed)
    units = [rnd.choice(C) for _ in range(n)]
    ext = [list(C[30:32]), list(C[40:42])]
    for k, s in enumerate((20, 100, 180, 260, 340)):
        units[s:s + 9] = CORE + ext[k % 2]
    units[180 + 3] = C[50]                              # a one-codon variant in one copy
    return A.encode(units)


class MenuCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        self.fp = os.path.join(self.dir, "session.fa")
        self.figs = []
        self._saved = {}

        def stub(obj, name, value):
            self._saved[(id(obj), name)] = (obj, name, getattr(obj, name))
            setattr(obj, name, value)
        stub(plt, "show", lambda *a, **k: None)
        stub(pb.BaseFigure, "show", lambda fig, *a, **k: self.figs.append(fig))
        for mod in (M, P, IO):
            for name in ("open_file_with_default_software", "open_pdf"):
                if hasattr(mod, name):
                    stub(mod, name, lambda *a, **k: None)
        stub(M, "_spawn_plot", lambda fn, *a, **k: fn(*a, **k))
        self.txt = codon_text()
        self.defvals = {'minxmlen': 4, 'maxxmlen': 60, 'mincorelen': 3, 'mincount': 2, 'pwr': 1.2,
                        'clr': False, 'reload': False, 'datadir': self.dir}
        with contextlib.redirect_stdout(io.StringIO()):
            self.strs = M.parsedata(self.txt, self.defvals, A)

    def tearDown(self):
        for obj, name, value in self._saved.values():
            setattr(obj, name, value)
        plt.close("all")
        self.tmp.cleanup()

    def run_fn(self, fn, answers):
        """Run fn with scripted answers; return (stdout, {new file: text})."""
        before = set(os.listdir(self.dir))
        self.figs.clear()
        M.safe_input = Script(answers)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            fn()
        files = {}
        for name in sorted(set(os.listdir(self.dir)) - before):
            path = os.path.join(self.dir, name)
            if name.endswith((".csv", ".html")):
                with open(path, encoding="utf-8") as fh:
                    files[name] = fh.read()
            else:
                files[name] = None
        return buf.getvalue(), files

    def figtext(self):
        return json.dumps([json.loads(f.to_json()) for f in self.figs], ensure_ascii=False)


class QueryTests(unittest.TestCase):
    def test_nucleotide_queries_pass_through(self):
        self.assertEqual(M._query_to_internal("GATTACA", None), "GATTACA")
        self.assertEqual(M._show("gattaca", None), "gattaca")

    def test_codon_queries_are_read_like_the_text(self):
        got = M._query_to_internal("atgaaa", A)
        self.assertEqual(A.decode_list(got), ["atg", "aaa"])

    def test_wildcard_is_kept_between_units(self):
        got = M._query_to_internal("atg.aaa", A, wildcard=True)
        self.assertEqual(got.count("."), 1)
        self.assertEqual(A.decode_list(got.replace(".", "")), ["atg", "aaa"])

    def test_letters(self):
        self.assertEqual(M._query_to_internal("the", latin()), "the")

    def test_alphabet_survives_a_session_round_trip(self):
        data = json.loads(json.dumps({"alphabet": A.to_dict()}))
        self.assertEqual(M._alphabet_of(data).to_dict(), A.to_dict())
        self.assertIsNone(M._alphabet_of({}))
        self.assertIs(M._alphabet_of({"alphabet": A}), A)


class WildcardQueryTests(unittest.TestCase):
    """``.`` in a query: one unit, and ``...`` is one codon, as in nucleotide text."""

    def setUp(self):
        self._safe = M.safe_input

    def tearDown(self):
        M.safe_input = self._safe

    def count(self, query):
        return len(M._query_to_internal(query, A, wildcard=True).replace(".", "")), \
            M._query_to_internal(query, A, wildcard=True).count(".")

    def test_three_dots_are_one_codon(self):
        self.assertEqual(self.count("act ... ... aga"), (2, 2))
        self.assertEqual(self.count("act...aga"), (2, 1))
        self.assertEqual(self.count("act......aga"), (2, 2))         # six dots: two codons

    def test_one_dot_is_one_codon_and_other_runs_are_one_per_dot(self):
        self.assertEqual(self.count("act . . aga"), (2, 2))
        self.assertEqual(self.count("act .. aga"), (2, 2))
        self.assertEqual(self.count("act .... aga"), (2, 4))

    def test_letters_keep_one_dot_per_letter(self):
        self.assertEqual(M._query_to_internal("t...e", latin(), wildcard=True), "t...e")

    def test_without_the_wildcard_flag_dots_are_dropped_as_before(self):
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(M._query_to_internal("act ... aga", A), A.encode(["act", "aga"]))

    def search(self, txt, query, rate="3"):
        M.safe_input = Script([query, rate, ""])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            M.search_input(txt, alphabet=A)
        return buf.getvalue()

    def test_search_with_mutations_takes_a_codon_wildcard(self):
        rnd = random.Random(2)
        pool = ["atg", "aaa", "ccc", "ggt", "tat", "gca", "tct"]
        units = [rnd.choice(pool) for _ in range(300)]
        units[40:48] = ["act", "ttt", "ggg", "aga", "aaa", "cca", "tct", "tat"]          # exact, other codons in the gap
        units[100:108] = ["act", "ccc", "aaa", "aga", "aaa", "ccc", "tct", "gca"]        # two codons differ outside the gap
        txt = A.encode(units)
        out = self.search(txt, "act ... ... aga aaa cca tct tat")
        self.assertIn("positions: [40]", out)
        self.assertIn("pos   100", out)
        self.assertIn("dist=2", out)
        self.assertIn("L=8", out)
        approx_line = [l for l in out.splitlines() if "pos   100" in l][0]
        self.assertIn("act ccc aaa aga aaa", approx_line)                        # the gap is shown plainly, not marked
        self.assertEqual(approx_line.count("\x1b[31m"), 2)                        # only the two real differences are red

    def test_a_query_of_only_wildcards_is_cancelled(self):
        self.assertEqual(self.search(A.encode(["atg", "aaa"] * 10), "... ..."), "")


class ReadingTests(unittest.TestCase):
    def setUp(self):
        self._show_menu, self._safe, self._open = M.show_menu, M.safe_input, M.openFile

    def tearDown(self):
        M.show_menu, M.safe_input, M.openFile = self._show_menu, self._safe, self._open

    def choose(self, data, choice, answers=()):
        M.show_menu = lambda *a, **k: choice
        M.safe_input = Script(answers)
        with contextlib.redirect_stdout(io.StringIO()):
            return M._choose_reading(data)

    def data(self, txt, **kw):
        d = {'file_path': 'x.fa', 'txt': txt, 'is_rna': M.contains_only_rna(txt), 'txtb': txt, 'dir': ''}
        d.update(kw)
        return d

    def test_automatic_leaves_nucleotide_text_alone(self):
        d = self.data("acgtacgtacgt")
        self.assertIs(self.choose(d, 1), d)

    def test_automatic_reads_other_text_as_letters(self):
        out = self.choose(self.data("thequickbrownfox"), 1)
        self.assertEqual(out["alphabet"].name, "Latin letters")
        self.assertEqual(out["txt"], "thequickbrownfox")
        self.assertEqual(out["reading"]["units"], 16)

    def test_codons_with_a_frame(self):
        out = self.choose(self.data("catgaaattt"), 2, answers=["1"])
        self.assertEqual(out["alphabet"].name, "codon")
        self.assertEqual(out["alphabet"].decode_list(out["txt"]), ["atg", "aaa", "ttt"])
        self.assertEqual(out["reading"]["frame"], 1)
        self.assertEqual(out["reading"]["dropped"], 1)
        self.assertFalse(out["is_rna"])

    def test_an_ambiguity_code_skips_the_codon_and_is_reported(self):
        out = self.choose(self.data("atgnnnaaattt"), 2, answers=[""])
        self.assertEqual(out["alphabet"].decode_list(out["txt"]), ["atg", "aaa", "ttt"])
        self.assertEqual(out["reading"]["skipped"], 1)

    def coding(self, codons, utr5="gcgcc", utr3="ggcgc"):
        return utr5 + "atg" + "gaa" * (codons - 1) + "taa" + utr3

    def test_a_coding_transcript_is_read_from_its_start_codon_by_default(self):
        txt = self.coding(400)
        out = self.choose(self.data(txt), 2, answers=[""])
        self.assertEqual(out["alphabet"].decode_list(out["txt"])[0], "atg")
        self.assertEqual(out["alphabet"].decode_list(out["txt"])[-1], "taa")
        self.assertEqual(out["reading"]["units"], 401)
        self.assertEqual(out["reading"]["orf"], {"start": 5, "end": 5 + 3 * 401, "codons": 400})
        out = self.choose(self.data(txt), 3, answers=["o"])         # also when translated
        self.assertEqual(out["alphabet"].decode(out["txt"]), "m" + "e" * 399 + "*")

    def test_a_frame_can_still_be_chosen_for_a_coding_transcript(self):
        out = self.choose(self.data(self.coding(400)), 2, answers=["0"])
        self.assertNotIn("orf", out["reading"])
        self.assertEqual(out["reading"]["frame"], 0)

    def test_a_short_open_reading_frame_is_not_the_default_and_a_tiny_one_is_refused(self):
        mid = self.data(self.coding(150))                          # could be chance: not the default
        out = self.choose(mid, 2, answers=[""])
        self.assertNotIn("orf", out["reading"])
        out = self.choose(mid, 2, answers=["o"])                   # but it can be asked for
        self.assertEqual(out["reading"]["orf"]["codons"], 150)
        tiny = self.data(self.coding(20))
        self.assertIsNone(self.choose(tiny, 2, answers=["o"]))     # not a coding region: back to the menu
        self.assertIsNotNone(self.choose(tiny, 2, answers=["1"]))  # a frame is still allowed

    def test_codons_translated_to_amino_acids(self):
        out = self.choose(self.data("atggcctaaggg"), 3, answers=[""])
        self.assertEqual(out["alphabet"].name, "amino acid")
        self.assertEqual(out["alphabet"].decode(out["txt"]), "ma*g")
        self.assertEqual(out["reading"]["translated"], "standard")
        self.assertEqual(out["reading"]["units"], 4)
        # the frame is applied before translating
        out = self.choose(self.data("catggcc"), 3, answers=["1"])
        self.assertEqual(out["alphabet"].decode(out["txt"]), "ma")

    def test_amino_acids_and_letters(self):
        out = self.choose(self.data("mkvlaxyg"), 4)           # x is not an amino acid
        self.assertEqual(out["alphabet"].name, "amino acid")
        self.assertEqual(out["alphabet"].decode(out["txt"]), "mkvlayg")
        self.assertEqual(out["reading"]["skipped"], 1)
        out = self.choose(self.data("thequick"), 5)
        self.assertEqual(out["alphabet"].name, "Latin letters")

    def test_alphabet_from_a_file(self):
        greek = Alphabet(name="Greek", symbols=tuple("αβγδ"))
        with tempfile.TemporaryDirectory() as d:
            path = greek.save(os.path.join(d, "greek.json"))
            M.openFile = lambda **k: path
            out = self.choose(self.data("abcd", raw="αβγδα βγ"), 6)
        self.assertEqual(out["alphabet"].name, "Greek")
        self.assertEqual(out["alphabet"].decode(out["txt"]), "αβγδαβγ")

    def test_cancel_and_sessions(self):
        self.assertIsNone(self.choose(self.data("acgt"), 7))
        self.assertIsNone(self.choose(self.data("acgt"), 0))
        session = self.data("acgt", corelist=["acg"], xmotifs=["acg"])
        self.assertIs(self.choose(session, 2), session)             # a saved session keeps its reading
        self.assertEqual(M._choose_reading(-1), -1)


class CachedTranscriptTests(unittest.TestCase):
    """A transcript that already has a saved session: the saved one is used, or the sequence is read again."""

    def setUp(self):
        self._safe, self._show = M.safe_input, M.show_menu
        self.tmp = tempfile.TemporaryDirectory()
        root = self.tmp.name
        self.gene = os.path.join(root, "BRCA")
        os.mkdir(self.gene)
        self.txt = "gcgcc" + "atg" + "gaa" * 399 + "taa" + "ggcgc"
        session = {"file_path": os.path.join(self.gene, "ENST1_X"), "txt": self.txt, "txtb": self.txt,
                   "is_rna": True, "corelist": ["gaa"], "xmotifs": ["gaagaa"], "dir": self.gene,
                   "alphabet": None, "reading": None}
        with open(os.path.join(self.gene, "ENST1_X.json"), "w", encoding="utf-8") as fh:
            json.dump(session, fh)

    def tearDown(self):
        M.safe_input, M.show_menu = self._safe, self._show
        self.tmp.cleanup()

    def load(self, answers):
        M.safe_input = Script(answers)
        with contextlib.redirect_stdout(io.StringIO()):
            return M.choose_enst(self.gene)

    def test_the_saved_session_is_the_default(self):
        got = self.load(["ENST1", ""])
        self.assertTrue(got.get("corelist"))
        self.assertNotIn("reread", got)

    def test_it_can_be_read_again_as_codons_and_is_saved_under_another_name(self):
        got = self.load(["ENST1", "2"])
        self.assertFalse(got.get("corelist"))
        self.assertTrue(got["reread"])
        M.show_menu = lambda *a, **k: 2
        M.safe_input = Script([""])
        with contextlib.redirect_stdout(io.StringIO()):
            out = M._choose_reading(got)
        self.assertEqual(out["reading"]["orf"]["codons"], 400)
        self.assertTrue(out["file_path"].endswith("ENST1_X_codon"))
        self.assertNotIn("reread", out)
        self.assertTrue(os.path.exists(os.path.join(self.gene, "ENST1_X.json")))      # the first one is untouched


class GatingTests(unittest.TestCase):
    def test_nucleotide_only_items_are_refused_for_other_alphabets(self):
        M.safe_input = Script([""])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertFalse(M._nucleotide_ok(A, False, "The sequence logo"))
            self.assertTrue(M._nucleotide_ok(None, True, "The sequence logo"))
            self.assertFalse(M._nucleotide_ok(None, False, "Hairpin search"))
        self.assertIn("needs a nucleotide sequence", buf.getvalue())
        self.assertIn("codon", buf.getvalue())


class MenuOverCodonsTests(MenuCase):
    def check(self, out, files, expect=None):
        stdout, = (out,)
        self.assertEqual(private_use(stdout), [], stdout)
        for name, text in files.items():
            if text is not None:
                self.assertEqual(private_use(text), [], name)
        if expect:
            self.assertIn(expect, stdout)

    def test_find_matches(self):
        out, files = self.run_fn(lambda: M.find_match_input(self.txt, self.strs, alphabet=A), [CORE_NT[:9] + ".", ""])
        self.check(out, files, A.show(A.encode(CORE[:3])))
        self.assertIn("Found 5 matches", out)          # "." is any one codon, so the variant matches too
        # exactly four codons reach the planted variant, which then no longer matches
        out, files = self.run_fn(lambda: M.find_match_input(self.txt, self.strs, alphabet=A), [CORE_NT[:12], ""])
        self.assertIn("Found 4 matches", out)

    def test_wildcard_search(self):
        q = CORE_NT[:9] + "." + CORE_NT[12:15]          # "any codon" where the variant sits
        out, files = self.run_fn(lambda: M.find_match_input(self.txt, self.strs, alphabet=A), [q, ""])
        self.check(out, files)
        self.assertIn("Found 5 matches", out)

    def test_search_with_mutations_marks_the_changed_codon(self):
        out, files = self.run_fn(lambda: M.search_input(self.txt, alphabet=A), [CORE_NT, "3", ""])
        self.check(out, files)
        self.assertIn(A.show(A.encode(CORE)), out)
        marked = A.show(A.encode(CORE[:3] + [C[50]] + CORE[4:]), mask=[0, 0, 0, 1, 0, 0, 0])
        self.assertIn(marked.split(" ")[3], out)                    # the changed codon, in upper case

    def test_motif_extensions(self):
        out, files = self.run_fn(lambda: M.motif_extensions_input(self.txt, self.fp, alphabet=A),
                                 [CORE_NT, "1", "", "", "", ""])
        self.check(out, files, "pos=")
        self.assertTrue(any(n.startswith("extensions_") for n in files))

    def test_print_core_and_coverage(self):
        out, files = self.run_fn(lambda: M.print_core_input(self.txt, self.strs, alphabet=A), [CORE_NT[:12], ""])
        self.check(out, files)
        out, files = self.run_fn(lambda: M.txt_coverage_input(self.txt, self.strs, alphabet=A), [CORE_NT[:12], "", ""])
        self.check(out, files, "length 4")

    def test_neighbors_plots_and_export(self):
        for fn in (M.neighbors_input, M.neighbors_condensed_input):
            out, files = self.run_fn(lambda: fn(self.fp, self.txt, self.strs, A), [CORE_NT, "", "", "", "", "", ""])
            self.check(out, files)
            self.assertEqual(private_use(self.figtext()), [])
            self.assertIn(A.show(A.encode(CORE)), self.figtext())
        out, files = self.run_fn(lambda: M.neighbors_condensed_export_input(self.fp, self.txt, self.strs, A),
                                 [CORE_NT, "", "", "", ""])
        self.check(out, files)
        name = next(n for n in files if n.endswith("_regions.csv"))
        self.assertIn(A.show(A.encode(CORE)), files[name])

    def test_kmers_and_coverage_plots(self):
        out, files = self.run_fn(lambda: M.kmers_input(self.fp, self.txt, self.defvals, A), ["3", "1", "", "z"])
        self.check(out, files)
        out, files = self.run_fn(lambda: M.coverage_input(self.txt, self.strs, self.fp, A), ["", "", "cov"])
        self.check(out, files)
        self.assertIn("cov.png", files)

    def test_sequence_hits_with_the_detailed_view(self):
        out, files = self.run_fn(lambda: M.sequence_hits_input(self.fp, self.txt, A),
                                 [CORE_NT, "", "3", "", "", "y", "", ""])
        self.check(out, files)
        html = files["session_detail.html"]
        self.assertIn("codons", html)
        self.assertNotIn("nucleotides", html)
        self.assertEqual(private_use(self.figtext()), [])

    def test_arcs(self):
        out, files = self.run_fn(lambda: M.self_similarity_arcs_input(self.txt, self.fp, A),
                                 [CORE_NT[:12], "", "", "", "arcs", ""])      # seed, range, rate, arc mode, file
        self.check(out, files)
        self.assertIn("arcs.png", files)
        self.assertEqual(plt.gca().get_xlabel(), "Position (codons)")

    def test_statistical_core_ranking(self):
        out, files = self.run_fn(lambda: M.statistical_core_input(self.fp, self.txt, self.strs, A),
                                 ["3", "", "", "", "", "", ""])
        self.check(out, files, "Candidates tested")
        text = files["session_ranked_cores_markov.csv"]
        first = text.splitlines()[1].split(",")[0]
        self.assertRegex(first, r"^[ACGT]{3}( [ACGT]{3})*$")          # decoded, upper-case codons

    def test_mutation_family(self):
        out, files = self.run_fn(lambda: M.mutation_family_input(self.fp, self.txt, self.strs, A),
                                 ["1", CORE_NT, "3", "", "", "", "", ""])
        self.check(out, files, "Motifs summarized")
        self.assertIn("hamming_variant_count", files["session_mutation_family_tests.csv"])

    def test_gapped_motif(self):
        out, files = self.run_fn(lambda: M.gapped_motif_input(self.fp, self.txt, A),
                                 [CORE_NT[:6], CORE_NT[12:18], "", "", ""])
        self.check(out, files, "observed=")
        self.assertIn(A.show(A.encode(CORE[:2]), upper=True) + "[gap:", out)

    def test_alignment_has_no_nucleotide_blast_score(self):
        out, files = self.run_fn(lambda: M.print_alignment_score(self.txt, A), ["20", "100", "20", "1", ""])
        self.check(out, files, "GLOBAL ALIGNMENT")
        self.assertIn("not computed", out)
        self.assertNotIn("E-value (vs", out)

    def test_markov_kmers_and_spacing(self):
        out, files = self.run_fn(lambda: M.markov_kmer_input(self.txt, self.fp, A), ["3", "", "", ""])
        self.check(out, files)
        out, files = self.run_fn(lambda: M.spacing_test_input(self.txt, A), [CORE_NT, "", ""])
        self.check(out, files, "Spacing analysis")
        self.assertIn("codons", out)
        self.assertNotIn(" nt", out)
        out, files = self.run_fn(lambda: M.batch_spacing_test_input(self.fp, self.txt, self.strs, A), ["", "", ""])
        self.check(out, files)

    def test_unit_profile_replaces_the_logo(self):
        out, files = self.run_fn(lambda: M.profile_input(self.txt, self.strs, self.fp, A),
                                 [CORE_NT[:12], "2", "", "", "", "prof", "", "y"])      # style: blank = the letter logo
        self.check(out, files, "occurrence(s) used")
        self.assertIn("prof.png", files)
        self.assertIn(" of ", out)                         # the enriched units are listed with their counts
        name = next(n for n in files if "_profile_" in n)
        self.assertEqual(private_use(files[name]), [])
        self.assertIn("count@-2", files[name])
        # the other style
        out, files = self.run_fn(lambda: M.profile_input(self.txt, self.strs, self.fp, A),
                                 [CORE_NT[:12], "2", "", "", "2", "heat", "", "n"])
        self.check(out, files, "occurrence(s) used")
        self.assertIn("heat.png", files)
        # a motif that is not in the text
        out, files = self.run_fn(lambda: M.profile_input(self.txt, self.strs, self.fp, A), ["gggggggggggg", ""])
        self.assertIn("not found", out)
        self.assertEqual(files, {})

    def test_session_summary(self):
        with contextlib.redirect_stdout(io.StringIO()):
            path = IO.init_summary(self.fp, self.strs["xmotifs"], self.strs["corelist"], self.txt, alphabet=A)
        text = open(path, encoding="utf-8").read()
        self.assertEqual(private_use(text), [])
        self.assertIn("q_markov", text)
        self.assertRegex(text.splitlines()[1].split(",")[0], r"^[a-z]{3}( [a-z]{3})*$")

    def test_parsedata_reports_units(self):
        d = dict(self.defvals, maxxmlen=5)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            M.parsedata(self.txt, d, A)
        self.assertIn("codons", buf.getvalue())
        self.assertNotIn(" nt)", buf.getvalue())


class ComparisonTests(MenuCase):
    def test_a_comparison_text_is_read_over_the_same_alphabet(self):
        data = {'file_path': 'b.fa', 'txt': "atgaaattt" * 3, 'is_rna': True, 'dir': ''}
        with contextlib.redirect_stdout(io.StringIO()):
            out = M._comparison_in_alphabet(data, A, {'frame': 0})
        self.assertEqual(out['alphabet'].name, "codon")
        self.assertEqual(out['alphabet'].decode_list(out['txt'])[:3], ["atg", "aaa", "ttt"])

    def test_a_translated_session_reads_the_comparison_text_the_same_way(self):
        aa = M.amino_acid()
        data = {'file_path': 'b.fa', 'txt': "atggcctaa", 'is_rna': True, 'dir': ''}
        with contextlib.redirect_stdout(io.StringIO()):
            out = M._comparison_in_alphabet(data, aa, {'frame': 0, 'translated': 'standard'})
        self.assertEqual(aa.decode(out['txt']), "ma*")
        self.assertEqual(out['reading']['translated'], 'standard')

    def test_a_session_with_another_alphabet_is_refused(self):
        other = {'file_path': 'b', 'txt': 'abc', 'corelist': ['ab'], 'xmotifs': ['ab'],
                 'alphabet': latin().to_dict()}
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertIsNone(M._comparison_in_alphabet(other, A, None))
        self.assertIn("cannot be compared", buf.getvalue())
        same = dict(other, alphabet=A.to_dict())
        self.assertIs(M._comparison_in_alphabet(same, A, None), same)


if __name__ == "__main__":
    unittest.main()

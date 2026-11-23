"""Robustness: wrong input, wrong alphabet, empty input, bad files, and what must not be trusted.

Each test names a way a user (or a file from someone else) can go wrong and checks that the
package says so, in words, and does not crash, hang, or quietly give an answer about something else.
"""

import contextlib
import io
import json
import os
import random
import tempfile
import unittest
import warnings

import rna_lexis.io as IO
import rna_lexis.menu as M
from rna_lexis import statistical as S
from rna_lexis import tokenizers as T
from rna_lexis.alphabet import MAX_ALPHABET_FILE_BYTES, Alphabet, codon, get_alphabet, latin, nucleotide
from rna_lexis.genetic_code import amino_acid
from rna_lexis.test_cli import main as cli_main


def caught(fn, *args, **kwargs):
    """Run fn; return (result, [warning categories], [warning messages])."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = fn(*args, **kwargs)
    return result, [x.category for x in w], [str(x.message) for x in w]


# ----------------------------------------------------------------------
# reading text: unknown letters, gaps, nothing at all
# ----------------------------------------------------------------------
class WrongKindOfInputTests(unittest.TestCase):
    def test_what_is_not_text_is_refused_with_a_clear_message(self):
        for bad in (None, 123, 4.5, b"acgt", [None], ["acg", 7], object()):
            with self.subTest(bad=bad), self.assertRaises(TypeError) as cm:
                T.tokenize(bad, nucleotide())
            self.assertIn("string", str(cm.exception))

    def test_the_alphabet_must_be_an_alphabet(self):
        with self.assertRaises(TypeError):
            T.tokenize("acgt", "nucleotide")

    def test_bad_options_are_named(self):
        for kw in ({"gap": 0}, {"frame": -1}, {"frame": 3}, {"width": 0}, {"width": "3"}, {"mode": "zzz"}):
            with self.subTest(kw=kw), self.assertRaises(ValueError):
                T.tokenize("atgaaa", codon(), **kw)

    def test_reading_codons_at_a_width_their_symbols_do_not_have_is_an_error_not_silence(self):
        with self.assertRaises(ValueError) as cm:
            T.tokenize("atgaaa", codon(), width=2)
        self.assertIn("3", str(cm.exception))


class NothingToReadTests(unittest.TestCase):
    def test_empty_and_blank_text_read_as_nothing_without_a_warning(self):
        for text in ("", "   \n\t ", []):
            for a in (nucleotide(), codon(), latin()):
                with self.subTest(text=text, alphabet=a.name):
                    tok, cats, _ = caught(T.tokenize, text, a)
                    self.assertEqual(tok.text, "")
                    self.assertEqual(tok.n_units, 0)
                    self.assertEqual(cats, [])

    def test_text_with_no_symbol_of_the_alphabet_warns_that_the_alphabet_may_be_wrong(self):
        tok, cats, msgs = caught(T.tokenize, "the quick brown fox jumps over the lazy dog", nucleotide())
        self.assertIn(T.AlphabetMismatchWarning, cats)
        self.assertIn(T.UnknownSymbolWarning, cats)
        self.assertTrue(any("right alphabet" in m for m in msgs))

    def test_the_statistics_return_nothing_for_nothing_and_do_not_fail(self):
        for text, a in (("", None), ("a", None), ("acg", None), ("", latin()), ("a", latin()), ("a" * 50, None)):
            with self.subTest(text=text[:5], alphabet=getattr(a, "name", None)):
                motif = "a" if a is not None else "ac"
                self.assertIsInstance(S.score_exact_motifs(text, [motif], alphabet=a), list)
                self.assertEqual(S.rank_core_candidates(text, [], alphabet=a), [])
                self.assertIsInstance(S.mutation_family_tests(text, [motif], alphabet=a), list)
                self.assertIsInstance(S.score_gapped_motif(text, "a", "c", alphabet=a), dict)
                if text.count("a") < 2:
                    self.assertEqual(S.unit_profile(text, "a", alphabet=a or nucleotide())["n"], 0)

    def test_a_motif_longer_than_the_text_is_not_found(self):
        self.assertEqual(S.score_exact_motifs("acg", ["acgtacgtacg"])[0]["exact_count"], 0)
        self.assertEqual(S.hamming_family_positions("acgtacgtacg", "acg", 3), [])


class GapTests(unittest.TestCase):
    """Whitespace is formatting, unless the alphabet declares a space; then a line break is a gap."""

    def test_without_a_declared_space_whitespace_is_formatting(self):
        self.assertEqual(T.tokenize("to be\nor", latin()).decode(), "tobeor")

    def test_with_a_declared_space_a_line_break_between_words_is_one_gap(self):
        a = latin().with_extras(" ")
        for text, expect in (("to be\nor", "to be or"), ("to be\r\nor not", "to be or not"),
                             ("x \n y", "x y"), ("tab\there", "tab here"), ("\n\nlead", " lead")):
            with self.subTest(text=text):
                self.assertEqual(T.tokenize(text, a).decode(), expect)

    def test_a_run_of_declared_spaces_is_read_as_written(self):
        a = latin().with_extras(" ")
        self.assertEqual(T.tokenize("a  b", a).decode(), "a  b")

    def test_in_fixed_mode_a_line_break_is_still_formatting(self):
        a = codon().with_extras(" ")
        self.assertEqual(a.decode_list(T.tokenize("atg\naaa", a).text), ["atg", "aaa"])
        self.assertEqual(a.decode_list(T.tokenize("atg aaa", a).text), ["atg", " ", "aaa"])

    def test_in_words_mode_whitespace_separates_and_cannot_be_a_symbol(self):
        a = Alphabet("syl", ("ba", "ka"), mode="words")
        self.assertEqual(T.tokenize("ba  ka\nba", a).decode(sep=" "), "ba ka ba")
        with self.assertRaises(ValueError) as cm:
            Alphabet("syl", ("ba", "ka"), extras=(" ",), mode="words")
        self.assertIn("separator", str(cm.exception))

    def test_segments_are_kept_apart_by_the_boundary(self):
        tok = T.tokenize(["ab", "cd"], latin(), gap=2)
        self.assertEqual(tok.text, "ab__cd")


# ----------------------------------------------------------------------
# the alphabet does not fit the text
# ----------------------------------------------------------------------
class AlphabetMismatchTests(unittest.TestCase):
    DNA = "acgtacgtttgacca" * 4

    def mismatch(self, text, alphabet, **kw):
        _, cats, msgs = caught(T.tokenize, text, alphabet, **kw)
        return [m for c, m in zip(cats, msgs) if c is T.AlphabetMismatchWarning] if len(cats) == len(msgs) else []

    def test_dna_read_with_latin_letters_or_amino_acids_is_flagged(self):
        for a in (latin(), amino_acid()):
            with self.subTest(alphabet=a.name):
                found = self.mismatch(self.DNA, a)
                self.assertEqual(len(found), 1)
                self.assertIn("DNA or RNA", found[0])
                self.assertIn("nucleotide", found[0])

    def test_english_read_as_nucleotides_and_protein_read_as_codons_are_flagged(self):
        self.assertEqual(len(self.mismatch("the quick brown fox jumps over the lazy dog", nucleotide())), 1)
        self.assertEqual(len(self.mismatch("MKVLAAGIVGLLLAQMKVLAAGIVGLLLAQ", codon())), 1)

    def test_words_read_with_the_wrong_separator_say_to_check_the_separator(self):
        a = Alphabet("syl", ("ba", "ka", "ta"), mode="words")
        found = self.mismatch("ba-ka-ta-ba-ka-ta-ba-ka-ta", a)
        self.assertEqual(len(found), 1)
        self.assertIn("separated by whitespace", found[0])

    def test_right_matches_are_not_flagged(self):
        english = "the quick brown fox jumps over the lazy dog and keeps running far away"
        self.assertEqual(self.mismatch(english, latin()), [])
        self.assertEqual(self.mismatch(self.DNA, nucleotide()), [])
        self.assertEqual(self.mismatch(self.DNA, codon()), [])
        self.assertEqual(self.mismatch("MKVLAAGIVGLLLAQMKVLAAGIVGLLLAQMKV", amino_acid()), [])

    def test_a_small_alphabet_that_includes_the_dna_letters_is_not_flagged(self):
        # IUPAC-style: 4 bases plus ambiguity codes; reading DNA with it is what it is for
        iupac = Alphabet("iupac", tuple("acgtnrykmswbdhv"))
        self.assertEqual(self.mismatch(self.DNA, iupac), [])

    def test_a_short_text_or_one_with_a_few_stray_symbols_is_not_flagged(self):
        self.assertEqual(self.mismatch("acgt", latin()), [])                    # too short to tell
        self.assertEqual(self.mismatch("acgtacgtacgtacg1acgtacgtacgt", nucleotide()), [])

    def test_the_menu_shows_the_warning_when_reading(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            internal, reading = M._read_text_as(self.DNA, latin())
        self.assertIsNotNone(internal)
        self.assertIn("looks like DNA", out.getvalue())

    def test_the_cli_shows_the_warning(self):
        with tempfile.TemporaryDirectory() as d:
            fa = os.path.join(d, "dna.fa")
            open(fa, "w").write(">x\n" + self.DNA + "\n")
            motifs = os.path.join(d, "m.txt")
            open(motifs, "w").write("acg\n")
            err, out = io.StringIO(), io.StringIO()
            with contextlib.redirect_stderr(err), contextlib.redirect_stdout(out):
                cli_main(["score-exact", "--input", fa, "--motifs", motifs, "--output", os.path.join(d, "o.csv"),
                          "--alphabet", "latin"])
            self.assertIn("looks like DNA", err.getvalue())


# ----------------------------------------------------------------------
# alphabets and alphabet files
# ----------------------------------------------------------------------
class AlphabetInputTests(unittest.TestCase):
    def test_wrong_kinds_of_field_are_named(self):
        bad = [
            dict(name=None, symbols=("a",)),
            dict(name="x", symbols=("a",), width="1"),
            dict(name="x", symbols=("a",), width=0),
            dict(name="x", symbols=("a",), width=True),
            dict(name="x", symbols=(1, 2)),
            dict(name="x", symbols=("a",), case_sensitive="yes"),
            dict(name="x", symbols=("a",), aliases={"u": 3}),
            dict(name="x", symbols=("a",), unit_name=7),
            dict(name="x", symbols=("a",), display_sep=7),
            dict(name="x", symbols=tuple("abcdefghij"), defaults={"minxmlen": "zz"}),
            dict(name="x", symbols=tuple("abcdefghij"), defaults={"mincount": -3}),
            dict(name="x", symbols=tuple("abcdefghij"), defaults={"mincorelen": 2.5}),
            dict(name="x", symbols=("a", "b"), codes={"a": ""}),
        ]
        for kw in bad:
            with self.subTest(kw=kw), self.assertRaises(ValueError):
                Alphabet(**kw)

    def test_an_unknown_default_key_is_kept(self):
        a = Alphabet("x", tuple("abcdefghij"), defaults={"note": "free text"})
        self.assertEqual(a.defaults["note"], "free text")


class AlphabetFileTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def file(self, content, name="a.json"):
        path = os.path.join(self.tmp.name, name)
        mode = "wb" if isinstance(content, bytes) else "w"
        with open(path, mode) as fh:
            fh.write(content)
        return path

    def test_a_file_that_is_not_an_alphabet_is_a_value_error_naming_the_file(self):
        cases = {
            "empty": "", "not json": "{{{", "a list": "[1, 2]", "null": "null", "a number": "5",
            "binary": b"\x00\xff\xfe\x01", "no symbols": '{"name": "x"}', "no name": '{"symbols": ["a"]}',
            "symbols null": '{"name": "x", "symbols": null}', "symbol an int": '{"name": "x", "symbols": [1]}',
            "aliases a list": '{"name": "x", "symbols": ["a"], "aliases": [1]}',
            "width text": '{"name": "x", "symbols": ["a"], "width": "a"}',
            "deep": "[" * 100000,
        }
        for label, content in cases.items():
            path = self.file(content)
            with self.subTest(case=label), self.assertRaises(ValueError) as cm:
                get_alphabet(path)
            self.assertIn("a.json", str(cm.exception))

    def test_a_file_that_is_too_large_is_refused_before_it_is_read(self):
        path = self.file("x")
        with open(path, "w") as fh:
            fh.write(" " * (MAX_ALPHABET_FILE_BYTES + 1))
        with self.assertRaises(ValueError) as cm:
            Alphabet.load(path)
        self.assertIn("too large", str(cm.exception))

    def test_a_good_file_still_loads_and_round_trips(self):
        a = Alphabet("greek", tuple("αβγδ"), extras=(" ",))
        path = a.save(os.path.join(self.tmp.name, "greek.json"))
        self.assertTrue(get_alphabet(path).same_symbols(a))

    def test_get_alphabet_takes_a_name_a_display_name_or_a_path(self):
        for name, expect in (("latin", "Latin letters"), ("Latin letters", "Latin letters"), ("  LATIN  ", "Latin letters"),
                             ("amino-acid", "amino acid"), ("amino acid", "amino acid"), ("Amino_Acid", "amino acid"),
                             ("codons", "codon"), ("nucleotide", "nucleotide")):
            with self.subTest(name=name):
                self.assertEqual(get_alphabet(name).name, expect)

    def test_get_alphabet_with_something_else_says_what_it_wants(self):
        for bad in (None, 5, "", "nonsense", self.tmp.name):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                get_alphabet(bad)


# ----------------------------------------------------------------------
# numbers and options
# ----------------------------------------------------------------------
class ParameterTests(unittest.TestCase):
    TXT = "acgtacgtacgtacgtttgacca" * 10

    def test_a_probability_outside_zero_one_is_refused(self):
        for alpha in (0, -0.1, 5, float("nan"), float("inf"), "0.05", None):
            with self.subTest(alpha=alpha), self.assertRaises(ValueError):
                S.score_exact_motifs(self.TXT, ["acg"], alpha=alpha)
        S.score_exact_motifs(self.TXT, ["acg"], alpha=1.0)

    def test_a_mutation_rate_above_one_is_refused(self):
        with self.assertRaises(ValueError):
            S.mutation_family_tests(self.TXT, ["acgtac"], mutr=2)
        with self.assertRaises(ValueError):
            S.mutation_family_tests(self.TXT, ["acgtac"], M=-1)

    def test_the_markov_order_must_be_a_whole_number_of_at_least_zero(self):
        for order in (-1, "a", 1.5, True):
            with self.subTest(order=order), self.assertRaises(ValueError):
                S.score_exact_motifs(self.TXT, ["acg"], markov_order=order)
        S.score_exact_motifs(self.TXT, ["acg"], markov_order=2)

    def test_an_order_with_more_contexts_than_letters_warns(self):
        _, cats, msgs = caught(S.score_exact_motifs, "acgt" * 5, ["ac"], markov_order=6)
        self.assertIn(UserWarning, cats)
        self.assertTrue(any("contexts" in m for m in msgs))

    def test_a_zero_pseudocount_is_refused_because_it_divides_by_zero(self):
        with self.assertRaises(ValueError):
            S.MarkovBackground("acgtacgt", pseudocount=0)

    def test_counts_and_lengths_must_be_whole_and_in_range(self):
        for kw in ({"min_count": 0}, {"min_count": 1.5}, {"min_xmotif_type_support": -1}):
            with self.subTest(kw=kw), self.assertRaises(ValueError):
                S.score_exact_motifs(self.TXT, ["acg"], **kw)
        with self.assertRaises(ValueError):
            S.enumerate_core_candidates(["acgtac"], candidate_min_len=9, candidate_max_len=5)

    def test_a_profile_window_that_would_fill_memory_is_refused(self):
        with self.assertRaises(ValueError) as cm:
            S.unit_profile(self.TXT, "acg", flank=999_999_999_999)
        self.assertIn("flank", str(cm.exception))
        for kw in ({"flank": -1}, {"mismatches": -1}, {"alpha": 2}):
            with self.subTest(kw=kw), self.assertRaises(ValueError):
                S.unit_profile(self.TXT, "acg", **kw)

    def test_the_menu_helper_takes_a_good_number_and_falls_back_on_a_bad_one(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(M._number("", 5, int, low=1), 5)
            self.assertEqual(M._number(" 7 ", 5, int, low=1), 7)
            self.assertEqual(M._number("7.0", 5, int, low=1), 7)
            self.assertEqual(M._number("2.5", 1.0, float, low=0), 2.5)
            for bad in ("abc", "nan", "inf", "1e999", "-3", "0", "3.5", "\x00", "5-3", "1,2"):
                self.assertEqual(M._number(bad, 5, int, low=1), 5, bad)
            self.assertEqual(M._number("0", 9.0, float, low=0, low_open=True), 9.0)
            self.assertEqual(M._number("99", 1.0, float, low=0, high=20), 1.0)
        self.assertIn("using the default", out.getvalue())


# ----------------------------------------------------------------------
# what is written into spreadsheets
# ----------------------------------------------------------------------
class SpreadsheetSafetyTests(unittest.TestCase):
    def test_csv_safe_only_changes_text_that_a_spreadsheet_would_run(self):
        for value in ("=1+1", "+cmd", "-cmd", "@SUM(A1)", "\tx", "=HYPERLINK(\"http://x\",\"y\")"):
            self.assertEqual(S.csv_safe(value), "'" + value)
        for value in ("ACGT", "-3.2", "+5", "1e-9", "", 12, 0.5, None, "ATG AAA"):
            self.assertEqual(S.csv_safe(value), value)

    def test_a_symbol_label_cannot_start_a_formula_in_a_results_file(self):
        evil = Alphabet("evil", ("=cmd", "+x", "-y", "@z", "ok"), mode="words")
        tok = T.tokenize("=cmd +x -y @z ok " * 6, evil)
        rows = S.score_exact_motifs(tok.text, [tok.text[:2]], alphabet=evil)
        with tempfile.TemporaryDirectory() as d:
            path = S.write_rows_csv(os.path.join(d, "r.csv"), rows, alphabet=evil)
            text = open(path, encoding="utf-8").read()
        self.assertNotIn("\n=CMD", text)
        self.assertIn("'=CMD", text)

    def test_nucleotide_results_are_written_exactly_as_before(self):
        rows = S.score_exact_motifs("acgtacgtacgtacgt", ["acg"])
        with tempfile.TemporaryDirectory() as d:
            path = S.write_rows_csv(os.path.join(d, "r.csv"), rows)
            text = open(path, encoding="utf-8").read()
        self.assertIn("ACG,", text)
        self.assertNotIn("'", text)

    def test_the_unit_profile_table_is_protected_too(self):
        evil = Alphabet("evil", ("=cmd", "ok"), mode="words")
        tok = T.tokenize("=cmd ok ok =cmd ok ok =cmd ok ok =cmd ok ok", evil)
        profile = S.unit_profile(tok.text, tok.text[:1], flank=1, alphabet=evil)
        with tempfile.TemporaryDirectory() as d:
            text = open(S.write_unit_profile_csv(os.path.join(d, "p.csv"), profile), encoding="utf-8").read()
        self.assertNotIn("\n=CMD", text)


# ----------------------------------------------------------------------
# files: FASTA, sessions, fetched identifiers
# ----------------------------------------------------------------------
class FastaTests(unittest.TestCase):
    def test_a_comment_line_before_the_header_does_not_make_the_headers_sequence(self):
        text = "; made by a tool\n; another comment\n>one gene\nacgt\nacgt\n>two\ntttt\n"
        headers, records = IO.split_fasta(text)
        self.assertEqual(headers, ["one gene", "two"])
        self.assertEqual([r.replace("\n", "") for r in records], ["acgtacgt", "tttt"])

    def test_plain_text_and_empty_text_are_not_fasta(self):
        for text in ("acgt\nacgt\n", "", "   \n", "; only a comment\n", "hello > world\n"):
            self.assertIsNone(IO.split_fasta(text), repr(text))

    def test_a_blank_line_and_a_comment_inside_a_record_are_dropped(self):
        headers, records = IO.split_fasta(">a\nac\n; note\n\ngt\n")
        self.assertEqual(records[0].replace("\n", ""), "acgt")

    def test_the_menu_loader_skips_a_leading_comment_and_all_headers(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "g.fa")
            open(path, "w").write("; comment, with letters acgt\n>ENST1 gene\nacgtacgt\n>ENST2\nacgtacgt\n")
            data = M._load_fasta_or_text(path)
        self.assertEqual(data["txt"], "acgtacgtacgtacgt")
        self.assertEqual(data["gene_name"], "ENST1")

    def test_a_bare_header_line_has_no_name_and_does_not_crash(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "g.fa")
            open(path, "w").write(">\nacgt\n")
            data = M._load_fasta_or_text(path)
        self.assertEqual(data["txt"], "acgt")
        self.assertNotIn("gene_name", data)


class SessionFileTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def path(self, content):
        p = os.path.join(self.tmp.name, "s.json")
        open(p, "w", encoding="utf-8").write(content if isinstance(content, str) else json.dumps(content))
        return p

    def good(self, **extra):
        a = latin()
        d = {"file_path": "x.fa", "txt": "abcabc_abc", "txtb": "", "is_rna": False, "dir": "",
             "corelist": ["abc"], "xmotifs": ["abca"], "stats": {}, "alphabet": a.to_dict(), "reading": {}}
        d.update(extra)
        return d

    def test_load_session_says_what_is_wrong_and_returns_none(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertIsNone(IO.load_session(self.path("[1, 2]")))
            self.assertIsNone(IO.load_session(self.path("{{{")))
            self.assertIsNone(IO.load_session(os.path.join(self.tmp.name, "missing.json")))
            self.assertIsNone(IO.load_session(self.tmp.name))
        self.assertIn("not a session file", buf.getvalue())
        self.assertIn("invalid format", buf.getvalue())
        self.assertIn("was not found", buf.getvalue())

    def test_a_consistent_session_has_no_problem(self):
        self.assertIsNone(IO.session_problem(self.good()))
        self.assertIsNone(IO.session_problem(self.good(alphabet=None, txt="acgt", corelist=[], xmotifs=[])))

    def test_each_kind_of_damage_is_named(self):
        cases = {
            "not an object": [],
            "no txt": {k: v for k, v in self.good().items() if k != "txt"},
            "txt a number": self.good(txt=5),
            "corelist not strings": self.good(corelist=[1, 2]),
            "corelist not a list": self.good(corelist="abc"),
            "alphabet broken": self.good(alphabet={"name": "x"}),
            "characters outside the alphabet": self.good(txt="abc"),
            "core outside the alphabet": self.good(corelist=["ab1"]),
        }
        for label, data in cases.items():
            with self.subTest(case=label):
                self.assertIsInstance(IO.session_problem(data), str)

    def test_the_menu_refuses_a_damaged_session_instead_of_looping_on_it(self):
        p = self.path(self.good(txt="abc"))
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertIsNone(M._load_session_checked(p))
        self.assertIn("Cannot use the session", buf.getvalue())
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertIsNotNone(M._load_session_checked(self.path(self.good())))


class FetchInputTests(unittest.TestCase):
    """Identifiers are put into a web address, so anything that could change the address is refused before any request."""

    def test_transcript_ids_that_are_not_ids_are_refused(self):
        for bad in ("ENST00000526704/../../x", "ENST1?type=x", "ENST1#frag", "ENST1 2", "ENST", "ENSTabc",
                    "ENST1.2.3", 5, None):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                IO.fetch_enst_cdna(bad)

    def test_genomic_ranges_that_are_not_ranges_are_refused(self):
        for args in (("chr10;x=1", 1, 5), ("10", 1, 5), ("chr10", -1, 5), ("chr10", 5, 5), ("chr10", 9, 2),
                     ("chr10", "1", 5), ("chr10", 1, 5, "hg38;x"), ("chr10", 0, IO.MAX_RANGE + 1),
                     ("chr10", True, 5), ("chr 10", 1, 5)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                IO.fetch_genomic_range(*args)

    def test_accessions_that_are_not_accessions_are_refused(self):
        for bad in ("EH38E1482203;x", "a b", "", None, "EH/../x"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                IO.fetch_encode_ccre(bad)

    def test_an_oversized_response_is_refused(self):
        class Big:
            def read(self, n):
                return b"x" * n
        with self.assertRaises(ValueError):
            IO._read_capped(Big())


# ----------------------------------------------------------------------
# the menu and the command line
# ----------------------------------------------------------------------
class MenuInputTests(unittest.TestCase):
    def test_a_pattern_that_is_not_a_pattern_is_reported_not_raised(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            for bad in ("[a", "(", "a{2,1}", "*", "a(b"):
                self.assertEqual(M._query_to_internal(bad, None), "", bad)
            self.assertEqual(M._query_to_internal("a.c", None), "a.c")
            self.assertEqual(M._query_to_internal("ac[gt]", None), "ac[gt]")
        self.assertIn("not a valid pattern", out.getvalue())

    def test_find_match_with_a_bad_pattern_says_so_and_asks_again(self):
        answers = iter(["[acg", "", "gta", ""])      # the pattern, "press any key", a good query, "press Enter"

        def scripted(prompt=""):
            try:
                return next(answers)
            except StopIteration:
                raise M.EOFSignal()          # what Ctrl+D does: the dialog ends instead of asking for ever
        saved = M.safe_input
        M.safe_input = scripted
        try:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                try:
                    M.find_match_input("acgtacgtacgt", {}, minlen=3, alphabet=None)
                except M.EOFSignal:
                    pass
        finally:
            M.safe_input = saved
        self.assertIn("not a valid pattern", out.getvalue())
        self.assertIn("Found 2 matches for gta", out.getvalue())

    def test_many_pairs_asks_first(self):
        saved = M.safe_input
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertTrue(M._confirm_many_pairs(100))                       # 4,950 pairs: no question
                M.safe_input = lambda prompt="": "n"
                self.assertFalse(M._confirm_many_pairs(3000))                     # 4.5 million pairs
                M.safe_input = lambda prompt="": "y"
                self.assertTrue(M._confirm_many_pairs(3000))
        finally:
            M.safe_input = saved


class CommandLineInputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def p(self, name):
        return os.path.join(self.tmp.name, name)

    def w(self, name, content):
        open(self.p(name), "w", encoding="utf-8").write(content)
        return self.p(name)

    def run_cli(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            cli_main(argv)
        return out.getvalue(), err.getvalue()

    def expect_error(self, argv, *words):
        with self.assertRaises(SystemExit) as cm:
            self.run_cli(argv)
        message = str(cm.exception)
        self.assertTrue(message.startswith("error:"), message)
        for word in words:
            self.assertIn(word, message)

    def test_a_missing_input_file_is_one_line_not_a_traceback(self):
        motifs = self.w("m.txt", "acg\n")
        self.expect_error(["score-exact", "--input", self.p("nope.fa"), "--motifs", motifs,
                           "--output", self.p("o.csv")], "nope.fa")

    def test_nucleotide_input_with_no_nucleotides_is_an_error_that_suggests_what_to_try(self):
        text = self.w("t.txt", "hello world is mine\n")      # no a, c, g, t or u in it
        motifs = self.w("m.txt", "acg\n")
        self.expect_error(["score-exact", "--input", text, "--motifs", motifs, "--output", self.p("o.csv")],
                          "no nucleotides", "--alphabet latin")

    def test_mostly_non_nucleotide_input_is_processed_as_before_with_a_warning_on_stderr(self):
        text = self.w("t.txt", "acgtacgtacgt " + "hello world this is mostly not dna " * 3 + "\n")
        motifs = self.w("m.txt", "acg\n")
        _, err = self.run_cli(["score-exact", "--input", text, "--motifs", motifs, "--output", self.p("o.csv")])
        self.assertIn("Is this nucleotide text", err)
        self.assertTrue(os.path.getsize(self.p("o.csv")) > 0)

    def test_a_fasta_file_that_begins_with_a_comment_is_read_without_its_headers(self):
        fa = self.w("c.fa", "; comment acgtacgtacgt\n>ENST header acgt\n" + "gattaca" * 6 + "\n")
        motifs = self.w("m.txt", "gattaca\n")
        self.run_cli(["score-exact", "--input", fa, "--motifs", motifs, "--output", self.p("o.csv"),
                      "--enrichment-threshold", "1"])
        self.assertIn(",6,", open(self.p("o.csv")).read().splitlines()[1])

    def test_bad_motif_files_are_named(self):
        fa = self.w("s.fa", ">x\n" + "gattaca" * 4 + "\n")
        out = self.p("o.csv")
        self.expect_error(["score-exact", "--input", fa, "--motifs", self.w("bad.json", "{{{"), "--output", out],
                          "not valid JSON")
        self.expect_error(["score-exact", "--input", fa, "--motifs", self.w("n.json", "5"), "--output", out],
                          "should hold")
        self.expect_error(["score-exact", "--input", fa, "--motifs", self.w("d.json", '{"a": 1}'),
                           "--motif-column", "zzz", "--output", out], "no field 'zzz'")
        self.expect_error(["score-exact", "--input", fa, "--motifs", self.w("t.csv", "seq,n\ngat,1\n"),
                           "--motif-column", "zzz", "--output", out], "no column 'zzz'")
        self.expect_error(["score-exact", "--input", fa, "--motifs", self.p("absent.txt"), "--output", out],
                          "absent.txt")

    def test_bad_options_are_one_line_errors_in_every_mode(self):
        fa = self.w("s.fa", ">x\n" + "gattaca" * 4 + "\n")
        motifs = self.w("m.txt", "gattaca\n")
        for extra in (["--alpha", "5"], ["--markov-order", "-1"], ["--min-count", "0"]):
            with self.subTest(extra=extra):
                self.expect_error(["score-exact", "--input", fa, "--motifs", motifs, "--output", self.p("o.csv")] + extra)
        self.expect_error(["score-exact", "--input", fa, "--motifs", motifs, "--output", self.p("o.csv"),
                           "--alphabet", "latin", "--frame", "2"], "frame")
        self.expect_error(["score-exact", "--input", fa, "--motifs", motifs, "--output", self.p("o.csv"),
                           "--alphabet", "nonsense"], "nonsense")


# ----------------------------------------------------------------------
# the answer to "does an interval contain this one" without comparing against all of them
# ----------------------------------------------------------------------
class IntervalIndexTests(unittest.TestCase):
    def test_same_answer_as_looking_at_every_interval(self):
        rnd = random.Random(8)
        for trial in range(60):
            intervals = []
            for _ in range(rnd.randint(0, 40)):
                a = rnd.randint(0, 200)
                intervals.append((a, a + rnd.randint(1, 60), "x"))
            index = S._IntervalIndex(intervals)
            for _ in range(80):
                s = rnd.randint(0, 220)
                e = s + rnd.randint(1, 30)
                self.assertEqual(index.contains(s, e), S._inside_any_interval(s, e, intervals))


# ----------------------------------------------------------------------
# Clear workspace keeps inputs; the logo leaves no intermediate file behind
# ----------------------------------------------------------------------
class ClearWorkspaceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.d = self.tmp.name
        for name in ("gene.fa", "gene.json", "other.fasta", "OTHER2.FNA", "third.fas", "gene_test_init.csv",
                     "plot.png", "notes.txt", "gene_logo.pdf"):
            open(os.path.join(self.d, name), "w").write("x")
        self.saved = M.safe_input

    def tearDown(self):
        M.safe_input = self.saved

    def run_clear(self, answers):
        it = iter(answers)
        M.safe_input = lambda prompt="": next(it, "")
        with contextlib.redirect_stdout(io.StringIO()) as out:
            done = M._clear_workspace(self.d, os.path.join(self.d, "gene.fa"))
        return done, out.getvalue()

    def test_inputs_sessions_and_other_fasta_files_are_kept_and_listed_as_kept(self):
        done, text = self.run_clear(["YES", ""])
        self.assertTrue(done)
        self.assertEqual(sorted(os.listdir(self.d)), ["OTHER2.FNA", "gene.fa", "gene.json", "other.fasta", "third.fas"])
        kept_part = text.split("KEPT")[1]
        for name in ("other.fasta", "OTHER2.FNA", "third.fas", "gene.json", "gene.fa"):
            self.assertIn(name, kept_part)

    def test_anything_but_yes_deletes_nothing(self):
        before = sorted(os.listdir(self.d))
        for answer in ("", "yes", "Y", "no"):
            done, _ = self.run_clear([answer, ""])
            self.assertFalse(done)
        self.assertEqual(sorted(os.listdir(self.d)), before)


class LogoIntermediateTests(unittest.TestCase):
    def test_the_sequences_go_to_a_temporary_folder_and_nothing_is_left_beside_the_output(self):
        import rna_lexis.plots as P
        seen = {}

        def fake_run(cmd, **kwargs):
            seen["cmd"] = cmd
            fa = cmd[cmd.index("-f") + 1]
            seen["fa_existed"] = os.path.isfile(fa)
            seen["fa"] = fa
            open(cmd[cmd.index("-o") + 1], "w").write("pdf")
        saved = (P.shutil.which, P.subprocess.run, P.open_pdf)
        P.shutil.which = lambda name: "/usr/bin/weblogo"
        P.subprocess.run = fake_run
        P.open_pdf = lambda *a, **k: None
        try:
            with tempfile.TemporaryDirectory() as d:
                mine = os.path.join(d, "mine.fa")
                open(mine, "w").write("my own file")
                out = os.path.join(d, "mine")                 # the user names the output like their input file
                with contextlib.redirect_stdout(io.StringIO()):
                    P.plot_logo("acg", 2, "ttacgaatacgcc" * 3, out, fmt="pdf")
                self.assertEqual(sorted(os.listdir(d)), ["mine.fa", "mine.pdf"])
                self.assertEqual(open(mine).read(), "my own file")      # not overwritten
        finally:
            P.shutil.which, P.subprocess.run, P.open_pdf = saved
        self.assertTrue(seen["fa_existed"])                              # WebLogo did get its input
        self.assertFalse(os.path.exists(seen["fa"]))                     # and it is gone now


if __name__ == "__main__":
    unittest.main()

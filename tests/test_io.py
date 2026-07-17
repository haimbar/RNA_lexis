import json
import os
import tempfile
import unittest
from unittest import mock

from rna_lexis.io import (
    _find_valid_sessions,
    _read_summary_hash,
    _summary_inputs_hash,
    _write_summary_hash,
    fetch_encode_ccre,
    fetch_genomic_range,
    is_valid_session,
    load_prefs,
    load_session,
    open_file_with_default_software,
    read_text,
    save_prefs,
    save_session,
)


class SessionRoundTripTests(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.base = os.path.join(self.tmpdir.name, "sess")
        # Includes a non-ASCII character to guard against the encoding
        # regression fixed in 0.1.12 (open() calls without encoding='utf-8'
        # break on Windows' default locale codepage for non-ASCII content).
        self.data = {
            "file_path": self.base,
            "txt": "acgtaccgu",
            "gene_name": "café-locus-ü",
            "is_rna": False,
            "corelist": ["acgt"],
            "xmotifs": ["acgtaccgu"],
            "dir": self.tmpdir.name,
            "stats": {},
            "txtb": "",
        }

    def test_round_trip_preserves_non_ascii_content(self):
        save_session(self.base, self.data)
        loaded = load_session(self.base + ".json")
        self.assertEqual(loaded, self.data)

    def test_is_valid_session_accepts_a_freshly_saved_session(self):
        save_session(self.base, self.data)
        result = is_valid_session(self.base + ".json")
        self.assertIsNotNone(result)
        self.assertEqual(result.txt, "acgtaccgu")

    def test_is_valid_session_rejects_malformed_json(self):
        bad_path = os.path.join(self.tmpdir.name, "bad.json")
        with open(bad_path, "w", encoding="utf-8") as f:
            f.write("not valid json")
        self.assertIsNone(is_valid_session(bad_path))

    def test_is_valid_session_rejects_missing_required_keys(self):
        incomplete_path = os.path.join(self.tmpdir.name, "incomplete.json")
        with open(incomplete_path, "w", encoding="utf-8") as f:
            json.dump({"txt": "acgt"}, f)
        self.assertIsNone(is_valid_session(incomplete_path))

    def test_load_session_returns_none_for_missing_file(self):
        self.assertIsNone(load_session(os.path.join(self.tmpdir.name, "nope.json")))

    def test_find_valid_sessions_lists_only_valid_ones(self):
        save_session(self.base, self.data)
        with open(os.path.join(self.tmpdir.name, "bad.json"), "w", encoding="utf-8") as f:
            f.write("not json")
        found = _find_valid_sessions(self.tmpdir.name)
        self.assertEqual([name for name, _ in found], ["sess.json"])


class ReadTextTests(unittest.TestCase):
    def test_strips_non_alphabetic_characters_and_lowercases(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "t.txt")
            with open(path, "w", encoding="utf-8") as f:
                f.write("ACGT\nacgu123!!")
            self.assertEqual(read_text(path), "acgtacgu")

    def test_keepspacing_collapses_separators_to_single_space(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "t.txt")
            with open(path, "w", encoding="utf-8") as f:
                f.write("acgt\n\nacgu")
            self.assertEqual(read_text(path, keepspacing=True), "acgt acgu")

    def test_raises_value_error_for_missing_file(self):
        with self.assertRaises(ValueError):
            read_text("/nonexistent/path/does/not/exist.txt")


class PrefsTests(unittest.TestCase):
    def test_round_trip_via_mocked_prefs_path(self):
        with tempfile.TemporaryDirectory() as d:
            prefs_path = os.path.join(d, "prefs.json")
            with mock.patch("rna_lexis.io._prefs_path", return_value=prefs_path):
                self.assertEqual(
                    load_prefs(), {"default_data_dir": "", "last_used_dir": ""}
                )
                save_prefs({"default_data_dir": "x", "last_used_dir": "y"})
                self.assertEqual(
                    load_prefs(), {"default_data_dir": "x", "last_used_dir": "y"}
                )

    def test_falls_back_to_defaults_when_config_dir_unavailable(self):
        # Regression test for the 0.1.12 fix: _prefs_path()'s os.makedirs()
        # can fail (read-only/restrictive config dir) and load_prefs() must
        # not crash — it should fall back to defaults instead.
        with mock.patch(
            "rna_lexis.io._prefs_path", side_effect=OSError("permission denied")
        ):
            self.assertEqual(
                load_prefs(), {"default_data_dir": "", "last_used_dir": ""}
            )


class SummaryHashTests(unittest.TestCase):
    def test_hash_round_trips_and_changes_with_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            csv_path = os.path.join(d, "x.csv")
            h1 = _summary_inputs_hash(["a"], ["b"], "acgt", 1 / 6, 4)
            _write_summary_hash(csv_path, h1)
            self.assertEqual(_read_summary_hash(csv_path), h1)

            h2 = _summary_inputs_hash(["a"], ["b"], "acgtt", 1 / 6, 4)
            self.assertNotEqual(h1, h2)

    def test_read_summary_hash_returns_none_when_missing(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertIsNone(_read_summary_hash(os.path.join(d, "nope.csv")))


class OpenFileWithDefaultSoftwareTests(unittest.TestCase):
    def test_does_not_raise_for_a_nonexistent_path(self):
        # Regression test for the 0.1.12 fix: this used to propagate
        # CalledProcessError/FileNotFoundError uncaught on macOS/Linux.
        with tempfile.TemporaryDirectory() as d:
            try:
                open_file_with_default_software(os.path.join(d, "nope.qqzz"))
            except Exception as exc:  # pragma: no cover - failure path
                self.fail(f"open_file_with_default_software raised: {exc!r}")


class FetchGenomicRangeTests(unittest.TestCase):
    def test_returns_lowercase_sequence(self):
        with mock.patch("rna_lexis.io._fetch_ucsc_json",
                         return_value={"dna": "ACGTacgtNN"}) as mocked:
            seq = fetch_genomic_range("chr10", 100, 110)
        self.assertEqual(seq, "acgtacgtnn")
        url = mocked.call_args[0][0]
        self.assertIn("chr10", url)
        self.assertIn("start=100", url)
        self.assertIn("end=110", url)

    def test_raises_when_no_sequence_returned(self):
        with mock.patch("rna_lexis.io._fetch_ucsc_json", return_value={"dna": ""}):
            with self.assertRaises(ValueError):
                fetch_genomic_range("chr10", 100, 110)


class FetchEncodeCcreTests(unittest.TestCase):
    def test_finds_accession_and_returns_annotated_sequence(self):
        track_response = {"encodeCcreCombined": [{
            "name": "EH38E1482203", "chrom": "chr10",
            "chromStart": 78974544, "chromEnd": 78974893,
            "ccre": "dELS,CTCF-bound",
            "description": "distal enhancer-like signature",
        }]}
        seq_response = {"dna": "AAATGTGGGG"}

        def side_effect(url):
            return track_response if "track=encodeCcreCombined" in url else seq_response

        with mock.patch("rna_lexis.io._fetch_ucsc_json", side_effect=side_effect):
            annotation, seq = fetch_encode_ccre("EH38E1482203", hint_chrom="chr10")
        self.assertIn("EH38E1482203", annotation)
        self.assertIn("chr10:78974544-78974893", annotation)
        self.assertEqual(seq, "aaatgtgggg")

    def test_raises_when_accession_not_found(self):
        with mock.patch("rna_lexis.io._fetch_ucsc_json",
                         return_value={"encodeCcreCombined": []}):
            with self.assertRaises(ValueError):
                fetch_encode_ccre("EH38ENOTFOUND", hint_chrom="chr10")


if __name__ == "__main__":
    unittest.main()

"""File I/O, session management, and network fetch utilities."""

import atexit
import hashlib
import os
import os.path
import re
import json
import glob
import platform
import subprocess
import sys
import tempfile
import urllib.request
from contextlib import ExitStack
from dataclasses import dataclass, asdict
from importlib import resources as _resources
from typing import List, Optional

import pandas as pd
from platformdirs import user_config_dir
from scipy.stats import binom as _binom

from rna_lexis.algorithms import find_with_mutations
from rna_lexis.statistical import csv_safe, score_exact_motifs

TEST_SUMMARY_SUFFIX = "_test_init.csv"

EXAMPLE_DATASETS = ("NORAD_human", "NORAD_mouse")

# Keeps any temp-extracted resource files alive for the life of the process
# (only matters if the package is ever installed zipped, which pip does not
# do by default, but importlib.resources.as_file() requires a context
# manager either way).
_resource_files = ExitStack()
atexit.register(_resource_files.close)


def example_dataset_path(name: str) -> str:
    """Return a real filesystem path to a bundled example FASTA dataset.

    Args:
        name: One of `EXAMPLE_DATASETS` (currently ``'NORAD_human'`` or
              ``'NORAD_mouse'``), without the `.fasta` extension.

    Returns:
        Path to the `.fasta` file, valid for the life of the process.
        Works whether the package is installed from a wheel, an editable
        install, or run directly from a repo checkout.

    Raises:
        ValueError: If `name` is not a known bundled dataset.
    """
    if name not in EXAMPLE_DATASETS:
        raise ValueError(f"Unknown example dataset {name!r}; choose from {EXAMPLE_DATASETS}")
    resource = _resources.files("rna_lexis").joinpath("data", f"{name}.fasta")
    path = _resource_files.enter_context(_resources.as_file(resource))
    return str(path)


def init_summary_path(fn: str) -> str:
    """Return the test-build summary CSV path for a session base path."""
    return f'{fn}{TEST_SUMMARY_SUFFIX}'


# ---------------------------------------------------------------------------
# Persistent user preferences
# ---------------------------------------------------------------------------

def _prefs_path() -> str:
    """Return the path to the RNA_lexis preferences JSON file.

    The directory is the platform-appropriate user config location:
      Linux:   ~/.config/RNA_lexis/prefs.json
      macOS:   ~/Library/Application Support/RNA_lexis/prefs.json
      Windows: C:\\Users\\<user>\\AppData\\Roaming\\RNA_lexis\\prefs.json
    """
    config_dir = user_config_dir("RNA_lexis")
    os.makedirs(config_dir, exist_ok=True)
    return os.path.join(config_dir, "prefs.json")


def load_prefs() -> dict:
    """Load persistent user preferences from the platform config directory.

    Returns a dict with at least the keys ``'default_data_dir'`` and
    ``'last_used_dir'``, both defaulting to ``''`` when absent.
    """
    defaults = {'default_data_dir': '', 'last_used_dir': ''}
    try:
        path = _prefs_path()
    except Exception:
        return defaults
    if not os.path.isfile(path):
        return defaults
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return {**defaults, **json.load(f)}
    except Exception:
        return defaults


def save_prefs(prefs: dict) -> None:
    """Persist user preferences to the platform config directory.

    Args:
        prefs: Dict containing preference keys such as ``'default_data_dir'``
               and ``'last_used_dir'``.  Unknown keys are preserved.
    """
    try:
        with open(_prefs_path(), 'w', encoding='utf-8') as f:
            json.dump(prefs, f, indent=4)
    except Exception as e:
        print(f"Warning: could not save preferences: {e}")


_SESSION_REQUIRED_KEYS = frozenset(
    {'file_path', 'txt', 'is_rna', 'corelist', 'xmotifs', 'dir', 'stats'}
)


@dataclass
class SessionData:
    """Structured representation of a saved RNA_explore session."""
    file_path: str
    txt: str
    is_rna: bool
    corelist: List[str]
    xmotifs: List[str]
    dir: str
    stats: dict
    txtb: str = ''
    alphabet: Optional[dict] = None
    reading: Optional[dict] = None

    @classmethod
    def from_dict(cls, d: dict) -> 'SessionData':
        """Construct a SessionData from a raw dict (e.g. loaded from JSON)."""
        return cls(
            file_path=d['file_path'],
            txt=d['txt'],
            is_rna=d['is_rna'],
            corelist=d['corelist'],
            xmotifs=d['xmotifs'],
            dir=d['dir'],
            stats=d['stats'],
            txtb=d.get('txtb', ''),
            alphabet=d.get('alphabet'),
            reading=d.get('reading'),
        )


def save_session(filename, data, debug=False):
    """Serialize a session to a JSON file.

    Args:
        filename: Path without the .json extension; the file is written as
                  ``<filename>.json``.
        data:     A SessionData instance or a plain dict containing the session
                  state.  SessionData objects are converted via dataclasses.asdict.
        debug:    Reserved for future diagnostic output (currently unused).
    """
    if isinstance(data, SessionData):
        data = asdict(data)
    with open(f"{filename}.json", "w", encoding='utf-8') as json_file:
        json.dump(data, json_file, indent=4)


MAX_SESSION_BYTES = 200 * 2**20
"""A session file larger than this is refused (a 200 MB sequence is already far beyond what the analysis can handle)."""


def load_session(filename):
    """Load a session dict from a JSON file.

    Args:
        filename: Full path to the JSON file (including the .json extension).

    Returns:
        The parsed session dict on success, or None (with a message printed) if the file is not found,
        is too large, cannot be read or decoded as JSON, or is JSON that is not a session (not an object).
    """
    try:
        if os.path.getsize(filename) > MAX_SESSION_BYTES:
            print(f"Error: '{filename}' is too large to be a session file.")
            return None
        with open(filename, 'r', encoding='utf-8') as file:
            data = json.load(file)
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return None
    except (json.JSONDecodeError, UnicodeDecodeError, RecursionError):
        print(f"Error: Failed to decode JSON from the file '{filename}' (invalid format).")
        return None
    except OSError as exc:
        print(f"Error: Could not read '{filename}': {exc}")
        return None
    if not isinstance(data, dict):
        print(f"Error: '{filename}' is not a session file (expected a JSON object, got {type(data).__name__}).")
        return None
    return data


def session_problem(data) -> Optional[str]:
    """What is wrong with a loaded session dict, or None if it is usable.

    Checks what the analysis relies on and a corrupt or hand-edited file would break: the text is
    a string; any alphabet is a valid description; and the text, motifs and core list contain only
    characters of that alphabet (and the segment boundary).  A session without an alphabet is nucleotide text,
    and then the lists must be lists of strings.
    """
    from rna_lexis.alphabet import Alphabet, BOUNDARY
    if not isinstance(data, dict):
        return "the session is not a JSON object"
    for key in ("txt", "file_path"):
        if not isinstance(data.get(key), str):
            return f"the session has no {key!r} (a string)"
    for key in ("corelist", "xmotifs"):
        v = data.get(key, [])
        if not isinstance(v, list) or not all(isinstance(x, str) for x in v):
            return f"{key!r} must be a list of strings"
    alphabet = data.get("alphabet")
    if alphabet is None:
        return None
    try:
        a = alphabet if isinstance(alphabet, Alphabet) else Alphabet.from_dict(alphabet)
    except ValueError as exc:
        return f"the session's alphabet is not valid: {exc}"
    allowed = set(a.chars) | {BOUNDARY}
    for what, text in (("text", data["txt"]), *((k, "".join(data.get(k, []))) for k in ("corelist", "xmotifs"))):
        bad = sorted(set(text) - allowed)
        if bad:
            return (f"the session's {what} has characters that are not in its alphabet "
                    f"{a.name!r} (first: {bad[:3]!r}); the file may be damaged or edited by hand")
    return None


def is_valid_session(filepath: str) -> Optional[SessionData]:
    """Return a SessionData if filepath is a valid RNA_explore session JSON, else None."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            d = json.load(f)
        if not isinstance(d, dict) or not _SESSION_REQUIRED_KEYS.issubset(d):
            return None
        return SessionData.from_dict(d)
    except Exception:
        return None


def _is_file_open_in_libreoffice(filepath):
    """Return True if *filepath* is currently locked by LibreOffice."""
    dirpath = os.path.dirname(os.path.abspath(filepath))
    basename = os.path.basename(filepath)
    return os.path.exists(os.path.join(dirpath, f'.~lock.{basename}#'))


def _try_raise_libreoffice_window(filepath):
    """Attempt to raise the LibreOffice window for *filepath* to the foreground.
    Returns True if a window-manager command succeeded."""
    stem = os.path.splitext(os.path.basename(filepath))[0]
    for cmd in (
        ['wmctrl', '-a', stem],
        ['xdotool', 'search', '--name', stem, 'windowactivate', '--sync'],
    ):
        try:
            if subprocess.run(cmd, capture_output=True, timeout=3).returncode == 0:
                return True
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
    return False


def _summary_inputs_hash(xm, cores, txt, mutr, M, alphabet=None):
    """Stable hash of the parameters that determine the content of the summary CSV."""
    items = [sorted(xm), sorted(set(cores) - set(xm)), txt, mutr, M]
    if alphabet is not None:
        items.append(alphabet.to_dict())     # a different alphabet changes the file
    payload = json.dumps(items, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()[:20]


def _read_summary_hash(csv_path):
    """Return the stored input hash for *csv_path*, or None if absent/unreadable."""
    try:
        with open(csv_path + '.chk', 'r', encoding='utf-8') as fh:
            return fh.read().strip() or None
    except Exception:
        return None


def _write_summary_hash(csv_path, h):
    """Persist *h* alongside *csv_path* for future change-detection."""
    try:
        with open(csv_path + '.chk', 'w', encoding='utf-8') as fh:
            fh.write(h)
    except Exception:
        pass


def open_file_with_default_software(filename):
    """Opens a file using the default application for the current OS.

    Best-effort: on failure (no default app registered, missing
    xdg-open/open binary, headless environment, etc.) this prints a
    message instead of raising, so a save-and-open action never crashes
    the caller when only the "open" half fails.
    """
    try:
        if platform.system() == 'Windows':
            os.startfile(filename) # For Windows
        elif platform.system() == 'Darwin':
            subprocess.run(['open', filename], check=True) # For macOS
        else:
            # On Linux, skip xdg-open if LibreOffice already has the file open.
            if _is_file_open_in_libreoffice(filename):
                _try_raise_libreoffice_window(filename)
            else:
                subprocess.run(['xdg-open', filename], check=True)
    except Exception as e:
        print(f"Could not open file automatically: {filename}")
        print(str(e))


def _close_in_excel(filepath):
    """On Windows, close the workbook in Excel if it is currently open.
    Returns True if the workbook was found and closed, False otherwise.
    Silently does nothing (returns False) when win32com is unavailable."""
    try:
        import win32com.client
        xl = win32com.client.GetActiveObject("Excel.Application")
        target = os.path.abspath(filepath).lower()
        for wb in xl.Workbooks:
            if wb.FullName.lower() == target:
                wb.Close(SaveChanges=False)
                return True
    except Exception:
        pass
    return False


def split_fasta(text: str):
    """Split FASTA text into ``(headers, records)``, or return None if the text is not FASTA.

    Text is FASTA when its first line that is not blank and not a ``;`` comment starts with ``>``.  Lines before
    that, and any ``;`` comment line inside, are dropped; each record is the lines under one header, joined.
    (Testing only whether the whole text starts with ``>`` took a file that begins with a comment line to be
    plain text, and its header lines were then read as sequence.)
    """
    lines = text.splitlines()
    first = next((i for i, ln in enumerate(lines) if ln.strip() and not ln.lstrip().startswith(";")), None)
    if first is None or not lines[first].startswith(">"):
        return None
    headers, records = [], []
    for ln in lines[first:]:
        if ln.startswith(">"):
            headers.append(ln[1:].strip())
            records.append([])
        elif not ln.lstrip().startswith(";"):
            records[-1].append(ln)
    return headers, ["\n".join(r) for r in records]


def read_text(filename: str, keepspacing=False) -> str:
    """Read a text file and return its contents as a lower-case alphabetic string.

    Args:
        filename:    Path to the input file.
        keepspacing: When True, consecutive non-alphabetic characters are
                     replaced by a single space so word boundaries are
                     preserved.  When False (default), all non-alphabetic
                     characters (including newlines) are removed entirely.

    Returns:
        Lower-case string containing only [a-z] (and spaces when keepspacing
        is True).

    Raises:
        ValueError: If the file does not exist.
    """
    if not os.path.isfile(filename):
        raise ValueError(f"File {filename} not found.")
    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        txt = f.readlines()
    txt = ''.join(txt)
    if keepspacing:
        txt = txt.replace("\n", " ").lower()
        txt = re.sub('[^a-z]+', ' ', txt)
        txt = re.sub(' +', ' ', txt)
    else:
        txt = txt.replace("\n", "").lower()
        txt = re.sub('[^a-z]+', '', txt)
    return txt


def open_pdf(filepath):
    """Opens a PDF file in the default system viewer."""
    try:
        if sys.platform == "win32":
            os.startfile(filepath)
        elif sys.platform == "darwin": # macOS
            subprocess.run(["open", filepath], check=False)
        else: # Linux
            subprocess.run(["xdg-open", filepath], check=False)
    except Exception as e:
        print(f"Could not open output file automatically: {filepath}")
        print(str(e))


def _find_valid_sessions(directory: str) -> list:
    """Return sorted list of (display_name, filepath) for valid session JSON files in directory."""
    result = []
    for f in sorted(glob.glob(os.path.join(directory, '*.json'))):
        if is_valid_session(f) is not None:
            result.append((os.path.basename(f), f))
    return result


MAX_RESPONSE_BYTES = 100 * 2**20


def _read_capped(response) -> bytes:
    """The body of a response, refusing one larger than MAX_RESPONSE_BYTES instead of filling memory with it."""
    body = response.read(MAX_RESPONSE_BYTES + 1)
    if len(body) > MAX_RESPONSE_BYTES:
        raise ValueError(f"the response is larger than {MAX_RESPONSE_BYTES // 2**20} MB")
    return body


def _fetch_ensembl_json(url: str):
    """Fetch a JSON response from the Ensembl REST API.

    Args:
        url: Full Ensembl REST endpoint URL.

    Returns:
        Parsed JSON response as a Python dict or list.

    Raises:
        urllib.error.URLError: On network errors or non-200 HTTP responses.
        json.JSONDecodeError: If the response body is not valid JSON.
    """
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=45) as response:
        return json.loads(_read_capped(response).decode("utf-8"))


def fetch_enst_cdna(enst: str):
    """Fetch a cDNA sequence and display name for an Ensembl transcript ID.

    Tries the versioned ID first (e.g. ENST00000526704.3), then falls back
    to the unversioned base ID.  Both the sequence and lookup metadata are
    retrieved from the Ensembl REST API.

    Args:
        enst: Ensembl transcript identifier, with or without version suffix
              (e.g. ``'ENST00000526704'`` or ``'ENST00000526704.3'``).

    Returns:
        Tuple (transcript_id, display_name, sequence) where sequence is a
        lower-case cDNA string.

    Raises:
        ValueError: If enst does not start with ``'ENST'``.
        RuntimeError: If the transcript cannot be fetched from the API.
    """
    if not isinstance(enst, str):
        raise ValueError("Transcript ID must be a string")
    enst = enst.strip()
    if not enst.upper().startswith("ENST"):
        raise ValueError("Transcript ID must start with ENST")
    # Letters, digits and one version suffix only: the ID is placed into a URL, so anything else
    # (a slash, '?', '#', a space) would change which address is requested.
    if not re.fullmatch(r"ENST[0-9]+(\.[0-9]+)?", enst, flags=re.I):
        raise ValueError("Transcript ID must look like ENST00000526704 or ENST00000526704.3")

    candidates = [enst]
    if "." in enst:
        candidates.append(enst.split(".")[0])

    last_error = None
    for tid in candidates:
        seq_url = f"https://rest.ensembl.org/sequence/id/{tid}?type=cdna"
        lookup_url = f"https://rest.ensembl.org/lookup/id/{tid}"
        try:
            seq_data = _fetch_ensembl_json(seq_url)
            seq = seq_data.get("seq", "")
            if not seq:
                raise ValueError(f"No sequence returned for {tid}")
            lookup = _fetch_ensembl_json(lookup_url)
            label = lookup.get("display_name", tid)
            return tid, label, seq.lower()
        except Exception as exc:
            last_error = exc
    raise RuntimeError(f"Could not fetch transcript {enst}: {last_error}")


# hg38 chromosome sizes, used to page through the UCSC track API for
# fetch_encode_ccre() (its `encodeCcreCombined` track endpoint requires an
# explicit chrom:start-end window rather than accepting a bare accession).
_HG38_SIZES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}


def _fetch_ucsc_json(url: str):
    """Fetch a JSON response from the UCSC genome-browser REST API.

    Args:
        url: Full UCSC API endpoint URL.

    Returns:
        Parsed JSON response as a Python dict.

    Raises:
        urllib.error.URLError: On network errors or non-200 HTTP responses.
        json.JSONDecodeError: If the response body is not valid JSON.
    """
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as response:
        return json.loads(_read_capped(response).decode("utf-8"))


MAX_RANGE = 10_000_000


def _check_genome(genome):
    """A UCSC assembly name (hg38, mm10, ...) is letters and digits; it is placed into a URL."""
    if not isinstance(genome, str) or not re.fullmatch(r"[A-Za-z0-9]{2,20}", genome):
        raise ValueError(f"genome must be an assembly name such as hg38, got {genome!r}")


def fetch_genomic_range(chrom: str, start: int, end: int, genome: str = "hg38"):
    """Fetch a raw genomic DNA sequence by coordinates via the UCSC REST API.

    Args:
        chrom:  Chromosome name (e.g. ``'chr10'``).
        start:  0-based start position.
        end:    0-based, exclusive end position.
        genome: UCSC genome assembly name (default ``'hg38'``).

    Returns:
        Lower-case DNA sequence string for the requested range.

    Raises:
        ValueError: If the API returns no sequence for the given range.
    """
    _check_genome(genome)
    if not isinstance(chrom, str) or not re.fullmatch(r"chr[0-9A-Za-z_]{1,40}", chrom):
        raise ValueError(f"chromosome must look like chr10 or chrX, got {chrom!r}")
    if isinstance(start, bool) or isinstance(end, bool) or not isinstance(start, int) \
            or not isinstance(end, int) or start < 0 or end <= start:
        raise ValueError(f"start and end must be whole numbers with 0 <= start < end, got {start!r}, {end!r}")
    if end - start > MAX_RANGE:
        raise ValueError(f"the range is {end - start:,} bases; the most that can be fetched is {MAX_RANGE:,}")
    seq_url = (f"https://api.genome.ucsc.edu/getData/sequence"
               f"?genome={genome};chrom={chrom};start={start};end={end}")
    seq = _fetch_ucsc_json(seq_url).get("dna", "")
    if not seq:
        raise ValueError(f"No sequence for {chrom}:{start}-{end} ({genome})")
    return seq.lower()


def fetch_encode_ccre(accession: str, genome: str = "hg38", hint_chrom: str = None):
    """Look up an ENCODE cCRE accession and fetch its forward-strand sequence.

    Searches the UCSC ``encodeCcreCombined`` track in 10 Mb windows across
    every chromosome for an item whose ``name`` matches `accession`, then
    fetches that region's DNA via `fetch_genomic_range()`.

    Args:
        accession:  ENCODE cCRE accession (e.g. ``'EH38E1482203'``).
        genome:     UCSC genome assembly name (default ``'hg38'``).
        hint_chrom: If given, this chromosome is searched first (speeds up
                    lookup when the region is already known, e.g. from a
                    nearby gene).

    Returns:
        Tuple (annotation_string, sequence) where sequence is a lower-case
        DNA string.

    Raises:
        ValueError: If the accession is not found, or has no sequence.
    """
    _check_genome(genome)
    if not isinstance(accession, str) or not re.fullmatch(r"[A-Za-z0-9]{4,30}", accession.strip()):
        raise ValueError(f"an ENCODE accession is letters and digits, like EH38E1482203, got {accession!r}")
    accession = accession.strip()
    chroms = list(_HG38_SIZES)
    if hint_chrom and hint_chrom in _HG38_SIZES:
        chroms = [hint_chrom] + [c for c in chroms if c != hint_chrom]
    window = 10_000_000
    item = None
    for chrom in chroms:
        for s0 in range(0, _HG38_SIZES[chrom], window):
            url = (f"https://api.genome.ucsc.edu/getData/track"
                   f"?genome={genome};track=encodeCcreCombined"
                   f";chrom={chrom};start={s0};end={min(s0 + window, _HG38_SIZES[chrom])}")
            try:
                for it in _fetch_ucsc_json(url).get("encodeCcreCombined", []):
                    if it.get("name") == accession:
                        item = it
                        break
            except Exception:
                pass
            if item:
                break
        if item:
            break
    if not item:
        raise ValueError(f"Accession {accession!r} not found in {genome}")

    chrom, start, end = item["chrom"], item["chromStart"], item["chromEnd"]
    ccre = item.get("ccre", "")
    annotation = (f">{accession} {item.get('description', accession)} ({ccre}) "
                  f"{chrom}:{start}-{end} strand=. genome={genome}")
    seq = fetch_genomic_range(chrom, start, end, genome=genome)
    return annotation, seq


def init_summary(fn, xm, cores, txt, mutr=1/6, M = 4, force=False, alphabet=None):
    '''Input: xmotifs, cores, and the gene sequence. Output: a
    report summarizing the properties of the sequences.

    With an :class:`rna_lexis.alphabet.Alphabet`, ``xm``, ``cores`` and ``txt``
    are internal strings: the statistical columns are computed over the declared
    alphabet (default Markov order and multiple-testing family for it), and the
    ``seq`` column shows each motif as the alphabet's labels.'''
    out_csv = init_summary_path(fn)
    current_hash = _summary_inputs_hash(xm, cores, txt, mutr, M, alphabet)
    if os.path.isfile(out_csv) and not force:
        try:
            with open(out_csv, 'r', encoding='utf-8', errors='ignore') as handle:
                header = handle.readline().strip().split(',')
        except Exception:
            header = []
        required_test_columns = {'expected_markov', 'enrichment_markov', 'p_markov', 'q_markov'}
        if required_test_columns.issubset(set(header)):
            stored_hash = _read_summary_hash(out_csv)
            if stored_hash is None or stored_hash == current_hash:
                try:
                    open_file_with_default_software(out_csv)
                except Exception:
                    pass
                return out_csv
            print(f"Input parameters changed; regenerating: {out_csv}")
        else:
            print(f"Existing summary lacks statistical columns; regenerating: {out_csv}")


    xm_cores = xm + list(set(cores) - set(xm))
    txt_lower = txt.lower()

    # Pre-compute core indices for O(1) lookup instead of linear search
    core_idx_map = {core: idx for idx, core in enumerate(cores)}

    # Collect records in list for batch DataFrame creation (O(1) instead of O(n²))
    records = []

    for i, seq in enumerate(xm_cores):
        is_xm   = 1 if i < len(xm) else 0
        is_core = 1 if seq in core_idx_map else 0
        idx     = i if is_xm else core_idx_map.get(seq, -1)

        L = len(seq)
        pos, approx, maxmut = find_with_mutations(seq, txt, mutr=mutr, M=M)
        N = len(pos)
        n_approx = len(approx)
        N_total = N + n_approx          # exact + approximate occurrences
        nmut = sum(d for _, _, d in approx)

        # Observed mutation rate per site across all occurrences
        pmut = nmut / (L * N_total) if N_total > 0 else 0.0

        # Binomial p-value: P(X <= nmut) where X ~ Bin(L*N_total, mutr).
        # Small p_stable means the sequence is significantly more conserved
        # than the null mutation rate (mutr), i.e., it is stable.
        p_stable = _binom.cdf(nmut, L * N_total, mutr) if N_total > 0 else 1.0

        records.append({
            'seq': seq,
            'xm': is_xm,
            'core': is_core,
            'len': L,
            'count': N,
            'n_approx': n_approx,
            'cover': L * N,
            'numt': nmut,
            'pmut': pmut,
            'p_stable': p_stable,
            'maxmut': maxmut,
            'idx': idx,
            'pos': pos
        })

    # Create DataFrame once from all records (not row-by-row append)
    df = pd.DataFrame(records)

    # Statistical layer: Markov enrichment/FDR columns, for nucleotide text or
    # for text read over an alphabet.
    if alphabet is not None or set(txt_lower).issubset(set('acgtu')):
        try:
            score_rows = score_exact_motifs(
                txt if alphabet is not None else txt_lower,
                xm_cores,
                xmotifs=xm,
                markov_order=1 if alphabet is None else None,
                enrichment_threshold=10.0,
                min_xmotif_type_support=1,
                alphabet=alphabet,
            )
            score_map = {row['motif'].lower(): row for row in score_rows}
            mapped_columns = {
                'nonoverlap_count': 'nonoverlap_count',
                'coverage_bp': 'coverage_bp',
                'area_score': 'area_score',
                'xmotif_type_support': 'xmotif_type_support',
                'inside_xmotif_count': 'inside_xmotif_count',
                'outside_xmotif_count': 'outside_xmotif_count',
                'expected_markov': 'expected_markov',
                'enrichment_markov': 'enrichment_markov',
                'p_markov': 'p_markov',
                'q_markov': 'q_markov',
                'statistically_supported': 'statistically_supported',
                'core_class': 'core_class',
                'rank_statistical': 'rank_statistical',
                'rank_coverage': 'rank_coverage',
            }
            for out_col, score_col in mapped_columns.items():
                df[out_col] = df['seq'].map(
                    lambda s: score_map.get(str(s).lower(), {}).get(score_col)
                )
        except Exception as exc:
            df['statistical_support_error'] = str(exc)

    if alphabet is not None:
        df['seq'] = df['seq'].map(lambda s: alphabet.show(s))
    if 'q_markov' in df.columns:
        df = df.sort_values(
            by=['statistically_supported', 'q_markov', 'p_stable'],
            ascending=[False, True, True],
            na_position='last',
        )
    else:
        df = df.sort_values(by='p_stable', ascending=True)
    for col in df.select_dtypes(include="object").columns:      # a symbol label must not start a spreadsheet formula
        df[col] = df[col].map(csv_safe)
    out_csv = init_summary_path(fn)
    tried = []
    for candidate in [
        out_csv,
        os.path.join(os.getcwd(), os.path.basename(out_csv)),
        os.path.join(os.path.expanduser("~"), os.path.basename(out_csv)),
        os.path.join(tempfile.gettempdir(), os.path.basename(out_csv)),
    ]:
        try:
            df.to_csv(candidate, index=False)
            _write_summary_hash(candidate, current_hash)
            print(f"Wrote summary to: {candidate}")
            open_file_with_default_software(candidate)
            return candidate
        except PermissionError:
            if platform.system() == 'Windows' and _close_in_excel(candidate):
                try:
                    df.to_csv(candidate, index=False)
                    _write_summary_hash(candidate, current_hash)
                    print(f"Wrote summary to: {candidate}")
                    open_file_with_default_software(candidate)
                    return candidate
                except Exception:
                    pass
            tried.append(candidate)
            continue
    raise PermissionError(f"Unable to write summary CSV. Tried: {tried}")

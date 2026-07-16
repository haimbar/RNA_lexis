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
from rna_lexis.statistical import score_exact_motifs

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


def load_session(filename):
    """Load a session dict from a JSON file.

    Args:
        filename: Full path to the JSON file (including the .json extension).

    Returns:
        The parsed session dict on success, or None if the file is not found
        or cannot be decoded as valid JSON.
    """
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            return(json.load(file))
    except FileNotFoundError:
        print(f"Error: The file '{filename}.json' was not found.")
    except json.JSONDecodeError:
        print(f"Error: Failed to decode JSON from the file (invalid format).")


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


def _summary_inputs_hash(xm, cores, txt, mutr, M):
    """Stable hash of the parameters that determine the content of the summary CSV."""
    payload = json.dumps(
        [sorted(xm), sorted(set(cores) - set(xm)), txt, mutr, M],
        sort_keys=True, default=str,
    )
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
        return json.loads(response.read().decode("utf-8"))


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
    enst = enst.strip()
    if not enst.upper().startswith("ENST"):
        raise ValueError("Transcript ID must start with ENST")

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


def init_summary(fn, xm, cores, txt, mutr=1/6, M = 4, force=False):
    '''Input: xmotifs, cores, and the gene sequence. Output: a
    report summarizing the properties of the sequences.'''
    out_csv = init_summary_path(fn)
    current_hash = _summary_inputs_hash(xm, cores, txt, mutr, M)
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

    # Statistical layer: Markov enrichment/FDR columns for RNA/DNA sequences only.
    if set(txt_lower).issubset(set('acgtu')):
        try:
            score_rows = score_exact_motifs(
                txt_lower,
                xm_cores,
                xmotifs=xm,
                markov_order=1,
                enrichment_threshold=10.0,
                min_xmotif_type_support=1,
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

    if 'q_markov' in df.columns:
        df = df.sort_values(
            by=['statistically_supported', 'q_markov', 'p_stable'],
            ascending=[False, True, True],
            na_position='last',
        )
    else:
        df = df.sort_values(by='p_stable', ascending=True)
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

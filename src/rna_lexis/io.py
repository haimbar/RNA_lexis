"""File I/O, session management, and network fetch utilities."""

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
from dataclasses import dataclass, asdict
from typing import List, Optional

import pandas as pd
from platformdirs import user_config_dir

from rna_lexis.algorithms import find_with_mutations


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
    path = _prefs_path()
    if not os.path.isfile(path):
        return defaults
    try:
        with open(path, 'r') as f:
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
        with open(_prefs_path(), 'w') as f:
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
    with open(f"{filename}.json", "w") as json_file:
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
        with open(filename, 'r') as file:
            return(json.load(file))
    except FileNotFoundError:
        print(f"Error: The file '{filename}.json' was not found.")
    except json.JSONDecodeError:
        print(f"Error: Failed to decode JSON from the file (invalid format).")


def is_valid_session(filepath: str) -> Optional[SessionData]:
    """Return a SessionData if filepath is a valid RNA_explore session JSON, else None."""
    try:
        with open(filepath, 'r') as f:
            d = json.load(f)
        if not isinstance(d, dict) or not _SESSION_REQUIRED_KEYS.issubset(d):
            return None
        return SessionData.from_dict(d)
    except Exception:
        return None


def open_file_with_default_software(filename):
    """Opens a file using the default application for the current OS."""
    if platform.system() == 'Windows':
        os.startfile(filename) # For Windows
    elif platform.system() == 'Darwin':
        subprocess.run(['open', filename], check=True) # For macOS
    else:
        # Assume Linux or other POSIX-like system
        # xdg-open is a common standard on many Linux distributions
        subprocess.run(['xdg-open', filename], check=True)


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


def init_summary(fn, xm, cores, txt, mutr=1/6, M = 4):
    '''Input: xmotifs, cores, and the gene sequence. Output: a
    report summarizing the properties of the sequences.'''
    out_csv = f'{fn}_init.csv'
    if os.path.isfile(out_csv):
        try:
            open_file_with_default_software(out_csv)
        except Exception:
            pass
        return out_csv

    xm_cores = xm + list(set(cores) - set(xm))
    txt_lower = txt.lower()

    # Pre-compute core indices for O(1) lookup instead of linear search
    core_idx_map = {core: idx for idx, core in enumerate(cores)}

    # Pre-compute which cores are in each xm (avoid nested loop)
    xm_core_containment = {i: any(core in xm[i] for core in cores)
                           for i in range(len(xm))}

    # Collect records in list for batch DataFrame creation (O(1) instead of O(n²))
    records = []

    for i, seq in enumerate(xm_cores):
        # Determine sequence type and index
        if i < len(xm):
            idx = i
            seqtype = 'xm'
            inccore = 1 if xm_core_containment[i] else 0
        else:
            seqtype = 'core'
            idx = core_idx_map.get(seq, -1)
            inccore = 0

        L = len(seq)
        pos, approx, maxmut = find_with_mutations(seq, txt, mutr=mutr, M=M)
        N = len(pos)
        # mutations = all matches within maxmut (exact + approx); nmut = total distance
        mutations = [txt_lower[p:p + L] for p in pos] + [m for _, m, _ in approx]
        nmut = sum(d for _, _, d in approx)

        # Compute pmut with safe denominator
        pmut = (nmut + 0.01) / (L * (N + len(mutations)) + 1)

        records.append({
            'seq': seq,
            'type': seqtype,
            'hascore': inccore,
            'len': L,
            'count': N,
            'cover': L * N,
            'numt': nmut,
            'pmut': pmut,
            'maxmut': maxmut,
            'idx': idx,
            'pos': pos
        })

    # Create DataFrame once from all records (not row-by-row append)
    df = pd.DataFrame(records)
    df = df.sort_values(by='pmut', ascending=True)
    out_csv = f'{fn}_init.csv'
    tried = []
    for candidate in [
        out_csv,
        os.path.join(os.getcwd(), os.path.basename(out_csv)),
        os.path.join(os.path.expanduser("~"), os.path.basename(out_csv)),
        os.path.join(tempfile.gettempdir(), os.path.basename(out_csv)),
    ]:
        try:
            df.to_csv(candidate, index=False)
            print(f"Wrote summary to: {candidate}")
            open_file_with_default_software(candidate)
            return candidate
        except PermissionError:
            if platform.system() == 'Windows' and _close_in_excel(candidate):
                try:
                    df.to_csv(candidate, index=False)
                    print(f"Wrote summary to: {candidate}")
                    open_file_with_default_software(candidate)
                    return candidate
                except Exception:
                    pass
            tried.append(candidate)
            continue
    raise PermissionError(f"Unable to write summary CSV. Tried: {tried}")

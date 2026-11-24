"""Convert between units and nucleotide text, and read/write codebooks.

The analysis works on units, so this module is needed only when the
nucleotide-level view is wanted: to compare a unit-level run with an
encoded-nucleotide run, to export a unit sequence as codon text, or to decode
a motif found in encoded text back to its units.

A *codebook* maps each unit (``"a"``, a syllable name, a codon) to a
fixed-width nucleotide code.  The assignment is stored in the codebook, so an
encoding can be reproduced without the seed that made it.

This is the general, package-side counterpart of GeneDB's ``eng2rna`` codec.
It works for any symbol set and any code width, not only a-z and 3 bases.
"""

from __future__ import annotations

import csv
import itertools
import random
import re
from typing import Dict, Mapping, Optional, Sequence, Union

from rna_lexis.alphabet import BOUNDARY, Alphabet
from rna_lexis.tokenizers import Tokenized, tokenize


# ----------------------------------------------------------------------
# codebooks
# ----------------------------------------------------------------------
def make_codebook(symbols: Sequence[str], width: int = 3, bases: str = "tcga",
                  seed: Optional[int] = None) -> Dict[str, str]:
    """Assign each symbol a distinct nucleotide code using a seeded permutation.

    All ``len(bases) ** width`` codes are shuffled and the first ``len(symbols)``
    are handed out in order, so the assignment is reproducible from ``seed``
    and not tied to the order of the symbols in any alphabet.

    Args:
        symbols: Labels to encode.
        width:   Code length in bases.
        bases:   The base characters.
        seed:    Seed for the shuffle; ``None`` uses system entropy (keep the
                 returned codebook, since the seed is then not recoverable).

    Raises:
        ValueError: if there are more symbols than codes.
    """
    symbols = list(symbols)
    if len(set(symbols)) != len(symbols):
        raise ValueError("symbols must be distinct")
    codes = ["".join(p) for p in itertools.product(bases, repeat=width)]
    if len(symbols) > len(codes):
        raise ValueError(
            f"{len(symbols)} symbols need more than the {len(codes)} codes "
            f"available from {len(bases)} bases at width {width}; raise the width")
    random.Random(seed).shuffle(codes)
    return dict(zip(symbols, codes))


def read_codebook_csv(path: str, symbol_col: str = "symbol", code_col: str = "code",
                      where: Optional[Mapping[str, str]] = None) -> Dict[str, str]:
    """Read ``{symbol: code}`` from a CSV file.

    Codes are lower-cased and ``u`` is read as ``t``.  ``where`` keeps only the
    rows whose columns match, for a file holding several codebooks, e.g.
    ``where={"bird": "Atama"}`` for GeneDB's ``birdsong_codebook.csv`` (columns
    ``bird, syllable, codon, count``) with ``symbol_col="syllable"`` and
    ``code_col="codon"``.

    Raises:
        ValueError: if a column is missing, or a symbol appears with two codes.
    """
    book: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        needed = [symbol_col, code_col] + list(where or {})
        missing = [c for c in needed if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: no column(s) {missing}; found {reader.fieldnames}")
        for row in reader:
            if where and any(row[k] != v for k, v in where.items()):
                continue
            sym, code = row[symbol_col], row[code_col].strip().lower().replace("u", "t")
            if sym in book and book[sym] != code:
                raise ValueError(
                    f"{path}: symbol {sym!r} has two codes ({book[sym]!r}, {code!r}); "
                    "use `where` to select one codebook")
            book[sym] = code
    return book


def write_codebook_csv(path: str, codebook: Mapping[str, str], symbol_col: str = "symbol",
                       code_col: str = "code") -> str:
    """Write ``{symbol: code}`` as a two-column CSV file and return the path."""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([symbol_col, code_col])
        for sym, code in codebook.items():
            writer.writerow([sym, code])
    return path


# ----------------------------------------------------------------------
# units <-> nucleotides
# ----------------------------------------------------------------------
def _as_internal(units: Union[str, Tokenized]) -> str:
    return units.text if isinstance(units, Tokenized) else units


def code_width(alphabet: Alphabet) -> int:
    """Length of the nucleotide codes of an alphabet."""
    if not alphabet.codes:
        raise ValueError(f"alphabet {alphabet.name!r} has no nucleotide codes")
    return len(next(iter(alphabet.codes.values())))


def encode_nt(units: Union[str, Tokenized], alphabet: Alphabet, extras: str = "error") -> str:
    """Internal unit string (or a :class:`Tokenized`) -> nucleotide text.

    Segment boundaries stay as ``"_"``.

    Args:
        units:    Internal string, or a Tokenized result.
        alphabet: The alphabet the units belong to; it must have codes.
        extras:   What to do with an extra symbol (space, digit), which has no
                  nucleotide code: ``"error"`` (default), ``"keep"`` the
                  character as it is, or ``"drop"`` it.
    """
    if extras not in ("error", "keep", "drop"):
        raise ValueError("extras must be 'error', 'keep' or 'drop'")
    code_width(alphabet)    # raises if there are no codes
    out = []
    for c in _as_internal(units):
        if c == BOUNDARY:
            out.append(BOUNDARY)
            continue
        label = alphabet.label(c)
        if label in alphabet.codes:
            out.append(alphabet.codes[label])
        elif extras == "keep":
            out.append(label)
        elif extras == "error":
            raise ValueError(f"symbol {label!r} has no nucleotide code; use extras='keep' or 'drop'")
    return "".join(out)


def decode_nt(nt: str, alphabet: Alphabet, frame: int = 0,
              on_unknown: str = "skip") -> Tokenized:
    """Nucleotide text -> units, the inverse of :func:`encode_nt`.

    Text is read in groups of the code width, starting ``frame`` bases in, and
    each group is looked up in the alphabet's codes.  Runs of ``"_"`` split the
    text into segments and restart the reading (one boundary character joins
    them again in the result).  Extra symbols kept by ``encode_nt(...,
    extras="keep")`` are read back as extras.

    To decode a motif that was found at position ``pos`` of the nucleotide
    text, use ``frame=frame_for_position(pos, code_width(alphabet))``.

    Returns:
        A :class:`Tokenized` for ``alphabet``.  Unknown groups and partial
        groups are handled as in :func:`rna_lexis.tokenizers.tokenize`.
    """
    width = code_width(alphabet)
    code_to_label = {code: label for label, code in alphabet.codes.items()}
    code_alpha = Alphabet(name=f"{alphabet.name} (codes)", symbols=tuple(code_to_label),
                          extras=alphabet.extras, mode="fixed", width=width,
                          aliases={"u": "t"})
    segments = [s for s in re.split(f"{re.escape(BOUNDARY)}+", nt)]
    tok = tokenize(segments, code_alpha, frame=frame, on_unknown=on_unknown)

    table = {code_alpha.char(code): alphabet.char(label) for code, label in code_to_label.items()}
    table.update({code_alpha.char(x): alphabet.char(x) for x in alphabet.extras})
    table[BOUNDARY] = BOUNDARY
    tok.text = "".join(table[c] for c in tok.text)
    tok.alphabet = alphabet
    return tok


def frame_for_position(pos: int, width: int = 3) -> int:
    """Characters to skip so that reading starting at ``pos`` lines up with a unit.

    ``pos`` is a 0-based position in nucleotide text whose units start at
    position 0.  A motif found at ``pos`` is read with this frame.
    """
    return (-pos) % width


# Open reading frames turn up by chance in any RNA, long non-coding ones included (the human lncRNA NORAD, 5,400
# nucleotides, has one of 108 codons).  Under MIN_CODING_CODONS the text is not offered as a coding region at all;
# from there to LIKELY_CODING_CODONS it is offered with a warning; from LIKELY_CODING_CODONS on it is the default.
MIN_CODING_CODONS = 100
LIKELY_CODING_CODONS = 300
_STOP_CODONS = ("taa", "tag", "tga")


def find_orf(nt: str):
    """The longest open reading frame of a transcript: ``ATG`` to the first stop codon, forward strand only.

    Looks in the three frames (``u`` is read as ``t``; case does not matter) and returns the longest stretch
    that starts at an ``atg`` and ends at a stop, with the stop included, as a dict
    ``{"start", "end", "frame", "codons"}``: ``nt[start:end]`` is the region, ``frame`` is ``start % 3``, and
    ``codons`` counts the codons before the stop.  A stretch with no stop before the text ends is not an
    open reading frame here.  None if there is none (or the text is empty).

    A transcript from a protein-coding gene has one that is hundreds or thousands of codons long
    (:data:`LIKELY_CODING_CODONS`); a long non-coding RNA does not.  This is not a gene finder: it does not look
    at the reverse strand, at introns, or at alternative start codons, and the longest ORF is not always the
    annotated one.
    """
    text = re.sub(r"[^a-z]", "", str(nt).lower().replace("u", "t"))
    best = None
    for frame in range(3):
        start = None
        for i in range(frame, len(text) - 2, 3):
            codon = text[i:i + 3]
            if start is None:
                if codon == "atg":
                    start = i
            elif codon in _STOP_CODONS:
                codons = (i - start) // 3
                if best is None or codons > best["codons"]:
                    best = {"start": start, "end": i + 3, "frame": start % 3, "codons": codons}
                start = None
    return best

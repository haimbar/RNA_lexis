"""Read text into an alphabet: text -> internal one-character-per-unit string.

:func:`tokenize` is the entry point.  It lower-cases the input (unless the
alphabet is case sensitive), applies the alphabet's aliases (``u`` -> ``t``),
splits the text into units, and returns a :class:`Tokenized` result that keeps
everything needed to map results back to the original text: where each unit
came from, what was skipped or dropped, and where a removal made two units
adjacent that were not adjacent in the source.

Three ways of reading, set by ``mode`` (the alphabet's own default is used
when none is given):

``"chars"``  every character is one unit (nucleotides, Latin letters).
``"fixed"``  ``width`` characters make one unit (codons).  A ``frame`` of 0, 1
             or 2 characters is dropped from the start of each segment.  Reading
             restarts after each extra symbol, and a stretch whose length is not
             a multiple of ``width`` loses its partial unit with a warning.
             A unit that is not a declared symbol (an ambiguity code such as
             ``n`` inside a codon) is an unknown *unit*: the whole group is
             handled at once, so the frame is never shifted by a bad character.
``"words"``  whitespace-separated labels (syllable names).

Whitespace is formatting and is removed unless the alphabet declares it as an
extra symbol.

Unknown symbols
---------------
``on_unknown="skip"`` (the default) removes them and issues an
:class:`UnknownSymbolWarning` saying how many, which, and where.
``on_unknown="error"`` raises :class:`UnknownSymbolError` at the first one.
Skipping makes the units on either side adjacent, so a motif found across a
removal may not exist in the original text; ``Tokenized.removal_points`` and
:meth:`Tokenized.spans_removal` let the caller flag such occurrences.
"""

from __future__ import annotations

import re
import warnings
from collections import Counter
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Union

from rna_lexis.alphabet import BOUNDARY, Alphabet


class UnknownSymbolError(ValueError):
    """A symbol that is not in the alphabet was found and on_unknown='error'."""


class UnknownSymbolWarning(UserWarning):
    """Unknown symbols were skipped."""


class PartialUnitWarning(UserWarning):
    """A stretch of text did not divide into whole units and a partial unit was dropped."""


@dataclass
class Skipped:
    """A symbol (or fixed-width group) that is not in the alphabet and was removed."""
    segment: int
    offset: int      # index in the original segment string
    text: str


@dataclass
class Dropped:
    """Characters dropped without being unknown: the reading frame or a partial unit."""
    segment: int
    offset: int
    text: str
    reason: str      # "frame" or "partial"


@dataclass
class Segment:
    """One input segment: where its units sit in :attr:`Tokenized.text`."""
    index: int
    start: int
    end: int
    source_length: int


@dataclass
class Tokenized:
    """Result of :func:`tokenize`.

    Attributes:
        alphabet:       The alphabet the text was read into.
        text:           Internal string, one character per unit, segments
                        joined by the boundary character ``"_"``.
        segments:       Where each segment sits in ``text``.
        offsets:        For each character of ``text``, the index in the
                        original segment string where that unit started;
                        -1 for a boundary character.
        skipped:        Unknown symbols that were removed.
        dropped:        Characters dropped for the reading frame or as a
                        partial unit.
        removal_points: Positions ``p`` in ``text`` where a removal made
                        ``text[p - 1]`` and ``text[p]`` adjacent although
                        they were not adjacent in the source.
    """
    alphabet: Alphabet
    text: str
    segments: List[Segment]
    offsets: List[int]
    skipped: List[Skipped] = field(default_factory=list)
    dropped: List[Dropped] = field(default_factory=list)
    removal_points: List[int] = field(default_factory=list)

    @property
    def n_units(self) -> int:
        """Number of units, boundaries not counted."""
        return sum(1 for c in self.text if c != BOUNDARY)

    def spans_removal(self, start: int, length: int) -> bool:
        """True if the window ``text[start:start+length]`` contains a removal point.

        A motif occurrence for which this is true joins two units that were
        separated in the source.
        """
        return any(start < p < start + length for p in self.removal_points)

    def decode(self, start: int = 0, stop: Optional[int] = None, sep: str = "") -> str:
        """The labels of ``text[start:stop]`` joined by ``sep``."""
        return self.alphabet.decode(self.text[start:stop], sep=sep)

    def project(self, projection) -> "Tokenized":
        """The same text over the projection's target alphabet (e.g. codons -> amino acids).

        Positions, segments, offsets and the record of what was skipped or dropped
        are unchanged: unit ``i`` of the result is unit ``i`` of this text.
        """
        if projection.source.to_dict() != self.alphabet.to_dict():
            raise ValueError(f"the projection is from {projection.source.name!r}, "
                             f"not from {self.alphabet.name!r}")
        return Tokenized(alphabet=projection.target, text=projection.apply(self.text),
                         segments=list(self.segments), offsets=list(self.offsets),
                         skipped=list(self.skipped), dropped=list(self.dropped),
                         removal_points=list(self.removal_points))


# ----------------------------------------------------------------------
# one segment
# ----------------------------------------------------------------------
class _SegmentReader:
    """Collects the units of one segment."""

    def __init__(self, alphabet: Alphabet, index: int, on_unknown: str):
        self.a = alphabet
        self.index = index
        self.on_unknown = on_unknown
        self.units: List[str] = []
        self.offsets: List[int] = []
        self.skipped: List[Skipped] = []
        self.dropped: List[Dropped] = []
        self.removals: List[int] = []

    # -- helpers --------------------------------------------------------
    def fold(self, s: str) -> str:
        s = s if self.a.case_sensitive else s.lower()
        return self.a.aliases.get(s, s)

    def emit(self, label: str, offset: int):
        self.units.append(self.a.char(label))
        self.offsets.append(offset)

    def unknown(self, offset: int, text: str):
        if self.on_unknown == "error":
            raise UnknownSymbolError(
                f"{text!r} at offset {offset} of segment {self.index} is not a symbol "
                f"of alphabet {self.a.name!r}")
        self.skipped.append(Skipped(self.index, offset, text))
        self.removals.append(len(self.units))

    # -- modes ----------------------------------------------------------
    def read_chars(self, seg: str):
        for i, ch in enumerate(seg):
            f = self.fold(ch)
            if f in self.a:
                self.emit(f, i)
            elif ch.isspace():
                continue                      # formatting
            else:
                self.unknown(i, ch)

    def read_words(self, seg: str):
        for m in re.finditer(r"\S+", seg):
            f = self.fold(m.group())
            if f in self.a:
                self.emit(f, m.start())
            else:
                self.unknown(m.start(), m.group())

    def read_fixed(self, seg: str, width: int, frame: int):
        extras = set(self.a.extras)
        state = {"run": [], "first": True}

        def flush():
            run = state["run"]
            if not run:
                return
            state["run"] = []
            if state["first"]:
                state["first"] = False
                if frame:
                    head, run = run[:frame], run[frame:]
                    if head:
                        self.dropped.append(Dropped(self.index, head[0][0],
                                                    "".join(c[2] for c in head), "frame"))
            whole = len(run) - len(run) % width
            for j in range(0, whole, width):
                chunk = run[j:j + width]
                label = "".join(c[1] for c in chunk)
                if label in self.a and label not in extras:
                    self.emit(label, chunk[0][0])
                else:
                    self.unknown(chunk[0][0], "".join(c[2] for c in chunk))
            if whole < len(run):
                tail = run[whole:]
                self.dropped.append(Dropped(self.index, tail[0][0],
                                            "".join(c[2] for c in tail), "partial"))
                self.removals.append(len(self.units))

        for i, ch in enumerate(seg):
            f = self.fold(ch)
            if f in extras:
                flush()
                self.emit(f, i)
            elif ch.isspace():
                continue                      # formatting
            else:
                state["run"].append((i, f, ch))
        flush()


# ----------------------------------------------------------------------
# entry point
# ----------------------------------------------------------------------
def tokenize(text: Union[str, Sequence[str]], alphabet: Alphabet, mode: Optional[str] = None,
             width: Optional[int] = None, frame: int = 0, on_unknown: str = "skip",
             gap: int = 1) -> Tokenized:
    """Read ``text`` into ``alphabet``.

    Args:
        text:       One string, or a sequence of strings (segments such as
                    song bouts or sentences).  Motifs must not span two
                    segments, so they are kept apart by boundary characters.
        alphabet:   The declared alphabet.
        mode:       ``"chars"``, ``"fixed"`` or ``"words"``; default from the
                    alphabet.
        width:      Characters per unit in ``"fixed"`` mode; default from the
                    alphabet.
        frame:      Characters to drop from the start of each segment in
                    ``"fixed"`` mode (the reading frame), 0 to ``width - 1``.
        on_unknown: ``"skip"`` (default, with a warning) or ``"error"``.
        gap:        Number of boundary characters between segments.

    Returns:
        A :class:`Tokenized` result.
    """
    if on_unknown not in ("skip", "error"):
        raise ValueError("on_unknown must be 'skip' or 'error'")
    mode = mode or alphabet.mode
    width = width or alphabet.width
    if gap < 1:
        raise ValueError("gap must be at least 1")
    if mode == "fixed" and not 0 <= frame < width:
        raise ValueError(f"frame must be between 0 and {width - 1}, got {frame}")
    if mode != "fixed" and frame:
        raise ValueError("frame applies only to 'fixed' mode")

    segments_in = [text] if isinstance(text, str) else list(text)
    out_text: List[str] = []
    offsets: List[int] = []
    segments: List[Segment] = []
    skipped: List[Skipped] = []
    dropped: List[Dropped] = []
    removal_points: List[int] = []

    pos = 0
    for k, seg in enumerate(segments_in):
        if k:
            out_text.append(BOUNDARY * gap)
            offsets.extend([-1] * gap)
            pos += gap
        rd = _SegmentReader(alphabet, k, on_unknown)
        if mode == "chars":
            rd.read_chars(seg)
        elif mode == "words":
            rd.read_words(seg)
        elif mode == "fixed":
            rd.read_fixed(seg, width, frame)
        else:
            raise ValueError(f"unknown mode {mode!r}")
        out_text.append("".join(rd.units))
        offsets.extend(rd.offsets)
        segments.append(Segment(k, pos, pos + len(rd.units), len(seg)))
        skipped.extend(rd.skipped)
        dropped.extend(rd.dropped)
        removal_points.extend(pos + r for r in rd.removals)
        pos += len(rd.units)

    result = Tokenized(alphabet=alphabet, text="".join(out_text), segments=segments,
                       offsets=offsets, skipped=skipped, dropped=dropped,
                       removal_points=removal_points)
    _warn(result)
    return result


def _warn(result: Tokenized):
    name = result.alphabet.name
    if result.skipped:
        counts = Counter(s.text for s in result.skipped)
        shown = ", ".join(f"{t!r} x{n}" for t, n in counts.most_common(5))
        more = "" if len(counts) <= 5 else f" and {len(counts) - 5} other(s)"
        first = result.skipped[0]
        warnings.warn(
            f"{len(result.skipped)} unknown symbol(s) skipped for alphabet {name!r}: "
            f"{shown}{more}; first at offset {first.offset} of segment {first.segment}. "
            "Units on either side of a removal are now adjacent "
            "(see Tokenized.removal_points).",
            UnknownSymbolWarning, stacklevel=3)
    partial = [d for d in result.dropped if d.reason == "partial"]
    if partial:
        warnings.warn(
            f"{len(partial)} stretch(es) did not divide into whole units for alphabet "
            f"{name!r}; the partial unit was dropped (first at offset {partial[0].offset} "
            f"of segment {partial[0].segment}).",
            PartialUnitWarning, stacklevel=3)

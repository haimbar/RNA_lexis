"""Statistical scoring helpers for RNA-Lexis-Stat.

This module keeps the original exact xmotif discovery intact and adds the
post-analysis layer discussed during the reviewer-response work:

* exact motifs are scored against a transcript-specific Markov background;
* candidate cores are filtered/ranked with enrichment p-values and FDR;
* mutation-tolerant Hamming families are accepted only when the expanded
  family remains enriched under the same background.

Every function that depends on the alphabet takes an optional ``alphabet``.
Left out, it is the nucleotide alphabet and behavior is as it has always been.
Given an :class:`rna_lexis.alphabet.Alphabet`, texts and motifs are *internal
strings* (one character per unit, from :func:`rna_lexis.tokenizers.tokenize`),
the alphabet's declared size is used for every pseudocount and count of
possible words whether or not all symbols occur, motifs may not span a
segment boundary (``"_"``), and wildcards in patterns are
:data:`rna_lexis.alphabet.WILDCARD`.  See ``calibration/README.md`` for how the
defaults for the Markov order and the multiple-testing family were chosen.
"""

from __future__ import annotations

import bisect
import csv
import itertools
import math
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Iterable, Iterator, Sequence

import numpy as np


from scipy.stats import binom, poisson

from rna_lexis.algorithms import cores, cover, find_boundary
from rna_lexis.alphabet import Alphabet, BOUNDARY, WILDCARD

try:
    from scipy.stats import false_discovery_control
except ImportError:  # pragma: no cover - compatibility for older SciPy
    false_discovery_control = None


ALPHABET = "acgt"

MAX_FLANK = 500
"""Most units on each side of a motif in a unit profile; the table has one column per position, and no figure can show more."""


def default_markov_order(alphabet_size: int, n_letters: int) -> int:
    """Default Markov order for an alphabet other than the nucleotide one.

    Order 1 when the text has at least one observation per transition on
    average (``n_letters >= alphabet_size ** 2``), otherwise order 0; never
    above 1.  Chosen by simulation (``calibration/README.md``): codons
    (K = 64, up to 1,500 codons) need order 0, letter alphabets (K about 26,
    1,000 letters and more) need order 1, and higher orders are not
    calibrated with a flat pseudocount.
    """

    return 1 if n_letters >= alphabet_size ** 2 else 0


def _require(name: str, value, *, low=None, high=None, low_open=False, high_open=False,
             integer=False):
    """Check a numeric parameter and return it; raise ``ValueError`` saying what is wrong.

    Without this a bad value (``alpha=5``, ``markov_order="a"``, a negative count) either failed deep
    inside the computation with a message about something else, or ran and returned meaningless numbers.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError(f"{name} must be a number, got {value!r}")
    if integer and float(value) != int(value):
        raise ValueError(f"{name} must be a whole number, got {value!r}")
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if low is not None and (value <= low if low_open else value < low):
        raise ValueError(f"{name} must be {'greater than' if low_open else 'at least'} {low}, got {value!r}")
    if high is not None and (value >= high if high_open else value > high):
        raise ValueError(f"{name} must be {'less than' if high_open else 'at most'} {high}, got {value!r}")
    return int(value) if integer else value


def _resolve_order(markov_order: int | None, alphabet, n_letters: int) -> int:
    if markov_order is not None:
        order = _require("markov_order", markov_order, low=0, integer=True)
        size = len(alphabet.chars) if isinstance(alphabet, Alphabet) else len(ALPHABET)
        if order and size ** order > max(n_letters, 1):
            warnings.warn(
                f"markov_order={order} gives {size}**{order} contexts for a text of {n_letters} "
                "units: most contexts are never seen, and the background is mostly the pseudocount",
                UserWarning, stacklevel=3)
        return order
    if isinstance(alphabet, Alphabet):
        return default_markov_order(alphabet.size, n_letters)
    return 1


_CHECKED: dict = {}
_CHECKED_MIN = 2000      # shorter texts are cheaper to check again than to remember
_CHECKED_KEEP = 4


def normalize_sequence(seq: str, alphabet=None) -> str:
    """Normalize a sequence for analysis.

    Without ``alphabet``: RNA/DNA text to the lower-case DNA alphabet (``u``
    becomes ``t``; anything else is dropped), as always.  With an
    :class:`Alphabet`: ``seq`` is an internal string; it is checked, not
    filtered, and the segment boundary is kept.  A character that is not an
    internal character of the alphabet raises ``ValueError``, so text that has
    not been read with :func:`rna_lexis.tokenizers.tokenize` is caught
    instead of silently losing symbols.

    Scoring looks up many motifs in the same long text and each lookup normalizes the text
    again.  The result for the last few long texts is remembered (by identity, so a different
    string is always checked afresh), which makes the repeated checks free.
    """

    if not isinstance(seq, str):
        raise TypeError(f"a sequence must be a string, got {type(seq).__name__}")
    if len(seq) >= _CHECKED_MIN:
        hit = _CHECKED.get((id(seq), id(alphabet)))
        if hit is not None and hit[0] is seq and hit[1] is alphabet:
            return hit[2]
        out = _normalize(seq, alphabet)
        if len(_CHECKED) >= _CHECKED_KEEP:
            _CHECKED.pop(next(iter(_CHECKED)))
        _CHECKED[(id(seq), id(alphabet))] = (seq, alphabet, out)    # holding seq keeps its id from being reused
        return out
    return _normalize(seq, alphabet)


def _normalize(seq: str, alphabet) -> str:
    if isinstance(alphabet, Alphabet):
        s = seq if alphabet.case_sensitive else seq.lower()
        if alphabet.mode == "chars" and alphabet.aliases:
            table = {ord(k): v for k, v in alphabet.aliases.items() if len(k) == 1}
            s = s.translate(table)
        bad = set(s) - set(alphabet.chars) - {BOUNDARY}
        if bad:
            raise ValueError(
                f"{sorted(bad)[:3]!r} is not an internal character of alphabet {alphabet.name!r}; "
                "read the text with rna_lexis.tokenizers.tokenize first")
        return s
    chars = ALPHABET if alphabet is None else alphabet
    s = seq.lower().replace("u", "t")
    if set(s) <= set(chars):
        return s
    return "".join(c for c in s if c in chars)


@dataclass(frozen=True)
class MarkovBackground:
    """Transcript-specific Markov background with pseudocount smoothing.

    ``alphabet`` is the nucleotide characters (default), a string of other
    characters, or an :class:`Alphabet`.  With an ``Alphabet`` the sequence is
    an internal string whose segment boundaries (``"_"``) are respected:
    boundaries are not counted as symbols, no transition is counted across one,
    and no word can occupy a window that spans one.  After construction
    ``alphabet`` is the string of internal characters, whatever was passed, so
    ``len(alphabet)`` is the declared alphabet size.
    """

    sequence: str
    order: int = 1
    pseudocount: float = 1.0
    alphabet: object = ALPHABET

    def __post_init__(self) -> None:
        if self.order < 0:
            raise ValueError("order must be non-negative")
        _require("pseudocount", self.pseudocount, low=0, low_open=True)
        alph_obj = self.alphabet if isinstance(self.alphabet, Alphabet) else None
        chars = alph_obj.chars if alph_obj is not None else self.alphabet
        seq = normalize_sequence(self.sequence, self.alphabet)
        segments = seq.split(BOUNDARY) if BOUNDARY in seq else [seq]
        letters = "".join(segments)
        object.__setattr__(self, "sequence", seq)
        object.__setattr__(self, "alphabet", chars)
        object.__setattr__(self, "_alphabet_obj", alph_obj)
        object.__setattr__(self, "_segment_lengths", tuple(len(s) for s in segments))
        object.__setattr__(self, "_n", len(letters))
        object.__setattr__(self, "_char_counts", Counter(letters))
        object.__setattr__(self, "_context_counts", Counter())
        object.__setattr__(self, "_transition_counts", Counter())

        if self.order > 0:
            context_counts = Counter()
            transition_counts = Counter()
            for segment in segments:
                for i in range(len(segment) - self.order):
                    ctx = segment[i : i + self.order]
                    nxt = segment[i + self.order]
                    context_counts[ctx] += 1
                    transition_counts[(ctx, nxt)] += 1
            object.__setattr__(self, "_context_counts", context_counts)
            object.__setattr__(self, "_transition_counts", transition_counts)

    def _norm(self, word: str) -> str:
        return normalize_sequence(word, self._alphabet_obj if self._alphabet_obj is not None else self.alphabet)

    def slots(self, length: int) -> int:
        """Number of positions where a word of ``length`` fits inside one segment."""

        return sum(max(0, n - length + 1) for n in self._segment_lengths)

    def probability(self, word: str) -> float:
        """Return Markov probability for one word."""

        word = self._norm(word)
        if not word or BOUNDARY in word:
            return 0.0

        alpha = len(self.alphabet)
        pc = self.pseudocount
        n = self._n
        p = (self._char_counts.get(word[0], 0) + pc) / (n + pc * alpha)
        if len(word) == 1:
            return p

        if self.order == 0:
            for char in word[1:]:
                p *= (self._char_counts.get(char, 0) + pc) / (n + pc * alpha)
            return p

        for i in range(1, len(word)):
            if i < self.order:
                # Short prefix while the requested context is not yet available.
                p *= (self._char_counts.get(word[i], 0) + pc) / (n + pc * alpha)
                continue
            ctx = word[i - self.order : i]
            den = self._context_counts.get(ctx, 0) + pc * alpha
            num = self._transition_counts.get((ctx, word[i]), 0) + pc
            p *= num / den
        return p

    def expected_count(self, word: str) -> float:
        """Expected overlapping count of word in this transcript."""

        word = self._norm(word)
        return self.slots(len(word)) * self.probability(word)


def literal_positions(seq: str, txt: str, alphabet=None) -> list[int]:
    """Return overlapping exact start positions for seq in txt.

    Positions index the normalized text.  With an :class:`Alphabet` that is the
    internal string itself (segment boundaries included), and an occurrence
    never spans a boundary because ``seq`` cannot contain one.
    """

    seq = normalize_sequence(seq, alphabet)
    txt = normalize_sequence(txt, alphabet)
    if BOUNDARY in seq:
        raise ValueError("a motif cannot contain the segment boundary")
    if not seq or len(seq) > len(txt):
        return []
    out = []
    start = 0
    while True:
        pos = txt.find(seq, start)
        if pos < 0:
            return out
        out.append(pos)
        start = pos + 1


def window_hamming(a: str, b: str) -> int:
    """Hamming distance between equal-length strings."""

    return sum(x != y for x, y in zip(a, b))


def hamming_family_positions(motif: str, txt: str, radius: int, alphabet=None) -> list[tuple[int, str, int]]:
    """Return all transcript windows within Hamming radius of motif.

    Windows that span a segment boundary are skipped.  Each hit is ``(start, window, distance)``, in order of start.
    Mismatches are counted one motif position at a time over all windows at once, so a long text costs
    ``len(motif)`` passes instead of ``len(text) * len(motif)`` steps.
    """

    motif = normalize_sequence(motif, alphabet)
    txt = normalize_sequence(txt, alphabet)
    L = len(motif)
    if radius < 0 or L == 0 or L > len(txt):
        return []
    n_windows = len(txt) - L + 1
    codes = np.frombuffer(txt.encode("utf-32-le"), dtype=np.uint32)
    distance = np.zeros(n_windows, dtype=np.int32)
    spans_boundary = np.zeros(n_windows, dtype=bool)
    boundary = ord(BOUNDARY)
    for j, ch in enumerate(motif):
        column = codes[j : j + n_windows]
        distance += column != ord(ch)
        spans_boundary |= column == boundary
    keep = np.flatnonzero((distance <= radius) & ~spans_boundary)
    return [(int(i), txt[i : i + L], int(distance[i])) for i in keep]


def nonoverlap_count(positions: Sequence[int], length: int) -> int:
    """Greedy non-overlapping occurrence count."""

    count = 0
    last_end = -1
    for pos in sorted(positions):
        if pos >= last_end:
            count += 1
            last_end = pos + length
    return count


def union_coverage_bp(positions: Sequence[int], length: int) -> int:
    """Number of transcript bases covered by at least one occurrence."""

    intervals = sorted((p, p + length) for p in positions)
    if not intervals:
        return 0
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
    total += cur_end - cur_start
    return total


def bh_adjust(pvalues: Sequence[float], m: int | float | None = None) -> list[float]:
    """Benjamini-Hochberg FDR correction.

    ``m`` is the total number of tests when it is larger than the number of
    p-values given (the missing tests are taken to have p = 1).  Left out, it is
    the number of p-values given, as usual.  Use it to correct over a search
    space that was only partly listed, e.g. all ``K ** k`` possible words of a
    length when only the observed ones have a p-value.
    """

    pvalues = [float(p) for p in pvalues]
    if not pvalues:
        return []
    if m is not None:
        if m < len(pvalues):
            raise ValueError("m cannot be smaller than the number of p-values")
        order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
        adjusted = [1.0] * len(pvalues)
        prev = math.inf
        for rank in range(len(order), 0, -1):
            i = order[rank - 1]
            prev = min(prev, pvalues[i] * m / rank)
            adjusted[i] = min(1.0, prev)
        return adjusted
    if false_discovery_control is not None:
        return [float(x) for x in false_discovery_control(pvalues, method="bh")]

    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [1.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        true_rank = n - rank + 1
        val = min(prev, pvalues[i] * n / true_rank)
        adjusted[i] = min(1.0, val)
        prev = val
    return adjusted


def _safe_enrichment(observed: int, expected: float) -> float:
    if expected <= 0:
        return math.inf if observed > 0 else 0.0
    return observed / expected


def xmotif_occurrence_intervals(txt: str, xmotifs: Iterable[str], alphabet=None) -> list[tuple[int, int, str]]:
    """Return intervals for every exact xmotif occurrence in txt."""

    intervals = []
    for xm in dict.fromkeys(normalize_sequence(x, alphabet) for x in xmotifs if x):
        L = len(xm)
        for pos in literal_positions(xm, txt, alphabet):
            intervals.append((pos, pos + L, xm))
    intervals.sort(key=lambda x: (x[0], x[1], x[2]))
    return intervals


def _inside_any_interval(start: int, end: int, intervals: Sequence[tuple[int, int, str]]) -> bool:
    return any(i0 <= start and end <= i1 for i0, i1, _ in intervals)


class _IntervalIndex:
    """Answers "does some interval contain ``[start, end]``?" without looking at every interval.

    The intervals are sorted by their start; for each prefix of that order the largest end is kept.  An
    interval that begins at or before ``start`` and ends at or after ``end`` exists exactly when the largest
    end among the intervals beginning at or before ``start`` reaches ``end``, so one binary search answers.
    That is the same answer as :func:`_inside_any_interval`, which compares against all of them (scoring a
    few thousand motif positions against a few thousand xmotif intervals took tens of seconds).
    """

    def __init__(self, intervals: Sequence[tuple[int, int, str]]):
        ordered = sorted((i0, i1) for i0, i1, _ in intervals)
        self._starts = [i0 for i0, _ in ordered]
        self._max_end = list(itertools.accumulate((i1 for _, i1 in ordered), max))

    def contains(self, start: int, end: int) -> bool:
        k = bisect.bisect_right(self._starts, start)
        return k > 0 and self._max_end[k - 1] >= end


def score_exact_motifs(
    txt: str,
    motifs: Iterable[str],
    *,
    xmotifs: Iterable[str] | None = None,
    markov_order: int | None = None,
    pseudocount: float = 1.0,
    alpha: float = 0.05,
    enrichment_threshold: float = 10.0,
    min_count: int = 2,
    min_xmotif_type_support: int = 0,
    area_power: float = 1.2,
    alphabet=None,
    family: str | None = None,
) -> list[dict]:
    """Score exact motif candidates with transcript-specific Markov/FDR support.

    Args:
        alphabet:     ``None`` (nucleotides, the default) or an
                      :class:`Alphabet`; with one, ``txt``, ``motifs`` and
                      ``xmotifs`` are internal strings.
        markov_order: Order of the background.  Left out it is 1 for
                      nucleotides and :func:`default_markov_order` for an
                      ``Alphabet`` (0 for codons, 1 for letters).
        family:       The multiple-testing family for the Benjamini-Hochberg
                      correction.  ``"candidates"`` corrects over the motifs
                      given, which are usually already selected for being
                      repeated, so it counts too few tests on a long text.
                      ``"bylength"`` corrects within each motif length ``L``
                      over all ``K ** L`` possible words of that length
                      (``K`` the declared alphabet size); it keeps false
                      positives at or below the nominal level in simulation.
                      Default: ``"candidates"`` for nucleotides, which keeps
                      existing results, and ``"bylength"`` for an ``Alphabet``.
    """

    if family is None:
        family = "candidates" if alphabet is None else "bylength"
    if family not in ("candidates", "bylength"):
        raise ValueError("family must be 'candidates' or 'bylength'")
    _require("alpha", alpha, low=0, high=1, low_open=True)
    _require("enrichment_threshold", enrichment_threshold, low=0)
    _require("min_count", min_count, low=1, integer=True)
    _require("min_xmotif_type_support", min_xmotif_type_support, low=0, integer=True)
    txt = normalize_sequence(txt, alphabet)
    motif_list = [m for m in dict.fromkeys(normalize_sequence(m, alphabet) for m in motifs if m)]
    if any(BOUNDARY in m for m in motif_list):
        raise ValueError("a motif cannot contain the segment boundary")
    order = _resolve_order(markov_order, alphabet, len(txt.replace(BOUNDARY, "")))
    bg = MarkovBackground(txt, order=order, pseudocount=pseudocount,
                          alphabet=alphabet if alphabet is not None else ALPHABET)
    xmotif_list = [normalize_sequence(x, alphabet) for x in (xmotifs or []) if x]
    xmotif_intervals = _IntervalIndex(xmotif_occurrence_intervals(txt, xmotif_list, alphabet)) if xmotif_list else None
    xmotif_set = set(xmotif_list)

    rows = []
    pvalues = []
    for motif in motif_list:
        L = len(motif)
        positions = literal_positions(motif, txt, alphabet)
        observed = len(positions)
        expected = bg.expected_count(motif)
        pvalue = float(poisson.sf(observed - 1, expected)) if expected > 0 else (0.0 if observed > 0 else 1.0)
        pvalues.append(pvalue)

        if xmotif_list:
            x_support = sum(1 for xm in xmotif_set if motif in xm)
            inside = sum(1 for p in positions if xmotif_intervals.contains(p, p + L))
        else:
            x_support = 0
            inside = 0
        outside = observed - inside

        rows.append({
            "motif": _shown(motif, alphabet),
            "length": L,
            "exact_count": observed,
            "nonoverlap_count": nonoverlap_count(positions, L),
            "coverage_bp": union_coverage_bp(positions, L),
            "area_score": nonoverlap_count(positions, L) * (L ** area_power),
            "xmotif_type_support": x_support,
            "inside_xmotif_count": inside,
            "outside_xmotif_count": outside,
            "expected_markov": expected,
            "enrichment_markov": _safe_enrichment(observed, expected),
            "p_markov": pvalue,
            "q_markov": 1.0,
            "markov_order": order,
        })

    if family == "candidates":
        qvalues = bh_adjust(pvalues)
    else:
        qvalues = [1.0] * len(pvalues)
        size = len(bg.alphabet)
        by_length: dict[int, list[int]] = defaultdict(list)
        for idx, row in enumerate(rows):
            by_length[row["length"]].append(idx)
        for length, idxs in by_length.items():
            adjusted = bh_adjust([pvalues[i] for i in idxs], m=float(size) ** length)
            for i, q in zip(idxs, adjusted):
                qvalues[i] = q
    for row, qvalue in zip(rows, qvalues):
        row["q_markov"] = qvalue
        row["statistically_supported"] = (
            row["exact_count"] >= min_count
            and row["q_markov"] <= alpha
            and row["enrichment_markov"] >= enrichment_threshold
            and row["xmotif_type_support"] >= min_xmotif_type_support
        )
        if not xmotif_list:
            row["core_class"] = "unclassified"
        elif row["statistically_supported"] and row["outside_xmotif_count"] == 0:
            row["core_class"] = "family_specific"
        elif row["statistically_supported"]:
            row["core_class"] = "generalized"
        else:
            row["core_class"] = "below_threshold"

    rows.sort(key=lambda r: (
        not r["statistically_supported"],
        r["q_markov"],
        -r["enrichment_markov"],
        -r["coverage_bp"],
        -r["length"],
        r["motif"],
    ))
    for idx, row in enumerate(rows, 1):
        row["rank_statistical"] = idx
    for idx, row in enumerate(sorted(rows, key=lambda r: (-r["coverage_bp"], r["q_markov"], r["motif"])), 1):
        row["rank_coverage"] = idx
    return rows


def enumerate_core_candidates(
    xmotifs: Iterable[str],
    *,
    candidate_min_len: int = 5,
    candidate_max_len: int = 18,
    min_xmotif_type_support: int = 2,
    alphabet=None,
) -> list[str]:
    """Enumerate substrings shared by at least two distinct xmotif sequences."""

    _require("candidate_min_len", candidate_min_len, low=1, integer=True)
    _require("candidate_max_len", candidate_max_len, low=candidate_min_len, integer=True)
    xmotif_list = [normalize_sequence(x, alphabet) for x in dict.fromkeys(xmotifs) if x]
    support: dict[str, set[str]] = defaultdict(set)
    for xm in xmotif_list:
        upper = min(candidate_max_len, len(xm))
        for L in range(candidate_min_len, upper + 1):
            for i in range(len(xm) - L + 1):
                support[xm[i : i + L]].add(xm)

    return sorted(
        (motif for motif, xs in support.items() if len(xs) >= min_xmotif_type_support),
        key=lambda m: (-len(m), m),
    )


def rank_core_candidates(
    txt: str,
    xmotifs: Iterable[str],
    *,
    candidate_min_len: int = 5,
    candidate_max_len: int = 18,
    min_xmotif_type_support: int = 2,
    alphabet=None,
    **score_kwargs,
) -> list[dict]:
    """Generate, score, and rank Markov-supported core motif candidates."""

    candidates = enumerate_core_candidates(
        xmotifs,
        candidate_min_len=candidate_min_len,
        candidate_max_len=candidate_max_len,
        min_xmotif_type_support=min_xmotif_type_support,
        alphabet=alphabet,
    )
    return score_exact_motifs(
        txt,
        candidates,
        xmotifs=xmotifs,
        min_xmotif_type_support=min_xmotif_type_support,
        alphabet=alphabet,
        **score_kwargs,
    )


def _shown(s: str, alphabet) -> str:
    """Text of a motif in a result row.

    Nucleotide strings are upper-cased, as they always have been.  Over an
    :class:`Alphabet` the internal string is kept as it is, so a row can be fed
    back into the functions here; :func:`write_rows_csv` decodes it for the file.
    """

    return s.upper() if alphabet is None else s


def _chars_of(alphabet) -> str:
    """The characters of an alphabet given as None (nucleotides), a string, or an Alphabet."""

    if alphabet is None:
        return ALPHABET
    if isinstance(alphabet, Alphabet):
        return alphabet.chars
    return alphabet


def _variant_count(length: int, radius: int, size: int = len(ALPHABET)) -> int:
    """Number of words within Hamming radius ``radius`` of a word of ``length`` over ``size`` symbols."""

    return sum(math.comb(length, d) * (size - 1) ** d for d in range(radius + 1))


def hamming_variants(motif: str, radius: int, *, max_variants: int = 500_000, alphabet=None) -> Iterator[str]:
    """Yield all words within Hamming radius of motif.

    The number of variants grows with the alphabet size: at radius 2 a word of
    6 codons has 15 * 63 ** 2 = 59,535 variants beyond the first-radius ones, so
    ``max_variants`` matters much more over 64 symbols than over 4.
    """

    motif = normalize_sequence(motif, alphabet)
    symbols = _chars_of(alphabet)
    if BOUNDARY in motif:
        raise ValueError("a motif cannot contain the segment boundary")
    count = _variant_count(len(motif), radius, len(symbols))
    if count > max_variants:
        raise ValueError(
            f"Hamming family has {count:,} variants; lower radius or raise max_variants"
        )
    yield motif
    for d in range(1, radius + 1):
        for positions in itertools.combinations(range(len(motif)), d):
            replacements = [
                [base for base in symbols if base != motif[pos]]
                for pos in positions
            ]
            for repl in itertools.product(*replacements):
                chars = list(motif)
                for pos, base in zip(positions, repl):
                    chars[pos] = base
                yield "".join(chars)


def hamming_family_probability(
    motif: str,
    radius: int,
    background: MarkovBackground,
    *,
    max_variants: int = 500_000,
) -> float:
    """Markov probability of the full Hamming family."""

    alphabet = background._alphabet_obj if background._alphabet_obj is not None else background.alphabet
    return sum(background.probability(v)
               for v in hamming_variants(motif, radius, max_variants=max_variants, alphabet=alphabet))


def mutation_family_tests(
    txt: str,
    motifs: Iterable[str],
    *,
    mutr: float = 1 / 6,
    M: int = 4,
    markov_order: int | None = None,
    pseudocount: float = 1.0,
    alpha: float = 0.05,
    enrichment_threshold: float = 5.0,
    expected_max: float = 5.0,
    min_family_count: int = 2,
    max_variants: int = 500_000,
    alphabet=None,
) -> list[dict]:
    """Evaluate Hamming-radius motif families with Markov enrichment and FDR.

    With an :class:`Alphabet` the texts are internal strings and a mismatch is a
    different unit.  The correction is over the rows listed (each motif at each
    radius), as before; that family has not been calibrated for other alphabets.
    """

    _require("mutr", mutr, low=0, high=1)
    _require("M", M, low=0, integer=True)
    _require("alpha", alpha, low=0, high=1, low_open=True)
    _require("enrichment_threshold", enrichment_threshold, low=0)
    _require("expected_max", expected_max, low=0, low_open=True)
    _require("min_family_count", min_family_count, low=1, integer=True)
    _require("max_variants", max_variants, low=1, integer=True)
    txt = normalize_sequence(txt, alphabet)
    motif_list = [m for m in dict.fromkeys(normalize_sequence(m, alphabet) for m in motifs if m)]
    markov_order = _resolve_order(markov_order, alphabet, len(txt.replace(BOUNDARY, "")))
    bg = MarkovBackground(txt, order=markov_order, pseudocount=pseudocount,
                          alphabet=alphabet if alphabet is not None else ALPHABET)
    size = len(bg.alphabet)
    rows = []
    pvalues = []

    for motif in motif_list:
        L = len(motif)
        max_radius = min(M, math.floor(L * mutr))
        for radius in range(max_radius + 1):
            hits = hamming_family_positions(motif, txt, radius, alphabet)
            hist = Counter(dist for _, _, dist in hits)
            family_count = len(hits)
            exact_count = hist.get(0, 0)
            approx_count = family_count - exact_count
            mismatch_burden = sum(dist * cnt for dist, cnt in hist.items())
            pmut = mismatch_burden / (L * family_count) if family_count else 0.0
            try:
                family_probability = hamming_family_probability(
                    motif, radius, bg, max_variants=max_variants
                )
                expected = bg.slots(L) * family_probability
                pvalue = float(poisson.sf(family_count - 1, expected)) if expected > 0 else (
                    0.0 if family_count > 0 else 1.0
                )
                variant_count = _variant_count(L, radius, size)
                variant_error = ""
            except ValueError as exc:
                expected = math.nan
                pvalue = 1.0
                variant_count = _variant_count(L, radius, size)
                variant_error = str(exc)

            pvalues.append(pvalue)
            rows.append({
                "motif": _shown(motif, alphabet),
                "length": L,
                "radius": radius,
                "exact_count": exact_count,
                "approx_count": approx_count,
                "family_count": family_count,
                "mismatch_histogram": ";".join(f"{d}:{hist.get(d, 0)}" for d in range(radius + 1)),
                "mismatch_burden": mismatch_burden,
                "pmut": pmut,
                "p_stable": float(binom.cdf(mismatch_burden, L * family_count, mutr)) if family_count else 1.0,
                "expected_family_markov": expected,
                "enrichment_family": _safe_enrichment(family_count, expected) if not math.isnan(expected) else math.nan,
                "p_family": pvalue,
                "q_family": 1.0,
                "markov_order": markov_order,
                "hamming_variant_count": variant_count,
                "variant_error": variant_error,
                "decision": "below_threshold",
            })

    qvalues = bh_adjust(pvalues)
    for row, qvalue in zip(rows, qvalues):
        row["q_family"] = qvalue
        accepted = (
            row["family_count"] >= min_family_count
            and row["q_family"] <= alpha
            and not math.isnan(row["enrichment_family"])
            and row["enrichment_family"] >= enrichment_threshold
            and row["expected_family_markov"] <= expected_max
        )
        if accepted and row["radius"] > 0 and row["approx_count"] > 0:
            row["decision"] = "mutation_supported"
        elif accepted:
            row["decision"] = "exact_or_specific_only"

    rows.sort(key=lambda r: (
        r["motif"],
        not (r["decision"] == "mutation_supported"),
        r["radius"],
        r["q_family"],
        -r["enrichment_family"] if not math.isnan(r["enrichment_family"]) else math.inf,
    ))
    return rows


def _decoded_row(row: dict, alphabet) -> dict:
    out = dict(row)
    for key in _TEXT_COLUMNS:
        if isinstance(out.get(key), str):
            out[key] = alphabet.show(out[key], upper=True)
    if "pattern" in out and "left_anchor" in row and "right_anchor" in row:
        out["pattern"] = (f"{out['left_anchor']}[gap:{row.get('min_gap')}-{row.get('max_gap')}]"
                          f"{out['right_anchor']}")
    return out


def best_mutation_family_per_motif(rows: Iterable[dict]) -> list[dict]:
    """Select the strongest accepted radius per motif for compact reporting."""

    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[row["motif"]].append(row)

    best = []
    decision_rank = {"mutation_supported": 0, "exact_or_specific_only": 1, "below_threshold": 2}
    for motif, motif_rows in grouped.items():
        motif_rows.sort(key=lambda r: (
            decision_rank.get(r["decision"], 9),
            -r["radius"],
            r["q_family"],
            -r["family_count"],
        ))
        best.append(motif_rows[0])
    best.sort(key=lambda r: (decision_rank.get(r["decision"], 9), r["q_family"], r["motif"]))
    return best


_TEXT_COLUMNS = ("motif", "left_anchor", "right_anchor", "gap_sequence", "matched_sequence", "core")


_FORMULA_STARTS = ("=", "+", "-", "@", "\t", "\r")


def csv_safe(value):
    """A cell value that a spreadsheet will show as text and never run as a formula.

    A cell that starts with ``=``, ``+``, ``-`` or ``@`` is a formula to Excel, LibreOffice and
    Google Sheets, so a symbol label from an alphabet file (``=HYPERLINK(...)``) written into a results table
    could run when the table is opened.  Such a string gets a leading apostrophe.  Numbers, and strings
    that are numbers (``-3.2``), are left alone, and so are the letters-only sequences of nucleotide text.
    """
    if isinstance(value, str) and value.startswith(_FORMULA_STARTS):
        try:
            float(value)
        except ValueError:
            return "'" + value
    return value


def write_rows_csv(path: str, rows: Sequence[dict], alphabet=None) -> str:
    """Write dict rows to CSV and return the path.

    With an :class:`Alphabet` the columns that hold sequences (``motif``,
    ``left_anchor``, ``right_anchor``, ``gap_sequence``, ``matched_sequence``,
    ``core`` and the ``pattern`` of a gapped motif) are decoded from internal
    strings to the alphabet's labels in upper case, e.g. ``ATG AAA TTT`` for
    codons.  Without one the rows are written as they are.
    """

    if alphabet is not None:
        rows = [_decoded_row(r, alphabet) for r in rows]
    if not rows:
        with open(path, "w", newline="", encoding="utf-8") as handle:
            handle.write("")
        return path

    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows({k: csv_safe(v) for k, v in row.items()} for row in rows)
    return path


def _allowed_pattern_chars(char: str, alphabet: str = ALPHABET) -> list[str]:
    return list(alphabet) if char in {"n", ".", "?", "*"} else [char]


def markov_pattern_probability(pattern: str, background: MarkovBackground) -> float:
    """Probability of a pattern with wildcard positions under a Markov model.

    Nucleotide background: fixed bases are written as A/C/G/T/U. Wildcards may
    be N, '.', '?', or '*' and match any nucleotide. This is used for gap-aware
    motifs such as ``TGTNNNTATA`` without enumerating every gap sequence.

    Background over an :class:`Alphabet`: the pattern is an internal string and
    the only wildcard is :data:`rna_lexis.alphabet.WILDCARD` (``n`` and ``.`` are
    ordinary symbols there); any other character that is not a symbol of the
    alphabet raises ``ValueError``.
    """

    if background._alphabet_obj is None:
        raw = pattern.lower().replace("u", "t")
        pattern = "".join(c if c in background.alphabet else "n" for c in raw)
        wildcards = {"n", ".", "?", "*"}
    else:
        bad = sorted({c for c in pattern if c not in background.alphabet and c != WILDCARD})
        if bad:
            raise ValueError(
                f"{bad[:3]!r} is neither a symbol of alphabet {background._alphabet_obj.name!r} "
                "nor the wildcard")
        wildcards = {WILDCARD}
    if not pattern:
        return 0.0

    def allowed(char: str) -> list[str]:
        return list(background.alphabet) if char in wildcards else [char]

    alpha = len(background.alphabet)
    pc = background.pseudocount
    n = background._n

    def base_prob(base: str) -> float:
        return (background._char_counts.get(base, 0) + pc) / (n + pc * alpha)

    def transition_prob(ctx: str, base: str) -> float:
        den = background._context_counts.get(ctx, 0) + pc * alpha
        num = background._transition_counts.get((ctx, base), 0) + pc
        return num / den

    order = background.order
    states: dict[str, float] = {}
    for base in allowed(pattern[0]):
        key = base[-order:] if order > 0 else ""
        states[key] = states.get(key, 0.0) + base_prob(base)

    for idx, char in enumerate(pattern[1:], 1):
        next_states: dict[str, float] = {}
        for state, state_prob in states.items():
            for base in allowed(char):
                if order == 0 or idx < order:
                    step_prob = base_prob(base)
                else:
                    step_prob = transition_prob(state[-order:], base)
                key = (state + base)[-order:] if order > 0 else ""
                next_states[key] = next_states.get(key, 0.0) + state_prob * step_prob
        states = next_states
    return sum(states.values())


def find_gapped_motif_hits(
    txt: str,
    left_anchor: str,
    right_anchor: str,
    *,
    min_gap: int = 0,
    max_gap: int = 30,
    alphabet=None,
) -> list[dict]:
    """Find exact anchor pairs separated by an unconstrained gap.

    A hit never spans a segment boundary.
    """

    txt = normalize_sequence(txt, alphabet)
    left = normalize_sequence(left_anchor, alphabet)
    right = normalize_sequence(right_anchor, alphabet)
    if not left or not right:
        raise ValueError("left_anchor and right_anchor must contain sequence letters")
    if min_gap < 0 or max_gap < min_gap:
        raise ValueError("gap range must satisfy 0 <= min_gap <= max_gap")

    hits = []
    for start in literal_positions(left, txt, alphabet):
        for gap in range(min_gap, max_gap + 1):
            right_start = start + len(left) + gap
            right_end = right_start + len(right)
            if right_end > len(txt):
                continue
            if txt[right_start:right_end] == right and BOUNDARY not in txt[start:right_end]:
                hits.append({
                    "start": start,
                    "end": right_end,
                    "left_anchor": _shown(left, alphabet),
                    "right_anchor": _shown(right, alphabet),
                    "gap_length": gap,
                    "gap_sequence": _shown(txt[start + len(left):right_start], alphabet),
                    "matched_sequence": _shown(txt[start:right_end], alphabet),
                })
    return hits


def score_gapped_motif(
    txt: str,
    left_anchor: str,
    right_anchor: str,
    *,
    min_gap: int = 0,
    max_gap: int = 30,
    markov_order: int | None = None,
    pseudocount: float = 1.0,
    alphabet=None,
) -> dict:
    """Score an anchor-gap-anchor motif as one statistical family."""

    txt = normalize_sequence(txt, alphabet)
    left = normalize_sequence(left_anchor, alphabet)
    right = normalize_sequence(right_anchor, alphabet)
    hits = find_gapped_motif_hits(
        txt,
        left,
        right,
        min_gap=min_gap,
        max_gap=max_gap,
        alphabet=alphabet,
    )
    markov_order = _resolve_order(markov_order, alphabet, len(txt.replace(BOUNDARY, "")))
    bg = MarkovBackground(txt, order=markov_order, pseudocount=pseudocount,
                          alphabet=alphabet if alphabet is not None else ALPHABET)
    wildcard = "n" if alphabet is None else WILDCARD
    expected = 0.0
    for gap in range(min_gap, max_gap + 1):
        pattern = left + (wildcard * gap) + right
        expected += bg.slots(len(pattern)) * markov_pattern_probability(pattern, bg)
    observed = len(hits)
    pvalue = float(poisson.sf(observed - 1, expected)) if expected > 0 else (0.0 if observed > 0 else 1.0)
    positions = [hit["start"] for hit in hits]
    lengths = [hit["end"] - hit["start"] for hit in hits]
    coverage = union_coverage_variable(positions, lengths)
    return {
        "pattern": f"{_shown(left, alphabet)}[gap:{min_gap}-{max_gap}]{_shown(right, alphabet)}",
        "left_anchor": _shown(left, alphabet),
        "right_anchor": _shown(right, alphabet),
        "min_gap": min_gap,
        "max_gap": max_gap,
        "observed_count": observed,
        "coverage_bp": coverage,
        "expected_markov": expected,
        "enrichment_markov": _safe_enrichment(observed, expected),
        "p_markov": pvalue,
        "q_markov": pvalue,
        "markov_order": markov_order,
    }


def unit_profile(
    txt: str,
    motif: str,
    *,
    flank: int = 5,
    mismatches: int = 0,
    alphabet=None,
    group_by=None,
    alpha: float = 0.01,
) -> dict:
    """What surrounds a motif: how often each symbol occurs at every position around its occurrences.

    The windows are the occurrences of ``motif`` (exact, or within ``mismatches``
    substitutions) with ``flank`` units on each side.  A window that would run off
    the text or span a segment boundary is left out.  For each position and symbol
    the result holds the count, the frequency, and the **log-odds against the
    text's own background**: ``log2(f / b)`` where ``b`` is the symbol's frequency
    in the whole text (the order-0 background the statistics use, with its
    pseudocount) and ``f`` is its frequency at that position, shrunk toward ``b``
    by one pseudo-observation so that a few occurrences cannot give extreme values.

    Log-odds alone would mislead with few windows: a symbol expected 0.1 times that
    appears once looks ninefold enriched, and one that never appears looks depleted,
    though neither is evidence of anything.  So every cell also gets a two-sided
    binomial p-value (the count against ``n`` windows at the background
    frequency), and ``significant`` marks the cells with ``p <= alpha``; the plot
    colors only those.  Depletion shows up only where a symbol was expected often
    enough to notice its absence.

    This is the view that replaces a sequence logo for alphabets of about 20 symbols
    or more: a logo measures information against a uniform background, so common
    symbols look informative when they are only common, and its height is capped at
    ``log2(K)`` bits, too little to read.  Log-odds against the background shows
    only what is enriched or depleted next to the motif, and scales to any size.

    Args:
        txt:        The text (an internal string when ``alphabet`` is given).
        motif:      The motif.
        flank:      Units on each side of the motif.
        mismatches: Hamming radius for the occurrences (0 = exact).
        alphabet:   ``None`` (nucleotides) or an :class:`Alphabet`.
        group_by:   Optional :class:`Projection` of the alphabet (for codons,
                    ``codon_to_amino_acid()``); rows are ordered by group, so
                    synonymous codons sit together.
        alpha:      A cell is ``significant`` when its p-value is at most this
                    (default 0.01).  The plot rings the strongest enriched cell of a
                    column only if it also passes a Bonferroni correction over all
                    cells.

    Returns:
        A dict with ``chars`` (the symbol rows, internal characters), ``labels``
        (the same for display), ``columns`` (position relative to the motif start,
        ``-flank`` to ``len(motif) + flank - 1``), ``motif_columns`` (the slice of
        columns that is the motif), ``counts``, ``frequency``, ``log_odds``,
        ``p_value`` and ``significant`` (lists of rows, one per symbol),
        ``alpha``, ``background``, ``n`` (windows used),
        ``skipped`` (occurrences left out), ``groups`` (``(label, first row, end row)``)
        and ``motif`` / ``motif_labels`` (the motif, internal and for display).
    """

    _require("flank", flank, low=0, high=MAX_FLANK, integer=True)
    _require("mismatches", mismatches, low=0, integer=True)
    _require("alpha", alpha, low=0, high=1, low_open=True)
    txt_n = normalize_sequence(txt, alphabet)
    motif_n = normalize_sequence(motif, alphabet)
    if not motif_n or BOUNDARY in motif_n:
        raise ValueError("the motif must be non-empty and cannot contain the segment boundary")
    L = len(motif_n)
    bg = MarkovBackground(txt_n, order=0, alphabet=alphabet if alphabet is not None else ALPHABET)
    chars = list(bg.alphabet)
    if group_by is not None and alphabet is not None:
        by_target = group_by.groups()
        labels_in_order = [lab for tgt in sorted(by_target) for lab in sorted(by_target[tgt])]
        chars = [alphabet.char(lab) for lab in labels_in_order]
    index = {c: i for i, c in enumerate(chars)}
    W = flank + L + flank
    counts = [[0] * W for _ in chars]

    starts = [p for p, _, _ in hamming_family_positions(motif_n, txt_n, mismatches, alphabet)]
    used = skipped = 0
    for p in starts:
        a, b = p - flank, p + L + flank
        if a < 0 or b > len(txt_n) or BOUNDARY in txt_n[a:b]:
            skipped += 1
            continue
        used += 1
        for j, c in enumerate(txt_n[a:b]):
            counts[index[c]][j] += 1

    background = [bg.probability(c) for c in chars]
    total = sum(background)
    background = [x / total for x in background]
    frequency, log_odds, p_value, significant = [], [], [], []
    for i in range(len(chars)):
        frow = [(counts[i][j] + background[i]) / (used + 1) for j in range(W)]
        frequency.append(frow)
        log_odds.append([math.log2(f / background[i]) for f in frow])
        prow = []
        for j in range(W):
            c, b = counts[i][j], background[i]
            if used == 0:
                prow.append(1.0)
            elif c >= used * b:
                prow.append(float(binom.sf(c - 1, used, b)))      # at least this many
            else:
                prow.append(float(binom.cdf(c, used, b)))         # at most this many
        p_value.append(prow)
        significant.append([p <= alpha for p in prow])

    groups = []
    if group_by is not None and alphabet is not None:
        start = 0
        for tgt in sorted(group_by.groups()):
            n = len(group_by.groups()[tgt])
            groups.append((tgt, start, start + n))
            start += n
    labels = [alphabet.show(c) if alphabet is not None else c for c in chars]
    return {"chars": chars, "labels": labels, "columns": list(range(-flank, L + flank)),
            "motif_columns": (flank, flank + L), "counts": counts, "frequency": frequency,
            "log_odds": log_odds, "p_value": p_value, "significant": significant, "alpha": alpha,
            "background": background, "n": used, "skipped": skipped,
            "groups": groups, "motif": motif_n,
            "motif_labels": [alphabet.show(c) if alphabet is not None else c for c in motif_n]}


def write_unit_profile_csv(path: str, profile: dict) -> str:
    """Write a :func:`unit_profile` as a table: one row per symbol, and for every
    position its count and its log-odds against the background."""

    header = ["symbol", "background"]
    for col in profile["columns"]:
        header += [f"count@{col}", f"log2_odds@{col}", f"p@{col}"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for i, label in enumerate(profile["labels"]):
            row = [csv_safe(label.upper()), f"{profile['background'][i]:.6g}"]
            for j in range(len(profile["columns"])):
                row += [profile["counts"][i][j], f"{profile['log_odds'][i][j]:.3f}",
                        f"{profile['p_value'][i][j]:.3g}"]
            writer.writerow(row)
    return path


def union_coverage_variable(positions: Sequence[int], lengths: Sequence[int]) -> int:
    """Coverage for intervals with different lengths."""

    intervals = sorted((p, p + l) for p, l in zip(positions, lengths))
    if not intervals:
        return 0
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
    total += cur_end - cur_start
    return total


def spacing_periodicity_test(
    positions: Sequence[int],
    seq_len: int,
) -> dict:
    """Test whether a motif's occurrences are more regularly spaced than chance.

    Two complementary tests are applied:

    **Gap cluster test** — estimates the candidate period T as the median
    consecutive gap, then counts how many gaps fall within δ = max(5, 5%·T)
    of T.  Under H0 each gap is Uniform[0, seq_len], so the count follows
    Bin(n_gaps, 2δ/seq_len).  The p-value is Bonferroni-corrected by n_gaps
    because T is estimated from the same data.

    **Rayleigh test (circular statistics)** — maps each position modulo T to
    an angle θ = 2π·(pos mod T)/T and tests whether the angles cluster on
    the circle.  The mean resultant length R ∈ [0,1] measures concentration;
    Z = m·R² is the test statistic.  A small p-value confirms a dominant
    periodic signal at period T (Mardia & Jupp, 2000 approximation).

    Parameters
    ----------
    positions : sequence of int
        Sorted or unsorted start positions of exact motif occurrences.
    seq_len : int
        Total length of the sequence (used as the gap null range).

    Returns
    -------
    dict with keys:
        m, gaps, mean_gap, std_gap, period,
        delta, k_near_T, gaps_near_T, p_cluster,
        rayleigh_r, rayleigh_z, p_rayleigh.
    """
    positions = sorted(positions)
    m = len(positions)
    gaps = [positions[i + 1] - positions[i] for i in range(m - 1)]

    mean_gap = sum(gaps) / len(gaps)
    std_gap = math.sqrt(sum((g - mean_gap) ** 2 for g in gaps) / len(gaps))

    # Candidate period: median of consecutive gaps
    sorted_gaps = sorted(gaps)
    n_g = len(sorted_gaps)
    if n_g % 2 == 1:
        period = sorted_gaps[n_g // 2]
    else:
        period = int(round((sorted_gaps[n_g // 2 - 1] + sorted_gaps[n_g // 2]) / 2.0))
    period = max(period, 1)

    # Gap cluster test: how many consecutive gaps fall within δ of T?
    # Under H0 each gap is Uniform[0, seq_len], so P(gap ∈ [T-δ, T+δ]) = 2δ/seq_len.
    # Bonferroni-corrected by n_g because T is estimated from the same data.
    delta = max(5, round(0.05 * period))
    gaps_near_T = [g for g in gaps if abs(g - period) <= delta]
    k_near_T = len(gaps_near_T)
    p_window = min(1.0, 2 * delta / seq_len)
    p_cluster_raw = float(1.0 - binom.cdf(k_near_T - 1, n_g, p_window))
    p_cluster = min(1.0, n_g * p_cluster_raw)

    # Rayleigh test at the candidate period
    TWO_PI = 2.0 * math.pi
    angles = [TWO_PI * (p % period) / period for p in positions]
    C = sum(math.cos(a) for a in angles) / m
    S = sum(math.sin(a) for a in angles) / m
    R = math.sqrt(C * C + S * S)
    Z = m * R * R
    p_rayleigh = math.exp(-Z) * (1.0 + (2.0 * Z - Z * Z) / (4.0 * m))
    p_rayleigh = max(0.0, min(1.0, p_rayleigh))

    return {
        'm': m,
        'gaps': gaps,
        'mean_gap': mean_gap,
        'std_gap': std_gap,
        'period': period,
        'delta': delta,
        'k_near_T': k_near_T,
        'gaps_near_T': gaps_near_T,
        'p_cluster': p_cluster,
        'rayleigh_r': R,
        'rayleigh_z': Z,
        'p_rayleigh': p_rayleigh,
    }


def shared_exact_motifs(
    txt_a: str,
    txt_b: str,
    cores_a: Sequence[str] | None = None,
    min_len: int = 6,
    n_top: int = 6,
    minxmlen: int = 7,
    maxxmlen: int = 60,
    mincorelen: int = 6,
    mincount: int = 2,
    alphabet=None,
) -> list[str]:
    """Select motifs with exact occurrences in both sequences.

    Ported from GeneDB `similarity.py::main()`'s `arc_candidates`
    construction: cores discovered jointly on the *concatenation*
    `txt_a + '_' + txt_b` (the '_' separator prevents spurious motifs that
    straddle the boundary) that have >= 1 exact hit in each half come
    first -- joint discovery surfaces shared motifs that neither sequence's
    own single-sequence core discovery necessarily finds on its own, since
    a motif's cores depend on its own recurrence pattern within that one
    sequence. `cores_a`, if given (e.g. sequence A's already-computed
    `corelist`), supplements this with any additional A-cores that also
    have an exact hit in B and weren't already found by the joint scan.
    All candidates are ranked by weighted coverage in sequence A
    (`algorithms.cover()`) descending.

    Note: only the shared-motif *list* is reproduced here, not GeneDB's
    concat-IoU/fold-enrichment/Markov-null statistics computed from the
    same joint discovery -- those numbers are a separate, out-of-scope
    analysis; this function only needs the motifs themselves.

    Args:
        txt_a:      Sequence A (used for the coverage ranking).
        txt_b:      Sequence B (a candidate must have >= 1 exact hit here).
        cores_a:    Optional pre-computed cores of sequence A (e.g. from an
                    already-parsed session) to supplement the joint scan.
        min_len:    Minimum motif length to consider (default 6).
        n_top:      Maximum number of motifs to return (default 6).
        minxmlen:   Passed to find_boundary() for the joint discovery scan.
        maxxmlen:   Passed to find_boundary() for the joint discovery scan.
        mincorelen: Passed to cores() for the joint discovery scan.
        mincount:   Passed to find_boundary() for the joint discovery scan.
        alphabet:   ``None`` (nucleotides) or an :class:`Alphabet`; with one,
                    both sequences are internal strings.

    Returns:
        Up to `n_top` motif strings, longest-coverage-in-A first.
    """
    concat_txt = f'{txt_a}_{txt_b}'
    concat_xmotifs = find_boundary(concat_txt, minlen=minxmlen, maxlen=maxxmlen,
                                    at_least=mincount)
    concat_cores = cores(concat_txt, concat_xmotifs, minclen=mincorelen)

    candidates = [c for c in concat_cores
                  if literal_positions(c, txt_a, alphabet) and literal_positions(c, txt_b, alphabet)]
    seen = set(candidates)
    for c in (cores_a or []):
        if c not in seen and literal_positions(c, txt_b, alphabet):
            candidates.append(c)
            seen.add(c)

    candidates = [c for c in candidates if len(c) >= min_len]
    candidates.sort(key=lambda m: cover(m, txt_a), reverse=True)
    return candidates[:n_top]

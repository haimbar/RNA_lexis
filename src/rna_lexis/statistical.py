"""Statistical scoring helpers for RNA-Lexis-Stat.

This module keeps the original exact xmotif discovery intact and adds the
post-analysis layer discussed during the reviewer-response work:

* exact motifs are scored against a transcript-specific Markov background;
* candidate cores are filtered/ranked with enrichment p-values and FDR;
* mutation-tolerant Hamming families are accepted only when the expanded
  family remains enriched under the same background.
"""

from __future__ import annotations

import csv
import itertools
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Iterable, Iterator, Sequence

from scipy.stats import binom, poisson

try:
    from scipy.stats import false_discovery_control
except ImportError:  # pragma: no cover - compatibility for older SciPy
    false_discovery_control = None


ALPHABET = "acgt"


@dataclass(frozen=True)
class MarkovBackground:
    """Transcript-specific Markov background with pseudocount smoothing."""

    sequence: str
    order: int = 1
    pseudocount: float = 1.0
    alphabet: str = ALPHABET

    def __post_init__(self) -> None:
        if self.order < 0:
            raise ValueError("order must be non-negative")
        seq = normalize_sequence(self.sequence)
        object.__setattr__(self, "sequence", seq)
        object.__setattr__(self, "_char_counts", Counter(seq))
        object.__setattr__(self, "_context_counts", Counter())
        object.__setattr__(self, "_transition_counts", Counter())

        if self.order > 0:
            context_counts = Counter()
            transition_counts = Counter()
            for i in range(len(seq) - self.order):
                ctx = seq[i : i + self.order]
                nxt = seq[i + self.order]
                context_counts[ctx] += 1
                transition_counts[(ctx, nxt)] += 1
            object.__setattr__(self, "_context_counts", context_counts)
            object.__setattr__(self, "_transition_counts", transition_counts)

    def probability(self, word: str) -> float:
        """Return Markov probability for one word."""

        word = normalize_sequence(word)
        if not word:
            return 0.0

        alpha = len(self.alphabet)
        pc = self.pseudocount
        n = len(self.sequence)
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

        word = normalize_sequence(word)
        slots = max(0, len(self.sequence) - len(word) + 1)
        return slots * self.probability(word)


def normalize_sequence(seq: str) -> str:
    """Normalize RNA/DNA text to lower-case DNA alphabet."""

    return "".join(c for c in seq.lower().replace("u", "t") if c in ALPHABET)


def literal_positions(seq: str, txt: str) -> list[int]:
    """Return overlapping exact start positions for seq in txt."""

    seq = normalize_sequence(seq)
    txt = normalize_sequence(txt)
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


def hamming_family_positions(motif: str, txt: str, radius: int) -> list[tuple[int, str, int]]:
    """Return all transcript windows within Hamming radius of motif."""

    motif = normalize_sequence(motif)
    txt = normalize_sequence(txt)
    L = len(motif)
    if radius < 0 or L == 0 or L > len(txt):
        return []
    hits = []
    for i in range(len(txt) - L + 1):
        window = txt[i : i + L]
        dist = window_hamming(window, motif)
        if dist <= radius:
            hits.append((i, window, dist))
    return hits


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


def bh_adjust(pvalues: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction."""

    pvalues = [float(p) for p in pvalues]
    if not pvalues:
        return []
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


def xmotif_occurrence_intervals(txt: str, xmotifs: Iterable[str]) -> list[tuple[int, int, str]]:
    """Return intervals for every exact xmotif occurrence in txt."""

    intervals = []
    for xm in dict.fromkeys(normalize_sequence(x) for x in xmotifs if x):
        L = len(xm)
        for pos in literal_positions(xm, txt):
            intervals.append((pos, pos + L, xm))
    intervals.sort(key=lambda x: (x[0], x[1], x[2]))
    return intervals


def _inside_any_interval(start: int, end: int, intervals: Sequence[tuple[int, int, str]]) -> bool:
    return any(i0 <= start and end <= i1 for i0, i1, _ in intervals)


def score_exact_motifs(
    txt: str,
    motifs: Iterable[str],
    *,
    xmotifs: Iterable[str] | None = None,
    markov_order: int = 1,
    pseudocount: float = 1.0,
    alpha: float = 0.05,
    enrichment_threshold: float = 10.0,
    min_count: int = 2,
    min_xmotif_type_support: int = 0,
    area_power: float = 1.2,
) -> list[dict]:
    """Score exact motif candidates with transcript-specific Markov/FDR support."""

    txt = normalize_sequence(txt)
    motif_list = [m for m in dict.fromkeys(normalize_sequence(m) for m in motifs if m)]
    bg = MarkovBackground(txt, order=markov_order, pseudocount=pseudocount)
    xmotif_list = [normalize_sequence(x) for x in (xmotifs or []) if x]
    xmotif_intervals = xmotif_occurrence_intervals(txt, xmotif_list) if xmotif_list else []

    rows = []
    pvalues = []
    for motif in motif_list:
        L = len(motif)
        positions = literal_positions(motif, txt)
        observed = len(positions)
        expected = bg.expected_count(motif)
        pvalue = float(poisson.sf(observed - 1, expected)) if expected > 0 else (0.0 if observed > 0 else 1.0)
        pvalues.append(pvalue)

        if xmotif_list:
            x_support = sum(1 for xm in set(xmotif_list) if motif in xm)
            inside = sum(1 for p in positions if _inside_any_interval(p, p + L, xmotif_intervals))
        else:
            x_support = 0
            inside = 0
        outside = observed - inside

        rows.append({
            "motif": motif.upper(),
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
            "markov_order": markov_order,
        })

    qvalues = bh_adjust(pvalues)
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
) -> list[str]:
    """Enumerate substrings shared by at least two distinct xmotif sequences."""

    xmotif_list = [normalize_sequence(x) for x in dict.fromkeys(xmotifs) if x]
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
    **score_kwargs,
) -> list[dict]:
    """Generate, score, and rank Markov-supported core motif candidates."""

    candidates = enumerate_core_candidates(
        xmotifs,
        candidate_min_len=candidate_min_len,
        candidate_max_len=candidate_max_len,
        min_xmotif_type_support=min_xmotif_type_support,
    )
    return score_exact_motifs(
        txt,
        candidates,
        xmotifs=xmotifs,
        min_xmotif_type_support=min_xmotif_type_support,
        **score_kwargs,
    )


def _variant_count(length: int, radius: int) -> int:
    return sum(math.comb(length, d) * (len(ALPHABET) - 1) ** d for d in range(radius + 1))


def hamming_variants(motif: str, radius: int, *, max_variants: int = 500_000) -> Iterator[str]:
    """Yield all words within Hamming radius of motif."""

    motif = normalize_sequence(motif)
    count = _variant_count(len(motif), radius)
    if count > max_variants:
        raise ValueError(
            f"Hamming family has {count:,} variants; lower radius or raise max_variants"
        )
    yield motif
    for d in range(1, radius + 1):
        for positions in itertools.combinations(range(len(motif)), d):
            replacements = [
                [base for base in ALPHABET if base != motif[pos]]
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

    return sum(background.probability(v) for v in hamming_variants(motif, radius, max_variants=max_variants))


def mutation_family_tests(
    txt: str,
    motifs: Iterable[str],
    *,
    mutr: float = 1 / 6,
    M: int = 4,
    markov_order: int = 1,
    pseudocount: float = 1.0,
    alpha: float = 0.05,
    enrichment_threshold: float = 5.0,
    expected_max: float = 5.0,
    min_family_count: int = 2,
    max_variants: int = 500_000,
) -> list[dict]:
    """Evaluate Hamming-radius motif families with Markov enrichment and FDR."""

    txt = normalize_sequence(txt)
    motif_list = [m for m in dict.fromkeys(normalize_sequence(m) for m in motifs if m)]
    bg = MarkovBackground(txt, order=markov_order, pseudocount=pseudocount)
    rows = []
    pvalues = []

    for motif in motif_list:
        L = len(motif)
        max_radius = min(M, math.floor(L * mutr))
        for radius in range(max_radius + 1):
            hits = hamming_family_positions(motif, txt, radius)
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
                expected = max(0, len(txt) - L + 1) * family_probability
                pvalue = float(poisson.sf(family_count - 1, expected)) if expected > 0 else (
                    0.0 if family_count > 0 else 1.0
                )
                variant_count = _variant_count(L, radius)
                variant_error = ""
            except ValueError as exc:
                expected = math.nan
                pvalue = 1.0
                variant_count = _variant_count(L, radius)
                variant_error = str(exc)

            pvalues.append(pvalue)
            rows.append({
                "motif": motif.upper(),
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


def write_rows_csv(path: str, rows: Sequence[dict]) -> str:
    """Write dict rows to CSV and return the path."""

    if not rows:
        with open(path, "w", newline="", encoding="utf-8") as handle:
            handle.write("")
        return path

    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _allowed_pattern_chars(char: str, alphabet: str = ALPHABET) -> list[str]:
    return list(alphabet) if char in {"n", ".", "?", "*"} else [char]


def markov_pattern_probability(pattern: str, background: MarkovBackground) -> float:
    """Probability of a pattern with wildcard positions under a Markov model.

    Fixed bases are written as A/C/G/T/U. Wildcards may be N, '.', '?', or '*'
    and match any nucleotide. This is used for gap-aware motifs such as
    ``TGTNNNTATA`` without enumerating every gap sequence.
    """

    raw = pattern.lower().replace("u", "t")
    pattern = "".join(c if c in background.alphabet else "n" for c in raw)
    if not pattern:
        return 0.0

    alpha = len(background.alphabet)
    pc = background.pseudocount
    n = len(background.sequence)

    def base_prob(base: str) -> float:
        return (background._char_counts.get(base, 0) + pc) / (n + pc * alpha)

    def transition_prob(ctx: str, base: str) -> float:
        den = background._context_counts.get(ctx, 0) + pc * alpha
        num = background._transition_counts.get((ctx, base), 0) + pc
        return num / den

    order = background.order
    states: dict[str, float] = {}
    for base in _allowed_pattern_chars(pattern[0], background.alphabet):
        key = base[-order:] if order > 0 else ""
        states[key] = states.get(key, 0.0) + base_prob(base)

    for idx, char in enumerate(pattern[1:], 1):
        next_states: dict[str, float] = {}
        for state, state_prob in states.items():
            for base in _allowed_pattern_chars(char, background.alphabet):
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
) -> list[dict]:
    """Find exact anchor pairs separated by an unconstrained gap."""

    txt = normalize_sequence(txt)
    left = normalize_sequence(left_anchor)
    right = normalize_sequence(right_anchor)
    if not left or not right:
        raise ValueError("left_anchor and right_anchor must contain sequence letters")
    if min_gap < 0 or max_gap < min_gap:
        raise ValueError("gap range must satisfy 0 <= min_gap <= max_gap")

    hits = []
    for start in literal_positions(left, txt):
        for gap in range(min_gap, max_gap + 1):
            right_start = start + len(left) + gap
            right_end = right_start + len(right)
            if right_end > len(txt):
                continue
            if txt[right_start:right_end] == right:
                hits.append({
                    "start": start,
                    "end": right_end,
                    "left_anchor": left.upper(),
                    "right_anchor": right.upper(),
                    "gap_length": gap,
                    "gap_sequence": txt[start + len(left):right_start].upper(),
                    "matched_sequence": txt[start:right_end].upper(),
                })
    return hits


def score_gapped_motif(
    txt: str,
    left_anchor: str,
    right_anchor: str,
    *,
    min_gap: int = 0,
    max_gap: int = 30,
    markov_order: int = 1,
    pseudocount: float = 1.0,
) -> dict:
    """Score an anchor-gap-anchor motif as one statistical family."""

    txt = normalize_sequence(txt)
    left = normalize_sequence(left_anchor)
    right = normalize_sequence(right_anchor)
    hits = find_gapped_motif_hits(
        txt,
        left,
        right,
        min_gap=min_gap,
        max_gap=max_gap,
    )
    bg = MarkovBackground(txt, order=markov_order, pseudocount=pseudocount)
    expected = 0.0
    for gap in range(min_gap, max_gap + 1):
        pattern = left + ("n" * gap) + right
        expected += max(0, len(txt) - len(pattern) + 1) * markov_pattern_probability(pattern, bg)
    observed = len(hits)
    pvalue = float(poisson.sf(observed - 1, expected)) if expected > 0 else (0.0 if observed > 0 else 1.0)
    positions = [hit["start"] for hit in hits]
    lengths = [hit["end"] - hit["start"] for hit in hits]
    coverage = union_coverage_variable(positions, lengths)
    return {
        "pattern": f"{left.upper()}[gap:{min_gap}-{max_gap}]{right.upper()}",
        "left_anchor": left.upper(),
        "right_anchor": right.upper(),
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

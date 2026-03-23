"""Sequence alignment with affine gap penalties (Gotoh global and local)."""

from dataclasses import dataclass
from typing import Tuple

NEG_INF = -10**15  # sufficiently small for typical scoring ranges


@dataclass
class AlignmentResult:
    mode: str                 # "global" or "local"
    aligned_a: str
    aligned_b: str
    markers: str
    score: int
    matches: int
    mismatches: int
    gaps: int


def score_sub(a: str, b: str, match: int, mismatch: int) -> int:
    """Return the substitution score for aligning character a against character b."""
    return match if a == b else mismatch


def make_markers(aligned_a: str, aligned_b: str) -> Tuple[str, int, int, int]:
    """Build an alignment marker string and summary counts.

    For each aligned column:
    * ``'|'`` — exact match
    * ``'.'`` — mismatch (both characters are non-gap)
    * ``' '`` — at least one side is a gap (``'-'``)

    Args:
        aligned_a: First aligned sequence (may contain ``'-'`` for gaps).
        aligned_b: Second aligned sequence (same length as aligned_a).

    Returns:
        Tuple (markers, matches, mismatches, gaps) where markers is the
        column-by-column symbol string and the counts sum to len(aligned_a).
    """
    marks = []
    matches = mismatches = gaps = 0
    for ca, cb in zip(aligned_a, aligned_b):
        if ca == "-" or cb == "-":
            marks.append(" ")
            gaps += 1
        elif ca == cb:
            marks.append("|")
            matches += 1
        else:
            marks.append(".")
            mismatches += 1
    return "".join(marks), matches, mismatches, gaps


def print_alignment(res: AlignmentResult, width: int = 80) -> None:
    """Print a formatted alignment to stdout.

    Outputs the mode header, score summary, and the two aligned sequences
    with their marker line in blocks of width columns.

    Args:
        res:   AlignmentResult from gotoh_global() or gotoh_local().
        width: Number of alignment columns per printed block (default 80).
    """
    print(f"{res.mode.upper()} ALIGNMENT (affine gaps)")
    print(f"Score: {res.score} | matches: {res.matches} | mismatches: {res.mismatches} | gaps: {res.gaps}\n")
    for i in range(0, len(res.aligned_a), width):
        print(res.aligned_a[i:i+width])
        print(res.markers[i:i+width])
        print(res.aligned_b[i:i+width])
        print()


def gotoh_global(a: str, b: str, match: int = 2, mismatch: int = -1,
                 gap_open: int = -5, gap_extend: int = -1) -> AlignmentResult:
    """
    Global alignment with affine gap penalties (Gotoh algorithm).

    Gap penalty model:
        gap(k) = gap_open + (k - 1) * gap_extend, for k >= 1
    """
    n, m = len(a), len(b)

    # DP matrices
    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]  # gap in B (A aligned to '-')
    Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]  # gap in A ('-' aligned to B)

    # Traceback pointers: store prev state (0=M,1=X,2=Y), or -1
    ptrM = [[-1] * (m + 1) for _ in range(n + 1)]
    ptrX = [[-1] * (m + 1) for _ in range(n + 1)]
    ptrY = [[-1] * (m + 1) for _ in range(n + 1)]

    # init
    M[0][0] = 0
    ptrM[0][0] = -1

    # Leading gaps: only X along first column, only Y along first row
    for i in range(1, n + 1):
        X[i][0] = gap_open + (i - 1) * gap_extend
        ptrX[i][0] = 1  # from X (extend) after first; doesn't matter much for edge
    for j in range(1, m + 1):
        Y[0][j] = gap_open + (j - 1) * gap_extend
        ptrY[0][j] = 2

    # Fill
    # Tie-break preference (if equal): M > X > Y (stable/consistent)
    def argmax3(v0, v1, v2) -> Tuple[int, int]:
        """Return (max_value, argmax_index) for three values with M>X>Y tie-breaking."""
        best = v0
        state = 0
        if v1 > best or (v1 == best and state != 0):
            best, state = v1, 1
        if v2 > best or (v2 == best and state not in (0, 1)):
            best, state = v2, 2
        return best, state

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = score_sub(a[i-1], b[j-1], match, mismatch)

            # M: come diagonally from any state
            candM0 = M[i-1][j-1]
            candM1 = X[i-1][j-1]
            candM2 = Y[i-1][j-1]
            best_prev, prev_state = argmax3(candM0, candM1, candM2)
            M[i][j] = best_prev + s
            ptrM[i][j] = prev_state

            # X: gap in B (advance i, keep j)
            # open from M or Y; extend from X
            open_from_M = M[i-1][j] + gap_open
            extend_from_X = X[i-1][j] + gap_extend
            open_from_Y = Y[i-1][j] + gap_open
            bestX, prevX = argmax3(open_from_M, extend_from_X, open_from_Y)
            X[i][j] = bestX
            ptrX[i][j] = prevX

            # Y: gap in A (advance j, keep i)
            open_from_M = M[i][j-1] + gap_open
            open_from_X = X[i][j-1] + gap_open
            extend_from_Y = Y[i][j-1] + gap_extend
            # Here argmax order is (from M, from X, from Y)
            bestY, prevY = argmax3(open_from_M, open_from_X, extend_from_Y)
            Y[i][j] = bestY
            ptrY[i][j] = prevY

    # Choose best ending state at (n,m)
    end_score, end_state = max((M[n][m], 0), (X[n][m], 1), (Y[n][m], 2), key=lambda x: x[0])

    # Traceback
    i, j = n, m
    state = end_state
    out_a, out_b = [], []

    while i > 0 or j > 0:
        if state == 0:  # M: diagonal
            prev = ptrM[i][j]
            out_a.append(a[i-1])
            out_b.append(b[j-1])
            i -= 1
            j -= 1
            state = prev
        elif state == 1:  # X: up (gap in B)
            prev = ptrX[i][j]
            out_a.append(a[i-1])
            out_b.append("-")
            i -= 1
            state = prev
        else:  # state == 2, Y: left (gap in A)
            prev = ptrY[i][j]
            out_a.append("-")
            out_b.append(b[j-1])
            j -= 1
            state = prev

        # At edges, state pointers can be odd due to initialization;
        # the coordinates guarantee progress, so we keep going.

    aligned_a = "".join(reversed(out_a))
    aligned_b = "".join(reversed(out_b))
    markers, matches, mismatches, gaps = make_markers(aligned_a, aligned_b)

    return AlignmentResult(
        mode="global",
        aligned_a=aligned_a,
        aligned_b=aligned_b,
        markers=markers,
        score=end_score,
        matches=matches,
        mismatches=mismatches,
        gaps=gaps
    )


def gotoh_local(a: str, b: str, match: int = 2, mismatch: int = -1,
                gap_open: int = -5, gap_extend: int = -1) -> AlignmentResult:
    """
    Local alignment with affine gap penalties (Gotoh-style Smith–Waterman).

    Gap penalty model:
        gap(k) = gap_open + (k - 1) * gap_extend, for k >= 1

    Local alignment resets negative scores to 0 and traceback stops at 0.
    """
    n, m = len(a), len(b)

    M = [[0] * (m + 1) for _ in range(n + 1)]
    X = [[0] * (m + 1) for _ in range(n + 1)]
    Y = [[0] * (m + 1) for _ in range(n + 1)]

    ptrM = [[-1] * (m + 1) for _ in range(n + 1)]
    ptrX = [[-1] * (m + 1) for _ in range(n + 1)]
    ptrY = [[-1] * (m + 1) for _ in range(n + 1)]

    best_score = 0
    best_i = best_j = 0
    best_state = 0

    def best_with_zero(v0, v1, v2) -> Tuple[int, int]:
        """Return (max_value, argmax_index) clamped to (0, -1) when max < 0 (local-alignment reset)."""
        best = v0
        state = 0
        if v1 > best:
            best, state = v1, 1
        if v2 > best:
            best, state = v2, 2
        if best < 0:
            return 0, -1
        return best, state

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = score_sub(a[i-1], b[j-1], match, mismatch)

            # M
            cand0 = M[i-1][j-1]
            cand1 = X[i-1][j-1]
            cand2 = Y[i-1][j-1]
            best_prev = max(cand0, cand1, cand2)
            new_val = best_prev + s
            if new_val <= 0:
                M[i][j] = 0
                ptrM[i][j] = -1
            else:
                # choose which state for pointer (tie-breaking: prefer M, then X, then Y)
                if cand0 == best_prev:
                    ptrM[i][j] = 0
                elif cand1 == best_prev:
                    ptrM[i][j] = 1
                else:
                    ptrM[i][j] = 2
                M[i][j] = new_val

            # X (gap in B): up
            vM = M[i-1][j] + gap_open
            vX = X[i-1][j] + gap_extend
            vY = Y[i-1][j] + gap_open
            bestX = max(vM, vX, vY)
            if bestX <= 0:
                X[i][j] = 0
                ptrX[i][j] = -1
            else:
                if vM == bestX:
                    ptrX[i][j] = 0
                elif vX == bestX:
                    ptrX[i][j] = 1
                else:
                    ptrX[i][j] = 2
                X[i][j] = bestX

            # Y (gap in A): left
            vM = M[i][j-1] + gap_open
            vX = X[i][j-1] + gap_open
            vY = Y[i][j-1] + gap_extend
            bestY = max(vM, vX, vY)
            if bestY <= 0:
                Y[i][j] = 0
                ptrY[i][j] = -1
            else:
                if vM == bestY:
                    ptrY[i][j] = 0
                elif vX == bestY:
                    ptrY[i][j] = 1
                else:
                    ptrY[i][j] = 2
                Y[i][j] = bestY

            # Track best over all states
            if M[i][j] > best_score:
                best_score, best_i, best_j, best_state = M[i][j], i, j, 0
            if X[i][j] > best_score:
                best_score, best_i, best_j, best_state = X[i][j], i, j, 1
            if Y[i][j] > best_score:
                best_score, best_i, best_j, best_state = Y[i][j], i, j, 2

    # Traceback from best cell until we hit score 0 (ptr = -1)
    i, j, state = best_i, best_j, best_state
    out_a, out_b = [], []

    while i > 0 or j > 0:
        if state == 0:
            if M[i][j] == 0 or ptrM[i][j] == -1:
                break
            prev = ptrM[i][j]
            out_a.append(a[i-1]); out_b.append(b[j-1])
            i -= 1; j -= 1
            state = prev
        elif state == 1:
            if X[i][j] == 0 or ptrX[i][j] == -1:
                break
            prev = ptrX[i][j]
            out_a.append(a[i-1]); out_b.append("-")
            i -= 1
            state = prev
        else:
            if Y[i][j] == 0 or ptrY[i][j] == -1:
                break
            prev = ptrY[i][j]
            out_a.append("-"); out_b.append(b[j-1])
            j -= 1
            state = prev

    aligned_a = "".join(reversed(out_a))
    aligned_b = "".join(reversed(out_b))
    markers, matches, mismatches, gaps = make_markers(aligned_a, aligned_b)

    return AlignmentResult(
        mode="local",
        aligned_a=aligned_a,
        aligned_b=aligned_b,
        markers=markers,
        score=best_score,
        matches=matches,
        mismatches=mismatches,
        gaps=gaps
    )

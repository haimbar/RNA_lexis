"""Simulation harness for choosing the codon-mode background (docs/scoping_codon_units.md, section 7).

The package's statistical scoring works on nucleotides only (statistical.py hardcodes
"acgt"), so this module re-implements, on one-character-per-codon text, exactly the
steps of the CLI's ``rank-cores`` pipeline that matter:

    xmotifs   = find_boundary(text, minlen, maxlen, at_least)         (the package's own function)
    candidates= substrings shared by >= 2 distinct xmotifs             (enumerate_core_candidates, generic)
    score     = Poisson upper tail of the observed count against the expected count under a
                background model, Benjamini-Hochberg across the candidates, then the same
                "statistically supported" rule (count, q, enrichment, xmotif support).

Only the *background model* is swapped.  ``selfcheck()`` shows that, on nucleotide text, the
generic scoring and the generic order-1 background reproduce ``score_exact_motifs`` and
``MarkovBackground`` exactly.
"""

from __future__ import annotations

import contextlib
import io
import itertools
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from rna_lexis.alphabet import codon                      # noqa: E402
from rna_lexis.algorithms import find_boundary            # noqa: E402
from rna_lexis.statistical import bh_adjust, _safe_enrichment   # noqa: E402

ALPHA = codon()
NCOD = ALPHA.size                                         # 64
CHARS = ALPHA.chars                                       # internal character of codon i
IDX = {c: i for i, c in enumerate(CHARS)}
LABELS = ALPHA.symbols                                    # 'aaa', 'aac', ... (index = 16*b1 + 4*b2 + b3, a<c<g<t)

# standard genetic code, indexed by our codon index
_TCAG = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODE = {"".join(p): aa for p, aa in zip(itertools.product("TCAG", repeat=3), _TCAG)}
AA_OF = [_CODE[lab.upper()] for lab in LABELS]
AAS = sorted(set(AA_OF) - {"*"})
AA_INDEX = {a: i for i, a in enumerate(AAS)}
AA_NUM = np.array([AA_INDEX.get(a, -1) for a in AA_OF])   # amino-acid index of each codon (-1 = stop)
NT_POS = np.array([[LABELS[c][j] for j in range(3)] for c in range(NCOD)])   # nucleotide at each codon position
NT_IDX = {"a": 0, "c": 1, "g": 2, "t": 3}
NT_OF = np.array([[NT_IDX[ch] for ch in LABELS[c]] for c in range(NCOD)])    # (64, 3)


def to_text(units) -> str:
    return "".join(CHARS[u] for u in units)


def to_units(word: str):
    return [IDX[c] for c in word]


# ----------------------------------------------------------------------
# data
# ----------------------------------------------------------------------
def load_genes(path, min_codons=1):
    """FASTA of CDS (stop removed) -> {symbol: np.array of codon indices}."""
    genes, cur = {}, None
    for line in open(path, encoding="utf-8"):
        if line.startswith(";") or not line.strip():
            continue
        if line.startswith(">"):
            cur = line[1:].split()[0]
            genes[cur] = []
        else:
            genes[cur].append(line.strip().lower())
    out = {}
    lab_index = {lab: i for i, lab in enumerate(LABELS)}
    for sym, parts in genes.items():
        seq = "".join(parts)
        units = [lab_index[seq[i:i + 3]] for i in range(0, len(seq) - 2, 3)]
        if len(units) >= min_codons:
            out[sym] = np.array(units, dtype=np.int16)
    return out


# ----------------------------------------------------------------------
# background models.  Each factory takes the codon-index array of the text under test and
# returns a Model whose batch(W) gives, for each row of the (n_words, k) index matrix W, the
# probability of that word at one position.  prob(w) is the single-word form.
# ----------------------------------------------------------------------
class Model:
    uses_counts = False

    def __init__(self, batch):
        self.batch = batch

    def prob(self, w):
        return float(self.batch(np.asarray([list(w)]))[0])


def _counts(units, K=NCOD):
    return np.bincount(units, minlength=K).astype(float)


def _nt_position_model(units, pc=0.5):
    """q(codon) = product over the 3 codon positions of the nucleotide frequency at that position."""
    q = np.ones(NCOD)
    for j in range(3):
        f = np.bincount(NT_OF[units, j], minlength=4).astype(float) + pc
        f /= f.sum()
        q *= f[NT_OF[:, j]]
    return q / q.sum()


def _order0(p):
    return Model(lambda W: np.prod(p[W], axis=1))


def _order1(first, trans_p):
    """first: (K,) probabilities of the first symbol; trans_p: (K, K) P(b | a)."""
    def batch(W):
        out = first[W[:, 0]].copy()
        for i in range(1, W.shape[1]):
            out *= trans_p[W[:, i - 1], W[:, i]]
        return out
    return Model(batch)


def _transitions(units, K):
    trans = np.zeros((K, K))
    np.add.at(trans, (units[:-1], units[1:]), 1)
    return trans


def bg_order0_flat(pc, K=NCOD):
    def make(units):
        return _order0((_counts(units, K) + pc) / (len(units) + pc * K))
    make.name = f"order0_flat_pc{pc:g}"
    return make


def bg_ntpos():
    def make(units):
        return _order0(_nt_position_model(units))
    make.name = "ntpos"
    return make


def bg_hybrid(s):
    """Codon usage shrunk toward the nucleotide-by-position model: (n_c + s*q_c) / (N + s)."""
    def make(units):
        q = _nt_position_model(units)
        return _order0((_counts(units) + s * q) / (len(units) + s))
    make.name = f"hybrid_s{s:g}"
    return make


def bg_order1_flat(pc, K=NCOD):
    """Order-1 Markov exactly as MarkovBackground: first symbol from the marginal, then transitions."""
    def make(units):
        n = len(units)
        trans = _transitions(units, K)
        ctx = trans.sum(axis=1)
        first = (_counts(units, K) + pc) / (n + pc * K)
        return _order1(first, (trans + pc) / (ctx[:, None] + pc * K))
    make.name = f"order1_flat_pc{pc:g}"
    return make


def bg_order1_backoff(s, s0=64):
    """Order-1 transitions shrunk toward the (hybrid) order-0 model: (n_ab + s*p0_b) / (n_a + s)."""
    def make(units):
        n = len(units)
        p0 = (_counts(units) + s0 * _nt_position_model(units)) / (n + s0)
        trans = _transitions(units, NCOD)
        ctx = trans.sum(axis=1)
        return _order1(p0, (trans + s * p0[None, :]) / (ctx[:, None] + s))
    make.name = f"order1_backoff_s{s:g}"
    return make


class MaskedModel(Model):
    """Order-0 usage estimated with the tested word's own occurrences removed from the counts.

    A repeated word raises the counts of its own codons and so the expected count it is judged
    against (self-inclusion).  Here each word is scored against usage estimated from the text
    without its occurrences: p_c = (n_c - cnt*m_c(w) + s*q_c) / (N - cnt*k + s), where m_c(w) is
    how many times codon c occurs in the word and cnt how many times the word occurs.
    """
    uses_counts = True

    def __init__(self, n_counts, N, q, s):
        self.n, self.N, self.q, self.s = n_counts, N, q, s

    def batch(self, W, cnt):
        k = W.shape[1]
        P = np.ones(len(W))
        for i in range(k):
            mult = (W == W[:, i:i + 1]).sum(axis=1)
            num = self.n[W[:, i]] - cnt * mult + self.s * self.q[W[:, i]]
            den = self.N - cnt * k + self.s
            P *= np.maximum(num, 1e-12) / den
        return P

    def prob(self, w):
        return float(self.batch(np.asarray([list(w)]), np.array([1]))[0])


def bg_masked_hybrid(s):
    def make(units):
        return MaskedModel(_counts(units), len(units), _nt_position_model(units), s)
    make.name = f"masked_hybrid_s{s:g}"
    return make


def bg_masked_flat(pc):
    def make(units):
        return MaskedModel(_counts(units), len(units), np.full(NCOD, 1.0 / NCOD), pc * NCOD)
    make.name = f"masked_order0_flat_pc{pc:g}"
    return make


def bg_reference(ref_usage):
    """A fixed reference codon usage (e.g. pooled from other genes), not estimated from the text."""
    p = ref_usage / ref_usage.sum()

    def make(units):
        return _order0(p)
    make.name = "reference_usage"
    return make


def bg_oracle(model):
    """The true generating model (an upper bound on what any estimated background can do)."""
    def make(units):
        return model
    make.name = "oracle"
    return make


# ----------------------------------------------------------------------
# discovery + scoring (the pipeline)
# ----------------------------------------------------------------------
def discover(text, minlen, mincand, at_least=2, maxlen=60, cmax=18, support=2):
    """xmotifs from find_boundary, then candidates shared by >= `support` distinct xmotifs."""
    with contextlib.redirect_stdout(io.StringIO()):
        xm = find_boundary(text, minlen=minlen, maxlen=maxlen, at_least=at_least)
    xmotifs = sorted(set(xm), key=lambda x: (-len(x), x))
    sup = defaultdict(set)
    for x in xmotifs:
        for L in range(mincand, min(cmax, len(x)) + 1):
            for i in range(len(x) - L + 1):
                sup[x[i:i + L]].add(x)
    cands = sorted((m for m, xs in sup.items() if len(xs) >= support), key=lambda m: (-len(m), m))
    return xmotifs, cands


def positions(word, text):
    out, start = [], 0
    while True:
        p = text.find(word, start)
        if p < 0:
            return out
        out.append(p)
        start = p + 1


def score(text, candidates, xmotifs, prob, alpha=0.05, enrichment=10.0, min_count=2, support=2):
    """The scoring step of score_exact_motifs with a pluggable per-position word probability."""
    n = len(text)
    xset = set(xmotifs)
    rows, pvals = [], []
    for m in candidates:
        obs = len(positions(m, text))
        slots = max(0, n - len(m) + 1)
        exp = slots * (prob.prob(to_units(m)) if isinstance(prob, Model) else prob(to_units(m)))
        p = float(poisson.sf(obs - 1, exp)) if exp > 0 else (0.0 if obs > 0 else 1.0)
        rows.append({"motif": m, "obs": obs, "exp": exp, "enrich": _safe_enrichment(obs, exp),
                     "x_support": sum(1 for x in xset if m in x), "p": p})
        pvals.append(p)
    for r, q in zip(rows, bh_adjust(pvals)):
        r["q"] = q
        r["fdr_only"] = r["obs"] >= min_count and q <= alpha and r["x_support"] >= support
        r["supported"] = r["fdr_only"] and r["enrich"] >= enrichment
    return rows


# ----------------------------------------------------------------------
# scoring stage only: every repeated word, three multiple-testing families
# ----------------------------------------------------------------------
KS = (2, 3, 4, 5, 6)


def kmer_table(units, ks=KS):
    """Every observed word of each length in ks: (k, W index matrix, counts)."""
    out = []
    for k in ks:
        if len(units) < k:
            continue
        win = np.lib.stride_tricks.sliding_window_view(np.asarray(units, dtype=np.int64), k)
        key = np.zeros(len(win), dtype=np.int64)
        for i in range(k):
            key = key * NCOD + win[:, i]
        uk, cnt = np.unique(key, return_counts=True)
        W = np.empty((len(uk), k), dtype=np.int64)
        rem = uk.copy()
        for i in range(k - 1, -1, -1):
            W[:, i] = rem % NCOD
            rem //= NCOD
        out.append((k, W, cnt))
    return out


def bh_vec(p, m=None):
    """Benjamini-Hochberg q-values for the p-values p, as if there were m tests in all
    (m >= len(p); the missing tests have p = 1)."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    m = n if m is None else m
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


FAMILIES = ("candidates", "observed", "possible", "possible_loose", "bylength", "bylength_loose")


def scoring_stage(units, model, alpha=0.05, enrichment=10.0, min_count=2, ks=KS):
    """Score every repeated word under a background and three multiple-testing families.

    candidates: BH over the words with count >= 2 only, as score_exact_motifs does with its
                candidate list (the list is already selected for being repeated);
    observed:   BH over every word observed at least once (the search space the text offers);
    possible:   BH as if every one of 64**k words of each length had been tested.

    Returns {family: (supported mask, word table)}; the table is (k, W, count, p).
    """
    tab = kmer_table(units, ks)
    ks_ = np.concatenate([np.full(len(c), k) for k, _, c in tab])
    cnt = np.concatenate([c for _, _, c in tab])
    p = np.empty(len(cnt))
    enrich = np.empty(len(cnt))
    pos = 0
    n = len(units)
    for k, W, c in tab:
        slots = max(0, n - k + 1)
        exp = slots * (model.batch(W, c) if model.uses_counts else model.batch(W))
        with np.errstate(divide="ignore", invalid="ignore"):
            pv = np.where(exp > 0, poisson.sf(c - 1, exp), np.where(c > 0, 0.0, 1.0))
            en = np.where(exp > 0, c / exp, np.inf)
        p[pos:pos + len(c)] = pv
        enrich[pos:pos + len(c)] = en
        pos += len(c)
    rep = cnt >= min_count
    res = {}
    q = np.ones(len(cnt))
    q[rep] = bh_vec(p[rep])
    res["candidates"] = rep & (q <= alpha) & (enrich >= enrichment)
    q = bh_vec(p)
    res["observed"] = rep & (q <= alpha) & (enrich >= enrichment)
    m_possible = float(sum(NCOD ** k for k in ks))
    q = bh_vec(p, m=m_possible)
    res["possible"] = rep & (q <= alpha) & (enrich >= enrichment)
    # a looser operating point, where a biased background shows up as false positives
    res["possible_loose"] = rep & (q <= 0.5) & (enrich >= 3.0)
    # corrected within each word length over that length's 64**k possible words
    qk = np.empty(len(cnt))
    pos = 0
    for k, W, c in tab:
        qk[pos:pos + len(c)] = bh_vec(p[pos:pos + len(c)], m=float(NCOD ** k))
        pos += len(c)
    res["bylength"] = rep & (qk <= alpha) & (enrich >= enrichment)
    res["bylength_loose"] = rep & (qk <= 0.5) & (enrich >= 3.0)
    words = [(k, W, c) for k, W, c in tab]
    return res, ks_, cnt, words


# ----------------------------------------------------------------------
# self-check against the package, on nucleotide text
# ----------------------------------------------------------------------
def selfcheck(seed=1):
    """Generic scoring/background == score_exact_motifs / MarkovBackground, on nucleotide text."""
    import random
    from rna_lexis.statistical import MarkovBackground, score_exact_motifs, rank_core_candidates

    rnd = random.Random(seed)
    txt = "".join(rnd.choice("acgt") for _ in range(400))
    for pos in (20, 150, 300):
        txt = txt[:pos] + "gattacagattaca" + txt[pos + 14:]
    nt = {c: i for i, c in enumerate("acgt")}
    units = np.array([nt[c] for c in txt])

    # 0. the vectorised order-1 model (alphabet size 4) is MarkovBackground
    m1 = bg_order1_flat(1.0, K=4)(units)
    # 1. order-1 background identical to MarkovBackground (alphabet size 4)
    bg = MarkovBackground(txt, order=1, pseudocount=1.0)
    n = len(units)
    c = np.bincount(units, minlength=4).astype(float)
    trans = np.zeros((4, 4))
    np.add.at(trans, (units[:-1], units[1:]), 1)
    ctx = trans.sum(axis=1)

    def prob(w):
        p = (c[w[0]] + 1) / (n + 4)
        for a, b in zip(w, w[1:]):
            p *= (trans[a, b] + 1) / (ctx[a] + 4)
        return p
    for word in ("gattaca", "acgt", "tttt", "gattacagatt"):
        w = [nt[ch] for ch in word]
        assert math.isclose((n - len(word) + 1) * prob(w), bg.expected_count(word), rel_tol=1e-12), word
        assert math.isclose((n - len(word) + 1) * m1.prob(w), bg.expected_count(word), rel_tol=1e-12), word

    # 2. generic discovery + scoring == the package's rank_core_candidates (same q-values and flags)
    with contextlib.redirect_stdout(io.StringIO()):
        xm_pkg = sorted(set(find_boundary(txt, minlen=7, maxlen=60, at_least=2)), key=lambda x: (-len(x), x))
    ref = {r["motif"].lower(): r for r in rank_core_candidates(txt, xm_pkg, markov_order=1)}
    cands = sorted(m for m in {x[i:i + L] for x in xm_pkg for L in range(5, min(18, len(x)) + 1)
                               for i in range(len(x) - L + 1)}
                   if sum(1 for x in xm_pkg if m in x) >= 2)
    mine = {}
    for r in _score_generic_nt(txt, cands, xm_pkg, prob, nt):
        mine[r["motif"]] = r
    assert set(mine) == set(ref), (len(mine), len(ref))
    for m, r in mine.items():
        assert math.isclose(r["q"], ref[m]["q_markov"], rel_tol=1e-9, abs_tol=1e-15), m
        assert r["supported"] == ref[m]["statistically_supported"], m
    return len(ref)


def _score_generic_nt(text, candidates, xmotifs, prob, nt):
    n = len(text)
    xset = set(xmotifs)
    rows, pvals = [], []
    for m in candidates:
        obs = len(positions(m, text))
        exp = max(0, n - len(m) + 1) * prob([nt[ch] for ch in m])
        p = float(poisson.sf(obs - 1, exp)) if exp > 0 else (0.0 if obs > 0 else 1.0)
        rows.append({"motif": m, "obs": obs, "enrich": _safe_enrichment(obs, exp),
                     "x_support": sum(1 for x in xset if m in x), "p": p})
        pvals.append(p)
    for r, q in zip(rows, bh_adjust(pvals)):
        r["q"] = q
        r["supported"] = r["obs"] >= 2 and q <= 0.05 and r["enrich"] >= 10.0 and r["x_support"] >= 2
    return rows


if __name__ == "__main__":
    n = selfcheck()
    print(f"selfcheck passed: generic scoring matches the package on {n} nucleotide candidates")

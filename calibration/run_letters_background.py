"""Calibration of the background order for letter alphabets (scoping doc, section 7).

Corpora (real text, read through rna_lexis.alphabet / tokenizers):
  antony   Mark Antony's oration, 26 declared Latin letters (GeneDB/MarkAntony/data/antony_oration.txt)
  homer    the opening of the Odyssey in polytonic Greek, reduced to its 24 base letters
           (diacritics removed, final sigma folded into sigma)

Language is not iid, so the question is which Markov order of background is right at which text
length.  Null texts are drawn from a Markov chain of known order m (0, 1 or 2) fitted to the corpus,
so the truth is known exactly.  Backgrounds are order 0..3 with a flat pseudocount of 1, estimated
from the text under test, exactly as MarkovBackground does.

Measured, per (corpus, true order m, N letters, background order):
  * KL per word (3 letters) and the probability mass of words underestimated 5x;
  * bias of the expected count on the words seen >= 2 times (the ones that can be called significant);
  * false positives at the usual and a loose operating point, family = every possible word;
  * recall of a planted phrase.

Usage: python3 run_letters_background.py [--reps 200] [--power-reps 100] [--quick]
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import os
import sys
import time
import unicodedata
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H                                           # noqa: E402  (bh_vec only)
from rna_lexis.alphabet import Alphabet                        # noqa: E402
from rna_lexis.tokenizers import tokenize                      # noqa: E402

HERE = Path(__file__).resolve().parent
GENEDB = Path(os.environ.get("RNA_LEXIS_GENEDB", "/Users/haim/Projects/GeneDB"))   # for antony_oration.txt
HOMER = Path(os.environ.get("RNA_LEXIS_HOMER", "/Users/haim/Desktop/Homer"))
TEXTS = Path(os.environ.get("RNA_LEXIS_TEXTS", HERE / "data" / "texts"))   # genesis.txt, constitution.txt
KS = (2, 3, 4, 5, 6, 7)
NS = [1000, 3000, 8000, 20000, 50000]
ORDERS = (0, 1, 2, 3)
TRUE_ORDERS = (0, 1, 2)
_G = {}


# ----------------------------------------------------------------------
# corpora
# ----------------------------------------------------------------------
def load_corpora():
    latin = Alphabet(name="Latin letters", symbols=tuple("abcdefghijklmnopqrstuvwxyz"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ant = tokenize(open(GENEDB / "MarkAntony/data/antony_oration.txt", encoding="utf-8").read().lower(), latin)
        raw = open(HOMER, encoding="utf-8").read()
        base = "".join(c for c in unicodedata.normalize("NFD", raw) if unicodedata.category(c) != "Mn").lower().replace("ς", "σ")
        greek = Alphabet(name="Greek letters", symbols=tuple("αβγδεζηθικλμνξοπρστυφχψω"))
        hom = tokenize(base, greek)
        extra = {}
        for name in ("genesis", "constitution"):
            path = TEXTS / f"{name}.txt"
            if path.is_file():
                extra[name] = tokenize(open(path, encoding="utf-8").read().lower(), latin)
    idx = lambda tok: np.array([tok.alphabet.chars.index(c) for c in tok.text], dtype=np.int64)
    out = {"antony": (latin.size, idx(ant)), "homer": (greek.size, idx(hom))}
    out.update({k: (latin.size, idx(v)) for k, v in extra.items()})
    return out


# ----------------------------------------------------------------------
# models
# ----------------------------------------------------------------------
class OrderModel:
    """Order-r Markov background with a flat pseudocount, as MarkovBackground (short prefix: marginal)."""

    def __init__(self, units, K, r, pc=1.0):
        self.K, self.r, self.pc = K, r, pc
        n = len(units)
        self.marg = (np.bincount(units, minlength=K) + pc) / (n + pc * K)
        if r:
            ctx = np.zeros(len(units) - r, dtype=np.int64)
            for j in range(r):
                ctx = ctx * K + units[j:len(units) - r + j]
            nxt = units[r:]
            self.T = np.zeros((K ** r, K))
            np.add.at(self.T, (ctx, nxt), 1)
            self.cc = self.T.sum(axis=1)

    def batch(self, W):
        P = self.marg[W[:, 0]].copy()
        for i in range(1, W.shape[1]):
            if i < self.r or self.r == 0:
                P *= self.marg[W[:, i]]
            else:
                ctx = np.zeros(len(W), dtype=np.int64)
                for j in range(self.r):
                    ctx = ctx * self.K + W[:, i - self.r + j]
                P *= (self.T[ctx, W[:, i]] + self.pc) / (self.cc[ctx] + self.pc * self.K)
        return P


class Truth:
    """Order-m Markov chain fitted to a corpus, with exact word probabilities."""

    def __init__(self, units, K, m, pc=0.5):
        self.K, self.m = K, m
        if m == 0:
            c = np.bincount(units, minlength=K) + pc
            self.P = (c / c.sum())[None, :]
        else:
            ctx = np.zeros(len(units) - m, dtype=np.int64)
            for j in range(m):
                ctx = ctx * K + units[j:len(units) - m + j]
            T = np.zeros((K ** m, K))
            np.add.at(T, (ctx, units[m:]), 1)
            T += pc
            self.P = T / T.sum(axis=1, keepdims=True)
        self.cum = np.cumsum(self.P, axis=1)
        # stationary distribution of the m-gram process (m <= 2)
        if m == 0:
            self.pi = self.P[0]
        else:
            S = K ** m
            M = np.zeros((S, S))
            for s in range(S):
                nxt = (s % (K ** (m - 1))) * K + np.arange(K) if m > 1 else np.arange(K)
                M[s, nxt] += self.P[s]
            pi = np.ones(S) / S
            for _ in range(2000):
                new = pi @ M
                if np.abs(new - pi).max() < 1e-13:
                    break
                pi = new
            self.pi = pi / pi.sum()

    def sample(self, N, rng):
        u = rng.random(N)
        out = np.empty(N, dtype=np.int64)
        K, m = self.K, self.m
        if m == 0:
            return np.searchsorted(self.cum[0], u).clip(0, K - 1)
        s = rng.choice(len(self.pi), p=self.pi)
        for j in range(m):
            out[m - 1 - j] = s % K
            s //= K
        ctx = 0
        for j in range(m):
            ctx = ctx * K + out[j]
        mod = K ** (m - 1) if m > 1 else 1
        for i in range(m, N):
            x = min(int(np.searchsorted(self.cum[ctx], u[i])), K - 1)
            out[i] = x
            ctx = (ctx % mod) * K + x if m > 1 else x
        return out

    def batch(self, W):
        K, m = self.K, self.m
        k = W.shape[1]
        if m == 0:
            return np.prod(self.pi[W], axis=1)
        pre = np.zeros(len(W), dtype=np.int64)
        for j in range(m):
            pre = pre * K + W[:, j]
        P = self.pi[pre].copy()
        for i in range(m, k):
            ctx = np.zeros(len(W), dtype=np.int64)
            for j in range(m):
                ctx = ctx * K + W[:, i - m + j]
            P *= self.P[ctx, W[:, i]]
        return P


def kmer_table(units, K, ks=KS):
    out = []
    for k in ks:
        if len(units) < k:
            continue
        win = np.lib.stride_tricks.sliding_window_view(np.asarray(units, dtype=np.int64), k)
        key = np.zeros(len(win), dtype=np.int64)
        for i in range(k):
            key = key * K + win[:, i]
        uk, cnt = np.unique(key, return_counts=True)
        W = np.empty((len(uk), k), dtype=np.int64)
        rem = uk.copy()
        for i in range(k - 1, -1, -1):
            W[:, i] = rem % K
            rem //= K
        out.append((k, W, cnt))
    return out


def score(units, K, model, tab=None):
    """Strict and loose calls, family = every possible word of each tested length."""
    tab = tab or kmer_table(units, K)
    n = len(units)
    cnt = np.concatenate([c for _, _, c in tab])
    p = np.empty(len(cnt))
    en = np.empty(len(cnt))
    pos = 0
    for k, W, c in tab:
        exp = max(0, n - k + 1) * model.batch(W)
        with np.errstate(divide="ignore", invalid="ignore"):
            p[pos:pos + len(c)] = np.where(exp > 0, poisson.sf(c - 1, exp), 0.0)
            en[pos:pos + len(c)] = np.where(exp > 0, c / exp, np.inf)
        pos += len(c)
    q = H.bh_vec(p, m=float(sum(K ** k for k in KS)))
    rep = cnt >= 2
    return {"strict": rep & (q <= 0.05) & (en >= 10.0), "loose": rep & (q <= 0.5) & (en >= 3.0)}, tab


def chosen_order(K, N, c=5):
    r = 0
    while K ** (r + 2) <= N / c and r + 1 <= max(ORDERS):
        r += 1
    return r


# ----------------------------------------------------------------------
def _init():
    _G["corp"] = load_corpora()
    _G["truth"] = {}
    for name, (K, u) in _G["corp"].items():
        for m in TRUE_ORDERS:
            t = Truth(u, K, m)
            big = t.sample(40000, np.random.default_rng([99, m, K]))
            _G["truth"][(name, m)] = (t, big)


def null_task(args):
    name, m, N, rep = args
    K, _ = _G["corp"][name]
    truth, big = _G["truth"][(name, m)]
    rng = np.random.default_rng([5, m, N, rep, K])
    units = truth.sample(N, rng)
    tab = kmer_table(units, K)
    out = {}
    for r in ORDERS:
        model = OrderModel(units, K, r)
        res, _ = score(units, K, model, tab)
        # quality on 3-letter words from the truth, and bias on repeated words
        W3 = np.lib.stride_tricks.sliding_window_view(big, 3)
        pt, pb = truth.batch(W3), model.batch(W3)
        kl = float(np.log(pt / pb).mean())
        t5 = float((pb < pt / 5).mean())
        lr = u2 = nn = 0
        for k, W, c in tab:
            if k not in (3, 4):
                continue
            mk = c >= 2
            if mk.any():
                a, b = model.batch(W[mk]), truth.batch(W[mk])
                ok = b > 0
                d = np.log(a[ok] / b[ok])
                lr += d.sum()
                u2 += (d < -np.log(2)).sum()
                nn += ok.sum()
        out[r] = {"strict": int(res["strict"].sum()), "loose": int(res["loose"].sum()),
                  "kl": kl, "t5": t5, "lr": lr, "u2": u2, "n": nn}
    return ("null", (name, m, N), out)


def power_task(args):
    name, m, N, Lp, copies, rep = args
    K, _ = _G["corp"][name]
    truth, big = _G["truth"][(name, m)]
    rng = np.random.default_rng([6, m, N, Lp * 10 + copies, rep, K])
    units = truth.sample(N, rng)
    s = int(rng.integers(0, len(big) - Lp))
    phrase = big[s:s + Lp].copy()
    for _ in range(1000):
        starts = np.sort(rng.choice(N - Lp + 1, size=copies, replace=False))
        if all(b - a >= Lp for a, b in zip(starts, starts[1:])):
            for st in starts:
                units[st:st + Lp] = phrase
            break
    tab = kmer_table(units, K)
    out = {}
    for r in ORDERS:
        res, _ = score(units, K, OrderModel(units, K, r), tab)
        pos, hit = 0, False
        for k, W, c in tab:
            if k == Lp:
                row = np.where((W == phrase).all(axis=1))[0]
                hit = bool(len(row) and res["strict"][pos + row[0]])
            pos += len(c)
        out[r] = hit
    return ("power", (name, m, N, Lp, copies), out)


def real_task(name):
    K, u = _G["corp"][name]
    out = {}
    for r in ORDERS:
        res, tab = score(u, K, OrderModel(u, K, r))
        out[r] = (int(res["strict"].sum()), int(res["loose"].sum()))
    return ("real", name, out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--power-reps", type=int, default=100)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    if a.quick:
        a.reps, a.power_reps = 4, 3
    corpora = list(load_corpora())
    print("corpora:", {c: len(load_corpora()[c][1]) for c in corpora}, flush=True)
    tasks = [(null_task, (c, m, N, r)) for c in corpora for m in TRUE_ORDERS for N in NS for r in range(a.reps)]
    # phrases that repeat the way words do: short, with several copies
    tasks += [(power_task, (c, m, N, Lp, cp, r)) for c in corpora for m in (0, 1, 2) for N in NS
              for Lp, cp in ((3, 6), (3, 12), (4, 4), (4, 8), (6, 2), (6, 3)) for r in range(a.power_reps)]
    tasks += [(real_task, c) for c in corpora]
    t0 = time.time()
    results = []
    with mp.Pool(a.procs, initializer=_init) as pool:
        jobs = [pool.apply_async(fn, (arg,)) for fn, arg in tasks]
        for i, j in enumerate(jobs):
            results.append(j.get())
            if (i + 1) % 500 == 0:
                print(f"  {i + 1}/{len(jobs)} tasks, {time.time() - t0:.0f}s", flush=True)
    print(f"done in {time.time() - t0:.0f}s")
    summarize(results, a)


def summarize(results, a):
    out = HERE / "results"
    tag = "quick_" if a.quick else ""
    agg = defaultdict(lambda: defaultdict(float))
    for r in results:
        if r[0] == "null":
            _, key, res = r
            for o, d in res.items():
                c = agg[key + (o,)]
                c["n"] += 1
                c["strict"] += d["strict"] > 0
                c["loose"] += d["loose"] > 0
                c["kl"] += d["kl"]
                c["t5"] += d["t5"]
                c["lr"] += d["lr"]
                c["u2"] += d["u2"]
                c["nn"] += d["n"]
    with open(out / f"{tag}letters_null.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["corpus", "true_order", "N", "background_order", "texts", "fp_strict", "fp_loose", "kl_per_word",
                    "mass_underestimated_5x", "repeated_mean_log_ratio", "repeated_frac_underestimated_2x"])
        for (name, m, N, o), c in sorted(agg.items()):
            nn = max(c["nn"], 1)
            w.writerow([name, m, N, o, int(c["n"]), round(c["strict"] / c["n"], 4), round(c["loose"] / c["n"], 4),
                        round(c["kl"] / c["n"], 4), round(c["t5"] / c["n"], 4), round(c["lr"] / nn, 4), round(c["u2"] / nn, 4)])
    pw = defaultdict(lambda: defaultdict(float))
    for r in results:
        if r[0] == "power":
            _, key, res = r
            for o, hit in res.items():
                c = pw[key + (o,)]
                c["n"] += 1
                c["hit"] += hit
    with open(out / f"{tag}letters_power.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["corpus", "true_order", "N", "phrase_letters", "copies", "background_order", "texts", "recall"])
        for (name, m, N, Lp, cp, o), c in sorted(pw.items()):
            w.writerow([name, m, N, Lp, cp, o, int(c["n"]), round(c["hit"] / c["n"], 4)])
    with open(out / f"{tag}letters_real.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["corpus", "background_order", "supported_words_strict", "supported_words_loose"])
        for r in results:
            if r[0] == "real":
                for o, (s, l) in sorted(r[2].items()):
                    w.writerow([r[1], o, s, l])
    print(f"wrote {out}/{tag}letters_*.csv")


if __name__ == "__main__":
    main()

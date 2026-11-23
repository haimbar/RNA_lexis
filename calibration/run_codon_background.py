"""Calibration of the codon-mode background (scoping doc, section 7).

Question: which background model, and which multiple-testing family, keeps false positives under
control on codon text while still recovering planted repeats?

Null texts (no planted phrases), N codons:
  A  independent codons drawn from the usage of a real human gene (the model an order-0
     codon-usage background assumes);
  B  amino acids from a first-order Markov chain fitted to all genes, codons chosen from the
     same gene's conditional usage P(codon | amino acid)  (structure beyond order 0).
Planted-repeat texts: null A plus r exact copies of a phrase of Lp codons.

Experiments
  scoring  every repeated word of 2..6 codons is scored under each background and three
           multiple-testing families (see harness.scoring_stage).  This isolates the
           background and the family from the discovery step.
  pipeline the CLI's full xmotif -> shared-substring -> score pipeline on the same null texts
           (shows how often the discovery step produces anything to test at all).
  nt       the package's own pipeline on iid nucleotide text, as a reference for what the
           existing defaults do.
  real     the scoring stage on real genes (descriptive: real genes contain real repeats).

Usage:  python3 run_codon_background.py [--reps 400] [--power-reps 200] [--quick]
Writes calibration/results/*.csv.
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H

HERE = Path(__file__).resolve().parent
SETTINGS = [(2, 2), (3, 2), (3, 3), (4, 3), (4, 4)]        # (xmotif minlen, candidate minlen) in codons
NS = [100, 300, 1500]
MIN_GENE = 250                                              # genes used as generating models
FAMILIES = H.FAMILIES

_G = {}                                                     # per-process data


def _init():
    genes = H.load_genes(HERE / "data" / "human_cds.fasta")
    pooled = sum(np.bincount(u, minlength=H.NCOD) for u in genes.values()).astype(float)
    T = np.zeros((len(H.AAS), len(H.AAS)))
    for u in genes.values():
        a = H.AA_NUM[u]
        ok = (a[:-1] >= 0) & (a[1:] >= 0)
        np.add.at(T, (a[:-1][ok], a[1:][ok]), 1)
    T = (T + 0.5) / (T + 0.5).sum(axis=1, keepdims=True)
    pi = np.ones(len(H.AAS)) / len(H.AAS)
    for _ in range(500):
        pi = pi @ T
    cond_pool = np.zeros(H.NCOD)
    for aa in range(len(H.AAS)):
        m = H.AA_NUM == aa
        cond_pool[m] = pooled[m] / pooled[m].sum()
    _G.update(genes=genes, pooled=pooled, T=T, pi=pi, cond_pool=cond_pool,
              names=sorted(s for s, u in genes.items() if len(u) >= MIN_GENE),
              all_names=sorted(genes))


def gene_models(sym):
    u = _G["genes"][sym]
    n = np.bincount(u, minlength=H.NCOD).astype(float)
    usage = n / n.sum()
    cond = np.zeros(H.NCOD)
    for aa in range(len(H.AAS)):
        m = H.AA_NUM == aa
        cond[m] = (n[m] + 5 * _G["cond_pool"][m]) / (n[m].sum() + 5)
    return n, usage, cond


def sample(kind, usage, cond, N, rng):
    if kind == "A":
        return rng.choice(H.NCOD, size=N, p=usage).astype(np.int16)
    T, pi = _G["T"], _G["pi"]
    aa = np.empty(N, dtype=int)
    aa[0] = rng.choice(len(H.AAS), p=pi)
    for i in range(1, N):
        aa[i] = rng.choice(len(H.AAS), p=T[aa[i - 1]])
    out = np.empty(N, dtype=np.int16)
    for a in range(len(H.AAS)):
        idx = np.where(aa == a)[0]
        if len(idx):
            cods = np.where(H.AA_NUM == a)[0]
            out[idx] = rng.choice(cods, size=len(idx), p=cond[cods] / cond[cods].sum())
    return out


def truth_model(kind, usage, cond):
    if kind == "A":
        return H.Model(lambda W: np.prod(usage[W], axis=1))
    T, pi = _G["T"], _G["pi"]

    def batch(W):
        a = H.AA_NUM[W]
        p = pi[np.maximum(a[:, 0], 0)] * np.prod(cond[W], axis=1)
        for i in range(1, W.shape[1]):
            p = p * T[np.maximum(a[:, i - 1], 0), np.maximum(a[:, i], 0)]
        return np.where((a < 0).any(axis=1), 0.0, p)
    return H.Model(batch)


def plant(units, phrase, copies, rng):
    N, L = len(units), len(phrase)
    for _ in range(1000):
        starts = np.sort(rng.choice(N - L + 1, size=copies, replace=False))
        if all(b - a >= L for a, b in zip(starts, starts[1:])):
            out = units.copy()
            for s in starts:
                out[s:s + L] = phrase
            return out
    raise RuntimeError("could not place the copies")


def backgrounds(sym, usage, cond, kind, with_oracle=True):
    n_g, _, _ = gene_models(sym)
    ref = _G["pooled"] - n_g + 0.5                          # leave-one-out pooled usage
    bgs = [H.bg_order0_flat(1.0), H.bg_order0_flat(0.1), H.bg_ntpos(),
           H.bg_hybrid(16), H.bg_hybrid(64), H.bg_hybrid(256), H.bg_hybrid(1024), H.bg_hybrid(4096),
           H.bg_hybrid(4), H.bg_order0_flat(0.5), H.bg_order0_flat(2.0),
           H.bg_masked_flat(1.0), H.bg_masked_hybrid(64), H.bg_masked_hybrid(256), H.bg_masked_hybrid(1024),
           H.bg_order1_flat(1.0), H.bg_order1_backoff(64), H.bg_reference(ref)]
    if with_oracle:
        bgs.append(H.bg_oracle(truth_model(kind, usage, cond)))
    return bgs


# ----------------------------------------------------------------------
# tasks
# ----------------------------------------------------------------------
def scoring_null_task(args):
    kind, N, rep = args
    rng = np.random.default_rng([{"A": 1, "B": 2}[kind], N, rep])
    sym = _G["names"][rng.integers(len(_G["names"]))]
    _, usage, cond = gene_models(sym)
    units = sample(kind, usage, cond, N, rng)
    out = {}
    for b in backgrounds(sym, usage, cond, kind):
        res, _, cnt, _ = H.scoring_stage(units, b(units.astype(int)))
        out[b.name] = {f: int(res[f].sum()) for f in FAMILIES}
    return ("snull", (kind, N), out)


def scoring_power_task(args):
    N, Lp, r, rep = args
    rng = np.random.default_rng([3, N, Lp * 10 + r, rep])
    sym = _G["names"][rng.integers(len(_G["names"]))]
    _, usage, cond = gene_models(sym)
    units = sample("A", usage, cond, N, rng)
    phrase = sample("A", usage, cond, Lp, rng)
    units = plant(units, phrase, r, rng)
    out = {}
    for b in backgrounds(sym, usage, cond, "A"):
        res, ks_, cnt, words = H.scoring_stage(units, b(units.astype(int)))
        pos = 0
        hit = {}
        for k, W, c in words:
            if k == Lp:
                row = np.where((W == phrase.astype(np.int64)).all(axis=1))[0]
                hit = {f: bool(res[f][pos + row[0]]) for f in FAMILIES} if len(row) else {f: False for f in FAMILIES}
            pos += len(c)
        out[b.name] = hit
    return ("spower", (N, Lp, r), out)


def quality_task(args):
    """How wrong are a background's word probabilities?  Compared with the truth on words drawn from it."""
    kind, N, rep = args
    rng = np.random.default_rng([{"A": 21, "B": 22}[kind], N, rep])
    sym = _G["names"][rng.integers(len(_G["names"]))]
    _, usage, cond = gene_models(sym)
    units = sample(kind, usage, cond, N, rng)
    truth = truth_model(kind, usage, cond)
    big = sample(kind, usage, cond, 20000, rng)                 # words drawn from the true model
    out = {}
    for b in backgrounds(sym, usage, cond, kind, with_oracle=False):
        model = b(units.astype(int))
        if model.uses_counts:
            continue
        row = {}
        for k in (2, 3, 4):
            W = np.lib.stride_tricks.sliding_window_view(big.astype(np.int64), k)
            pt, pb = truth.batch(W), model.batch(W)
            with np.errstate(divide="ignore"):
                lr = np.log(pt / pb)
            row[k] = (float(lr.mean()), float((pb < pt / 5).mean()), float((pb < pt / 2).mean()))
        out[b.name] = row
    return ("quality", (kind, N), out)


def repeated_words_task(args):
    """Bias of the expected count for the words that can actually be called significant:
    words seen >= 2 times in a null text, compared with their true expected count."""
    kind, N, rep = args
    rng = np.random.default_rng([{"A": 31, "B": 32}[kind], N, rep])
    sym = _G["names"][rng.integers(len(_G["names"]))]
    _, usage, cond = gene_models(sym)
    units = sample(kind, usage, cond, N, rng)
    truth = truth_model(kind, usage, cond)
    out = {}
    tab = H.kmer_table(units, (2, 3, 4))
    for b in backgrounds(sym, usage, cond, kind, with_oracle=False):
        model = b(units.astype(int))
        lr, under2, under5, n = 0.0, 0, 0, 0
        for k, W, c in tab:
            m = c >= 2
            if not m.any():
                continue
            Wm, cm = W[m], c[m]
            pb = model.batch(Wm, cm) if model.uses_counts else model.batch(Wm)
            pt = truth.batch(Wm)
            ok = pt > 0
            r = np.log(pb[ok] / pt[ok])
            lr += r.sum()
            under2 += (r < -np.log(2)).sum()
            under5 += (r < -np.log(5)).sum()
            n += ok.sum()
        out[b.name] = (lr, under2, under5, n)
    return ("repeated", (kind, N), out)


def pipeline_null_task(args):
    kind, N, rep = args
    rng = np.random.default_rng([{"A": 11, "B": 12}[kind], N, rep])
    sym = _G["names"][rng.integers(len(_G["names"]))]
    _, usage, cond = gene_models(sym)
    units = sample(kind, usage, cond, N, rng)
    text = H.to_text(units)
    probs = {b.name: b(units.astype(int)) for b in backgrounds(sym, usage, cond, kind)}
    res = {}
    for mx, cm in SETTINGS:
        xm, cands = H.discover(text, mx, cm)
        row = {"n_cand": len(cands)}
        for name, model in probs.items():
            rows = H.score(text, cands, xm, model) if cands else []
            row[name] = sum(r["supported"] for r in rows)
        res[(mx, cm)] = row
    return ("pipeline", (kind, N), res)


def nt_task(args):
    """The package's own pipeline on iid nucleotide text (the existing defaults)."""
    import contextlib
    import io
    from rna_lexis.statistical import rank_core_candidates
    from rna_lexis.algorithms import find_boundary
    n_nt, rep = args
    rng = np.random.default_rng([4, n_nt, rep])
    sym = _G["names"][rng.integers(len(_G["names"]))]
    cod = np.bincount(_G["genes"][sym], minlength=H.NCOD)
    nt = np.zeros(4)
    for c in range(H.NCOD):
        for j in range(3):
            nt[H.NT_OF[c, j]] += cod[c]
    nt /= nt.sum()
    txt = "".join("acgt"[i] for i in rng.choice(4, size=n_nt, p=nt))
    out = {}
    for label, (mlen, at_least, cmin) in {"cli_default(7,2,5)": (7, 2, 5), "paper(8,3,7)": (8, 3, 7)}.items():
        with contextlib.redirect_stdout(io.StringIO()):
            xm = sorted(set(find_boundary(txt, minlen=mlen, maxlen=60, at_least=at_least)), key=lambda x: (-len(x), x))
            rows = rank_core_candidates(txt, xm, candidate_min_len=cmin, markov_order=1) if xm else []
        out[label] = {"n_cand": len(rows), "supported": sum(r["statistically_supported"] for r in rows)}
    return ("nt", n_nt, out)


def real_task(sym):
    units = _G["genes"][sym][:1500].astype(np.int64)
    if len(units) < 100:
        return ("real", sym, None)
    _, usage, cond = gene_models(sym)
    out = {}
    for b in backgrounds(sym, usage, cond, "A", with_oracle=False):
        res, _, _, _ = H.scoring_stage(units, b(units))
        out[b.name] = {f: int(res[f].sum()) for f in FAMILIES}
    return ("real", sym, out)


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=400)
    ap.add_argument("--power-reps", type=int, default=200)
    ap.add_argument("--pipeline-reps", type=int, default=150)
    ap.add_argument("--nt-reps", type=int, default=300)
    ap.add_argument("--quick", action="store_true", help="tiny run to check the plumbing")
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    if a.quick:
        a.reps, a.power_reps, a.pipeline_reps, a.nt_reps = 6, 4, 4, 6
    print(f"selfcheck: {H.selfcheck()} candidates match the package", flush=True)

    _init()
    real_syms = [s for s in _G["all_names"] if len(_G["genes"][s]) >= 100]
    tasks = [(scoring_null_task, (k, N, r)) for k in "AB" for N in NS for r in range(a.reps)]
    tasks += [(scoring_power_task, (N, Lp, r, rep)) for N in NS for Lp in (3, 4, 6) for r in (2, 3, 4)
              for rep in range(a.power_reps)]
    tasks += [(pipeline_null_task, (k, N, r)) for k in "AB" for N in NS for r in range(a.pipeline_reps)]
    tasks += [(quality_task, (k, N, r)) for k in "AB" for N in NS for r in range(a.reps // 2)]
    tasks += [(repeated_words_task, (k, N, r)) for k in "AB" for N in NS for r in range(a.reps)]
    tasks += [(nt_task, (n, rep)) for n in (900, 4500) for rep in range(a.nt_reps)]
    tasks += [(real_task, s) for s in real_syms]

    t0 = time.time()
    results = []
    with mp.Pool(a.procs, initializer=_init) as pool:
        jobs = [pool.apply_async(fn, (arg,)) for fn, arg in tasks]
        for i, j in enumerate(jobs):
            results.append(j.get())
            if (i + 1) % 1000 == 0:
                print(f"  {i + 1}/{len(jobs)} tasks, {time.time() - t0:.0f}s", flush=True)
    print(f"done in {time.time() - t0:.0f}s")
    summarize(results, a)


def summarize(results, a):
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    tag = "quick_" if a.quick else ""

    # ---- scoring: null ----
    agg = defaultdict(lambda: defaultdict(float))
    for r in results:
        if r[0] == "snull":
            _, (kind, N), res = r
            for bg, d in res.items():
                for f in FAMILIES:
                    c = agg[(kind, N, bg, f)]
                    c["texts"] += 1
                    c["fp"] += d[f] > 0
                    c["nsup"] += d[f]
    with open(out / f"{tag}scoring_null.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["null", "N", "background", "family", "texts", "fp_rate", "mean_supported"])
        for (kind, N, bg, f), d in sorted(agg.items()):
            w.writerow([kind, N, bg, f, int(d["texts"]), round(d["fp"] / d["texts"], 4), round(d["nsup"] / d["texts"], 4)])

    # ---- scoring: power ----
    pw = defaultdict(lambda: defaultdict(float))
    for r in results:
        if r[0] == "spower":
            _, key, res = r
            for bg, d in res.items():
                for f in FAMILIES:
                    c = pw[(key, bg, f)]
                    c["texts"] += 1
                    c["rec"] += d[f]
    with open(out / f"{tag}scoring_power.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["N", "phrase_codons", "copies", "background", "family", "texts", "recall"])
        for ((N, Lp, rr), bg, f), d in sorted(pw.items()):
            w.writerow([N, Lp, rr, bg, f, int(d["texts"]), round(d["rec"] / d["texts"], 4)])

    # ---- background quality ----
    ql = defaultdict(lambda: defaultdict(float))
    for r in results:
        if r[0] == "quality":
            _, (kind, N), res = r
            for bg, rows in res.items():
                for k, (kl, t5, t2) in rows.items():
                    c = ql[(kind, N, bg, k)]
                    c["n"] += 1
                    c["kl"] += kl
                    c["t5"] += t5
                    c["t2"] += t2
    with open(out / f"{tag}background_quality.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["null", "N", "background", "word_length", "texts", "kl_per_word", "mass_underestimated_5x", "mass_underestimated_2x"])
        for (kind, N, bg, k), c in sorted(ql.items()):
            w.writerow([kind, N, bg, k, int(c["n"]), round(c["kl"] / c["n"], 4), round(c["t5"] / c["n"], 4), round(c["t2"] / c["n"], 4)])

    # ---- bias on repeated words ----
    rp = defaultdict(lambda: defaultdict(float))
    for r in results:
        if r[0] == "repeated":
            _, (kind, N), res = r
            for bg, (lr, u2, u5, n) in res.items():
                c = rp[(kind, N, bg)]
                c["lr"] += lr
                c["u2"] += u2
                c["u5"] += u5
                c["n"] += n
    with open(out / f"{tag}repeated_word_bias.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["null", "N", "background", "repeated_words", "mean_log_expected_ratio", "frac_underestimated_2x", "frac_underestimated_5x"])
        for (kind, N, bg), c in sorted(rp.items()):
            n = max(c["n"], 1)
            w.writerow([kind, N, bg, int(c["n"]), round(c["lr"] / n, 4), round(c["u2"] / n, 4), round(c["u5"] / n, 4)])

    # ---- pipeline null ----
    pl = defaultdict(lambda: defaultdict(float))
    for r in results:
        if r[0] == "pipeline":
            _, (kind, N), res = r
            for (mx, cm), row in res.items():
                c = pl[(kind, N, mx, cm)]
                c["texts"] += 1
                c["with_cand"] += row["n_cand"] > 0
                c["cand"] += row["n_cand"]
                for bg, v in row.items():
                    if bg != "n_cand":
                        c[bg] += v > 0
    with open(out / f"{tag}pipeline_null.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["null", "N", "xmotif_minlen", "cand_minlen", "texts", "frac_with_candidates",
                    "mean_candidates", "fp_rate_by_background"])
        for (kind, N, mx, cm), c in sorted(pl.items()):
            bgs = {k: round(v / c["texts"], 4) for k, v in c.items() if k not in ("texts", "with_cand", "cand")}
            w.writerow([kind, N, mx, cm, int(c["texts"]), round(c["with_cand"] / c["texts"], 4),
                        round(c["cand"] / c["texts"], 3), ";".join(f"{k}={v}" for k, v in bgs.items())])

    # ---- nucleotide reference ----
    nt = defaultdict(lambda: defaultdict(float))
    for r in results:
        if r[0] == "nt":
            _, n_nt, res = r
            for label, d in res.items():
                c = nt[(n_nt, label)]
                c["texts"] += 1
                c["with_cand"] += d["n_cand"] > 0
                c["fp"] += d["supported"] > 0
    with open(out / f"{tag}nucleotide_reference.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["n_nt", "settings(xmotif_minlen,at_least,cand_minlen)", "texts", "frac_with_candidates", "fp_rate"])
        for (n_nt, label), c in sorted(nt.items()):
            w.writerow([n_nt, label, int(c["texts"]), round(c["with_cand"] / c["texts"], 4), round(c["fp"] / c["texts"], 4)])

    # ---- real genes ----
    rl = defaultdict(lambda: defaultdict(float))
    for r in results:
        if r[0] == "real" and r[2] is not None:
            for bg, d in r[2].items():
                for f in FAMILIES:
                    c = rl[(bg, f)]
                    c["genes"] += 1
                    c["any"] += d[f] > 0
                    c["nsup"] += d[f]
    with open(out / f"{tag}real_genes.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["background", "family", "genes", "frac_genes_with_supported_word", "mean_supported_words"])
        for (bg, f), c in sorted(rl.items()):
            w.writerow([bg, f, int(c["genes"]), round(c["any"] / c["genes"], 4), round(c["nsup"] / c["genes"], 3)])
    print(f"wrote {out}/{tag}*.csv")


if __name__ == "__main__":
    main()

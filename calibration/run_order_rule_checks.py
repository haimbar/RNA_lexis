"""Checks of the Markov-order rule where the earlier calibration did not reach.

The rule (statistical.default_markov_order): order 1 when the text has at least K**2 units, else order 0.
Earlier simulations covered codons up to 1,500 codons (K = 64, so K**2 = 4,096) and letters (K = 26).  This script
adds (a) codon text of 4,000 to 20,000 codons, where the rule switches to order 1, and (b) amino acids
(K = 21 declared, K**2 = 441) from 100 to 8,000 residues.

For each text, backgrounds of order 0 and 1 (flat pseudocount 1) are compared on the words seen at least twice:
mean log(expected / true expected) (negative = anti-conservative) and the false-positive rate of the strict rule with
the per-length correction.

Usage: python3 run_order_rule_checks.py [--reps 60]
Writes results/order_rule_checks.csv.
"""
import argparse
import csv
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))
sys.path.insert(0, str(HERE))
import harness as H                                              # noqa: E402
import run_codon_background as CB                                # noqa: E402
import run_letters_background as LB                              # noqa: E402

AA_K = 20
_S = {}


def _init():
    CB._init()
    genes = CB._G["genes"]
    res = [H.AA_NUM[u][H.AA_NUM[u] >= 0] for u in genes.values()]
    protein = np.concatenate(res).astype(np.int64)
    _S["protein"] = protein
    _S["truth"] = {m: LB.Truth(protein, AA_K, m) for m in (0, 1)}


def _bylength_calls(units, K, declared, model):
    """Strict calls with the per-length correction over declared**L possible words (as score_exact_motifs)."""
    tab = LB.kmer_table(units, K, ks=(2, 3, 4, 5))
    n = len(units)
    from scipy.stats import poisson
    hits = 0
    for k, W, c in tab:
        exp = max(0, n - k + 1) * model.batch(W)
        with np.errstate(divide="ignore", invalid="ignore"):
            p = np.where(exp > 0, poisson.sf(c - 1, exp), 0.0)
            en = np.where(exp > 0, c / exp, np.inf)
        rep = c >= 2
        if not rep.any():
            continue
        q = H.bh_vec(p[rep] if False else p, m=float(declared) ** k)
        hits += int(((q <= 0.05) & (en >= 10.0) & rep).sum())
    return hits


def task(args):
    kind, N, rep = args
    if kind.startswith("codon"):
        rng = np.random.default_rng([41, {"codonA": 1, "codonB": 2}[kind], N, rep])
        sym = CB._G["names"][rng.integers(len(CB._G["names"]))]
        _, usage, cond = CB.gene_models(sym)
        units = CB.sample(kind[-1], usage, cond, N, rng).astype(np.int64)
        K, declared = 64, 64
        truth = CB.truth_model(kind[-1], usage, cond)
        models = {0: H.bg_order0_flat(1.0)(units), 1: H.bg_order1_flat(1.0)(units)}
        pt = truth.batch
    else:
        m = 0 if kind == "aminoA" else 1
        rng = np.random.default_rng([42, m, N, rep])
        truth = _S["truth"][m]
        units = truth.sample(N, rng)
        K, declared = AA_K, 21
        models = {0: LB.OrderModel(units, K, 0), 1: LB.OrderModel(units, K, 1)}
        pt = truth.batch
    out = {}
    tab = LB.kmer_table(units, K, ks=(3, 4))
    for order, model in models.items():
        lr = u2 = nn = 0
        for k, W, c in tab:
            mk = c >= 2
            if mk.any():
                a, b = model.batch(W[mk]), pt(W[mk])
                ok = b > 0
                d = np.log(a[ok] / b[ok])
                lr += d.sum()
                u2 += (d < -np.log(2)).sum()
                nn += ok.sum()
        out[order] = (lr, u2, nn, _bylength_calls(units, K, declared, model))
    return kind, N, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    sizes = {"codonA": (4000, 8000, 20000), "codonB": (4000, 8000, 20000),
             "aminoA": (100, 300, 1000, 3000, 8000), "aminoB": (100, 300, 1000, 3000, 8000)}
    tasks = [(k, N, r) for k, ns in sizes.items() for N in ns for r in range(a.reps)]
    t0 = time.time()
    with mp.Pool(a.procs, initializer=_init) as pool:
        results = pool.map(task, tasks, chunksize=4)
    print(f"done in {time.time() - t0:.0f}s")
    agg = defaultdict(lambda: defaultdict(float))
    for kind, N, out in results:
        for order, (lr, u2, nn, hits) in out.items():
            c = agg[(kind, N, order)]
            c["n"] += 1
            c["lr"] += lr
            c["u2"] += u2
            c["nn"] += nn
            c["fp"] += hits > 0
    path = HERE / "results" / "order_rule_checks.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["null", "N", "background_order", "texts", "repeated_words", "mean_log_expected_ratio",
                    "frac_underestimated_2x", "fp_rate"])
        for (kind, N, order), c in sorted(agg.items()):
            nn = max(c["nn"], 1)
            w.writerow([kind, N, order, int(c["n"]), int(c["nn"]), round(c["lr"] / nn, 3),
                        round(c["u2"] / nn, 3), round(c["fp"] / c["n"], 3)])
    print(open(path, encoding="utf-8").read())


if __name__ == "__main__":
    main()

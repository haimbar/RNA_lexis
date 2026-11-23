"""Calibration of the discovery settings per alphabet (scoping doc, sections 7 and 10, phase 5).

Question: for codons, Latin letters and amino acids, which minimum xmotif length, minimum core
(candidate) length and minimum occurrence count keep the false-positive rate of the whole pipeline
at or below 5% per null text, while still recovering planted core families?

The pipeline is the package's own, as the CLI runs it:

    xmotifs = find_boundary(text, minlen=mx, maxlen=60, at_least=mc)
    rows    = rank_core_candidates(text, xmotifs, candidate_min_len=cm, min_xmotif_type_support=2,
                                   min_count=mc, alphabet=A)            # order and family by the calibrated rule

Null texts (nothing planted), N units, from a known generating model:
  codon   independent codons from a real gene's usage (A) or amino-acid Markov + codon choice (B);
  letters an order-1 Markov chain fitted to a real English text (Antony, Genesis, the Constitution);
  amino   independent residues, or a residue Markov chain, fitted to the translated genes.

Planted texts: a core of Lc units with `copies` copies, each followed by one of two different
2-unit extensions (the pipeline needs two distinct xmotifs that share the core).  A planted core
is "recovered" when a statistically supported motif contains at least 3 consecutive units of it.

Usage: python3 run_pipeline_defaults.py [--null-reps 80] [--plant-reps 40] [--quick]
Writes results/pipeline_defaults.csv and results/chosen_defaults.json.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
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
from rna_lexis.algorithms import find_boundary                   # noqa: E402
from rna_lexis.alphabet import codon, latin                      # noqa: E402
from rna_lexis.genetic_code import amino_acid                    # noqa: E402
from rna_lexis.statistical import rank_core_candidates           # noqa: E402

# (xmotif minimum length, candidate minimum length): candidate length cannot exceed xmotif length
PAIRS = [(3, 2), (3, 3), (4, 3), (4, 4), (5, 3), (5, 4), (5, 5), (6, 4), (6, 5), (6, 6)]
MINCOUNTS = (2, 3)
SIZES = {"codon": (100, 300, 1500), "letters": (1000, 3000, 8000), "amino": (100, 300, 1500)}
CORES = {"codon": (3, 4, 6), "letters": (4, 6, 8), "amino": (3, 4, 6)}
COPIES = (4, 6)
ALPHA = {"codon": codon, "letters": latin, "amino": amino_acid}
_S = {}


def _init():
    CB._init()
    LB_corp = LB.load_corpora()
    _S["corp"] = {k: v for k, v in LB_corp.items() if k in ("antony", "genesis", "constitution")}
    _S["truth"] = {name: LB.Truth(u, K, 1) for name, (K, u) in _S["corp"].items()}
    _S["alpha"] = {k: f() for k, f in ALPHA.items()}
    # residue chain from the translated genes: P(aa) and the Markov-1 transitions
    T, pi = CB._G["T"], CB._G["pi"]
    _S["aa_T"], _S["aa_pi"] = T, pi
    aas = [a.lower() for a in H.AAS]                              # same order as the chain's rows
    _S["aa_labels"] = aas


# ----------------------------------------------------------------------
# texts
# ----------------------------------------------------------------------
def make_units(kind, N, rng):
    """Integer units for a null text of `kind`: returns (units, chars) with chars[i] = internal char of unit i."""
    a = _S["alpha"]
    if kind.startswith("codon"):
        sym = CB._G["names"][rng.integers(len(CB._G["names"]))]
        _, usage, cond = CB.gene_models(sym)
        units = CB.sample(kind[-1], usage, cond, N, rng).astype(int)
        return units, list(a["codon"].chars)
    if kind.startswith("letters"):
        name = ("antony", "genesis", "constitution")[rng.integers(3)]
        units = _S["truth"][name].sample(N, rng).astype(int)
        return units, list(a["letters"].chars[:26])
    # amino acids: a Markov-1 chain over the 20 residues (kind 'amino' = iid from the stationary frequencies)
    T, pi = _S["aa_T"], _S["aa_pi"]
    labels = _S["aa_labels"]
    chars = [a["amino"].char(x) for x in labels]
    if kind == "amino_iid":
        units = rng.choice(len(labels), size=N, p=pi).astype(int)
    else:
        units = np.empty(N, dtype=int)
        units[0] = rng.choice(len(labels), p=pi)
        for i in range(1, N):
            units[i] = rng.choice(len(labels), p=T[units[i - 1]])
    return units, chars


def plant(units, Lc, copies, rng, K):
    """Plant a core with two alternating 2-unit extensions; return (units, core)."""
    units = units.copy()
    N = len(units)
    core = rng.integers(K, size=Lc)
    ext = [rng.integers(K, size=2), rng.integers(K, size=2)]
    width = Lc + 2
    for _ in range(1000):
        starts = np.sort(rng.choice(N - width, size=copies, replace=False))
        if all(b - a >= width for a, b in zip(starts, starts[1:])):
            for k, s in enumerate(starts):
                units[s:s + Lc] = core
                units[s + Lc:s + width] = ext[k % 2]
            return units, core
    raise RuntimeError("could not place the copies")


def run_pipeline(text, A, mx, cm, mc):
    with contextlib.redirect_stdout(io.StringIO()):
        xm = sorted(set(find_boundary(text, minlen=mx, maxlen=60, at_least=mc)), key=lambda x: (-len(x), x))
        rows = rank_core_candidates(text, xm, candidate_min_len=cm, candidate_max_len=18,
                                    min_xmotif_type_support=2, min_count=mc, alphabet=A) if xm else []
    return [r["motif"] for r in rows if r["statistically_supported"]]


def evaluate(text, A, core_s=None):
    """Run every setting on one text; return {(mx, cm, mc): hit} (any supported / planted core recovered)."""
    out = {}
    for mc in MINCOUNTS:
        for mx in sorted({p[0] for p in PAIRS}):
            with contextlib.redirect_stdout(io.StringIO()):
                xm = sorted(set(find_boundary(text, minlen=mx, maxlen=60, at_least=mc)), key=lambda x: (-len(x), x))
            for (px, cm) in PAIRS:
                if px != mx:
                    continue
                if xm:
                    with contextlib.redirect_stdout(io.StringIO()):
                        rows = rank_core_candidates(text, xm, candidate_min_len=cm, candidate_max_len=18,
                                                    min_xmotif_type_support=2, min_count=mc, alphabet=A)
                    sup = [r["motif"] for r in rows if r["statistically_supported"]]
                else:
                    sup = []
                if core_s is None:
                    out[(mx, cm, mc)] = len(sup) > 0
                else:
                    grams = {core_s[i:i + 3] for i in range(len(core_s) - 2)}
                    out[(mx, cm, mc)] = any(any(g in m for g in grams) for m in sup)
    return out


def _kind_id(kind):
    """A stable number for a text kind (str hashes change from run to run)."""
    return sum(ord(c) * 31 ** i for i, c in enumerate(kind)) % 100000


def null_task(args):
    alph, kind, N, rep = args
    rng = np.random.default_rng([21, _kind_id(kind), N, rep])
    units, chars = make_units(kind, N, rng)
    text = "".join(chars[u] for u in units)
    return ("null", alph, kind, N, evaluate(text, _S["alpha"][alph]))


def plant_task(args):
    alph, kind, N, Lc, copies, rep = args
    rng = np.random.default_rng([22, _kind_id(kind), N, Lc * 10 + copies, rep])
    units, chars = make_units(kind, N, rng)
    units, core = plant(units, Lc, copies, rng, len(chars))
    text = "".join(chars[u] for u in units)
    core_s = "".join(chars[u] for u in core)
    return ("plant", alph, kind, N, (Lc, copies), evaluate(text, _S["alpha"][alph], core_s))


KINDS = {"codon": ("codonA", "codonB"), "letters": ("letters",), "amino": ("amino_iid", "amino_markov")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null-reps", type=int, default=80)
    ap.add_argument("--plant-reps", type=int, default=40)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    if a.quick:
        a.null_reps, a.plant_reps = 3, 2
    tasks = []
    for alph in SIZES:
        for kind in KINDS[alph]:
            for N in SIZES[alph]:
                tasks += [(null_task, (alph, kind, N, r)) for r in range(a.null_reps)]
                for Lc in CORES[alph]:
                    for cp in COPIES:
                        tasks += [(plant_task, (alph, kind, N, Lc, cp, r)) for r in range(a.plant_reps)]
    print(f"{len(tasks)} tasks", flush=True)
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
    null = defaultdict(lambda: defaultdict(float))
    plant_ = defaultdict(lambda: defaultdict(float))
    for r in results:
        if r[0] == "null":
            _, alph, kind, N, res = r
            for key, hit in res.items():
                c = null[(alph, kind, N, key)]
                c["n"] += 1
                c["fp"] += hit
        else:
            _, alph, kind, N, (Lc, cp), res = r
            for key, hit in res.items():
                c = plant_[(alph, kind, N, Lc, cp, key)]
                c["n"] += 1
                c["hit"] += hit
    with open(out / f"{tag}pipeline_defaults.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["alphabet", "kind", "N", "core_units", "copies", "xmotif_minlen", "cand_minlen", "min_count",
                    "texts", "rate"])
        for (alph, kind, N, key), c in sorted(null.items()):
            w.writerow([alph, kind, N, "", "", *key, int(c["n"]), round(c["fp"] / c["n"], 4)])
        for (alph, kind, N, Lc, cp, key), c in sorted(plant_.items()):
            w.writerow([alph, kind, N, Lc, cp, *key, int(c["n"]), round(c["hit"] / c["n"], 4)])
    print(f"wrote {out}/{tag}pipeline_defaults.csv")
    if not a.quick:
        json.dump(choose(null, plant_), open(out / "chosen_defaults.json", "w"), indent=2)
        print("wrote chosen_defaults.json")


def choose(null, plant_, max_fp=0.05):
    """Per alphabet: the setting with the best mean recall among those whose FP rate is at most max_fp at
    every length and in every null (ties: longer minimum lengths, then higher minimum count)."""
    chosen = {}
    for alph in SIZES:
        scores = {}
        keys = {k[3] for k in null if k[0] == alph}
        for key in keys:
            fps = [c["fp"] / c["n"] for (al, kind, N, k), c in null.items() if al == alph and k == key]
            worst = max(fps)
            rec = [c["hit"] / c["n"] for (al, kind, N, Lc, cp, k), c in plant_.items() if al == alph and k == key]
            scores[key] = (worst, float(np.mean(rec)))
        ok = {k: v for k, v in scores.items() if v[0] <= max_fp}
        if not ok:
            chosen[alph] = {"error": "no setting keeps false positives at or below the limit",
                            "worst_fp_by_setting": {str(k): v[0] for k, v in scores.items()}}
            continue
        best = max(ok, key=lambda k: (round(ok[k][1], 3), k[0], k[1], k[2]))
        chosen[alph] = {"xmotif_minlen": best[0], "cand_minlen": best[1], "min_count": best[2],
                        "worst_fp": ok[best][0], "mean_recall": ok[best][1],
                        "all_settings": {str(k): {"worst_fp": v[0], "mean_recall": v[1]} for k, v in sorted(scores.items())}}
    return chosen


if __name__ == "__main__":
    main()

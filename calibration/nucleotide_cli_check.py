"""False-positive rate of `rna_lexis_stat_cli rank-cores` on iid nucleotide text, for several settings.

Runs the same two calls as the CLI (`_discover_xmotifs`, then `rank_core_candidates` with the CLI's
arguments), on random text whose composition is that of a human coding sequence, and reports the share
of texts in which at least one core is called "statistically supported".  Nothing is planted, so every
call is a chance finding.

Usage: python3 nucleotide_cli_check.py [--reps 300]
Writes results/nucleotide_cli_check.csv.
"""
import argparse
import contextlib
import csv
import io
import multiprocessing as mp
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))
sys.path.insert(0, str(HERE))
import harness as H                                              # noqa: E402
from rna_lexis.statistical import rank_core_candidates           # noqa: E402
from rna_lexis.test_cli import _discover_xmotifs                 # noqa: E402

# (label, --min-xmotif-len, --min-count, --candidate-min-len); the first is the CLI default
SETTINGS = [("CLI defaults", 7, 2, 5), ("stricter", 8, 3, 7)]
LENGTHS = [900, 2000, 4500]
_comp = {}


def _init():
    genes = H.load_genes(HERE / "data" / "human_cds.fasta")
    nt = np.zeros(4)
    for u in genes.values():
        for c, n in zip(*np.unique(u, return_counts=True)):
            for j in range(3):
                nt[H.NT_OF[c, j]] += n
    _comp["p"] = nt / nt.sum()


def task(args):
    n, rep = args
    rng = np.random.default_rng([7, n, rep])
    txt = "".join("acgt"[i] for i in rng.choice(4, size=n, p=_comp["p"]))
    out = {}
    for label, mx, mc, cm in SETTINGS:
        with contextlib.redirect_stdout(io.StringIO()):
            xm = _discover_xmotifs(txt, mx, 60, mc)
            rows = rank_core_candidates(txt, xm, candidate_min_len=cm, candidate_max_len=18,
                                        min_xmotif_type_support=2, markov_order=1, alpha=0.05,
                                        enrichment_threshold=10.0, min_count=mc) if xm else []
        out[label] = (len(rows) > 0, any(r["statistically_supported"] for r in rows))
    return n, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=300)
    a = ap.parse_args()
    tasks = [(n, r) for n in LENGTHS for r in range(a.reps)]
    agg = {}
    with mp.Pool(max(1, (os.cpu_count() or 2) - 1), initializer=_init) as pool:
        for n, out in pool.imap_unordered(task, tasks):
            for label, (cand, sup) in out.items():
                d = agg.setdefault((n, label), [0, 0, 0])
                d[0] += 1
                d[1] += cand
                d[2] += sup
    path = HERE / "results" / "nucleotide_cli_check.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["n_nt", "settings", "min_xmotif_len", "min_count", "candidate_min_len", "texts",
                    "frac_with_candidates", "frac_with_supported_core"])
        for (n, label), (t, c, s) in sorted(agg.items()):
            mx, mc, cm = next((x, y, z) for l, x, y, z in SETTINGS if l == label)
            w.writerow([n, label, mx, mc, cm, t, round(c / t, 4), round(s / t, 4)])
    print(open(path, encoding="utf-8").read())


if __name__ == "__main__":
    main()

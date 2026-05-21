"""Command-line entry point for RNA-Lexis-Stat workflows."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from rna_lexis.algorithms import cores, find_boundary
from rna_lexis.statistical import (
    best_mutation_family_per_motif,
    find_gapped_motif_hits,
    mutation_family_tests,
    normalize_sequence,
    rank_core_candidates,
    score_gapped_motif,
    score_exact_motifs,
    write_rows_csv,
)


def _read_sequence(path: str) -> str:
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    if text.lstrip().startswith(">"):
        text = "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))
    return normalize_sequence(text)


def _read_list(path: str, column: str | None = None) -> list[str]:
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    if path.lower().endswith(".json"):
        data = json.loads(text)
        if column:
            values = data[column]
        else:
            values = data.get("corelist") or data.get("xmotifs") or []
        return [normalize_sequence(str(x)) for x in values if str(x).strip()]

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return []
    if "\t" in lines[0] or "," in lines[0]:
        sep = "\t" if "\t" in lines[0] else ","
        header = [h.strip() for h in lines[0].split(sep)]
        idx = header.index(column) if column and column in header else 0
        return [normalize_sequence(row.split(sep)[idx]) for row in lines[1:] if row.split(sep)[idx].strip()]
    return [normalize_sequence(line.split()[0]) for line in lines]


def _discover_xmotifs(sequence: str, min_xmotif_len: int, max_xmotif_len: int, min_count: int) -> list[str]:
    xmotifs = find_boundary(sequence, minlen=min_xmotif_len, maxlen=max_xmotif_len, at_least=min_count)
    return sorted(set(xmotifs), key=lambda x: (-len(x), x))


def cmd_score_exact(args: argparse.Namespace) -> None:
    sequence = _read_sequence(args.input)
    motifs = _read_list(args.motifs, args.motif_column)
    xmotifs = _read_list(args.xmotifs, args.xmotif_column) if args.xmotifs else []
    rows = score_exact_motifs(
        sequence,
        motifs,
        xmotifs=xmotifs,
        markov_order=args.markov_order,
        alpha=args.alpha,
        enrichment_threshold=args.enrichment_threshold,
        min_count=args.min_count,
        min_xmotif_type_support=args.min_xmotif_type_support,
    )
    write_rows_csv(args.output, rows)
    print(f"Wrote {len(rows)} exact motif score rows to {args.output}")


def cmd_rank_cores(args: argparse.Namespace) -> None:
    sequence = _read_sequence(args.input)
    if args.xmotifs:
        xmotifs = _read_list(args.xmotifs, args.xmotif_column)
    else:
        xmotifs = _discover_xmotifs(sequence, args.min_xmotif_len, args.max_xmotif_len, args.min_count)
    rows = rank_core_candidates(
        sequence,
        xmotifs,
        candidate_min_len=args.candidate_min_len,
        candidate_max_len=args.candidate_max_len,
        min_xmotif_type_support=args.min_xmotif_type_support,
        markov_order=args.markov_order,
        alpha=args.alpha,
        enrichment_threshold=args.enrichment_threshold,
        min_count=args.min_count,
    )
    write_rows_csv(args.output, rows)
    print(f"Wrote {len(rows)} ranked core candidate rows to {args.output}")

    if args.legacy_cores_output:
        legacy = cores(sequence, xmotifs, minclen=args.candidate_min_len)
        write_rows_csv(args.legacy_cores_output, [{"core": c.upper()} for c in legacy])
        print(f"Wrote {len(legacy)} legacy core rows to {args.legacy_cores_output}")


def cmd_mutation_families(args: argparse.Namespace) -> None:
    sequence = _read_sequence(args.input)
    motifs = _read_list(args.motifs, args.motif_column)
    rows = mutation_family_tests(
        sequence,
        motifs,
        mutr=args.mutr,
        M=args.max_mismatches,
        markov_order=args.markov_order,
        alpha=args.alpha,
        enrichment_threshold=args.enrichment_threshold,
        expected_max=args.expected_max,
        min_family_count=args.min_family_count,
        max_variants=args.max_variants,
    )
    write_rows_csv(args.output, rows)
    print(f"Wrote {len(rows)} motif-radius test rows to {args.output}")

    if args.best_output:
        best = best_mutation_family_per_motif(rows)
        write_rows_csv(args.best_output, best)
        print(f"Wrote {len(best)} best-radius rows to {args.best_output}")


def cmd_gapped_motif(args: argparse.Namespace) -> None:
    sequence = _read_sequence(args.input)
    score = score_gapped_motif(
        sequence,
        args.left,
        args.right,
        min_gap=args.min_gap,
        max_gap=args.max_gap,
        markov_order=args.markov_order,
    )
    hits = find_gapped_motif_hits(
        sequence,
        args.left,
        args.right,
        min_gap=args.min_gap,
        max_gap=args.max_gap,
    )
    rows = [{**score, **hit} for hit in hits]
    if not rows:
        rows = [score]
    write_rows_csv(args.output, rows)
    print(f"Wrote {len(hits)} gapped motif hit row(s) to {args.output}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rna_lexis_stat_cli",
        description="RNA-Lexis-Stat workflow: Markov/FDR scoring, mutation families, and gapped motifs.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def add_common_score_options(p: argparse.ArgumentParser) -> None:
        p.add_argument("--markov-order", type=int, default=1)
        p.add_argument("--alpha", type=float, default=0.05)
        p.add_argument("--enrichment-threshold", type=float, default=10.0)
        p.add_argument("--min-count", type=int, default=2)

    exact = sub.add_parser("score-exact", help="Score provided exact motifs with Markov/FDR statistics.")
    exact.add_argument("--input", required=True, help="FASTA or plain sequence file.")
    exact.add_argument("--motifs", required=True, help="Motif list/CSV/TSV/session JSON.")
    exact.add_argument("--motif-column", default=None)
    exact.add_argument("--xmotifs", default=None, help="Optional xmotif list/CSV/TSV/session JSON.")
    exact.add_argument("--xmotif-column", default=None)
    exact.add_argument("--min-xmotif-type-support", type=int, default=0)
    exact.add_argument("--output", required=True)
    add_common_score_options(exact)
    exact.set_defaults(func=cmd_score_exact)

    rank = sub.add_parser("rank-cores", help="Generate and rank shared core candidates from xmotifs.")
    rank.add_argument("--input", required=True, help="FASTA or plain sequence file.")
    rank.add_argument("--xmotifs", default=None, help="Optional xmotif list/CSV/TSV/session JSON.")
    rank.add_argument("--xmotif-column", default=None)
    rank.add_argument("--candidate-min-len", type=int, default=5)
    rank.add_argument("--candidate-max-len", type=int, default=18)
    rank.add_argument("--min-xmotif-type-support", type=int, default=2)
    rank.add_argument("--min-xmotif-len", type=int, default=7)
    rank.add_argument("--max-xmotif-len", type=int, default=60)
    rank.add_argument("--output", required=True)
    rank.add_argument("--legacy-cores-output", default=None)
    add_common_score_options(rank)
    rank.set_defaults(func=cmd_rank_cores)

    mut = sub.add_parser("mutation-families", help="Test Hamming-radius motif families with Markov/FDR support.")
    mut.add_argument("--input", required=True, help="FASTA or plain sequence file.")
    mut.add_argument("--motifs", required=True, help="Motif list/CSV/TSV/session JSON.")
    mut.add_argument("--motif-column", default=None)
    mut.add_argument("--mutr", type=float, default=1 / 6)
    mut.add_argument("--max-mismatches", type=int, default=4)
    mut.add_argument("--markov-order", type=int, default=1)
    mut.add_argument("--alpha", type=float, default=0.05)
    mut.add_argument("--enrichment-threshold", type=float, default=5.0)
    mut.add_argument("--expected-max", type=float, default=5.0)
    mut.add_argument("--min-family-count", type=int, default=2)
    mut.add_argument("--max-variants", type=int, default=500_000)
    mut.add_argument("--output", required=True)
    mut.add_argument("--best-output", default=None)
    mut.set_defaults(func=cmd_mutation_families)

    gap = sub.add_parser("gapped-motif", help="Search an anchor-gap-anchor motif and score it.")
    gap.add_argument("--input", required=True, help="FASTA or plain sequence file.")
    gap.add_argument("--left", required=True, help="Left exact anchor sequence.")
    gap.add_argument("--right", required=True, help="Right exact anchor sequence.")
    gap.add_argument("--min-gap", type=int, default=0)
    gap.add_argument("--max-gap", type=int, default=30)
    gap.add_argument("--markov-order", type=int, default=1)
    gap.add_argument("--output", required=True)
    gap.set_defaults(func=cmd_gapped_motif)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    args.func(args)


if __name__ == "__main__":
    main()

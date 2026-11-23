"""Command-line entry point for RNA-Lexis-Stat workflows."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import warnings
from pathlib import Path

from rna_lexis.algorithms import cores, find_boundary
from rna_lexis.alphabet import get_alphabet
from rna_lexis.tokenizers import UnknownSymbolError, tokenize
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


_BUILTIN_ALPHABETS = ("nucleotide", "codon", "Latin letters")


def _alphabet(args):
    """The Alphabet named by --alphabet (a built-in name or an alphabet file), or None for nucleotides."""
    return get_alphabet(args.alphabet) if getattr(args, "alphabet", None) else None


def _tokenize(text, alphabet, args, frame=0):
    """Read text into the alphabet, printing what was skipped or dropped to stderr.

    The reading frame belongs to the sequence only; motifs and anchors are read from their start.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        tok = tokenize(text, alphabet, frame=frame, on_unknown=args.on_unknown)
    for w in caught:
        print(f"warning: {w.message}", file=sys.stderr)
    return tok


def _read_sequence(path: str, alphabet=None, args=None) -> str:
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    if alphabet is None:
        if text.lstrip().startswith(">"):
            text = "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))
        return normalize_sequence(text)
    # Over an alphabet each FASTA record is a segment, so no motif spans two records.
    if text.lstrip().startswith(">"):
        records = [r for r in re.split(r"^>.*$", text, flags=re.M) if r.strip()]
    else:
        records = [text]
    if alphabet.name in _BUILTIN_ALPHABETS:        # letters only, as when the menu loads text
        records = [re.sub(r"[^A-Za-z]+", "", r).lower() for r in records]
    tok = _tokenize(records, alphabet, args, frame=args.frame)
    if tok.n_units == 0:
        sys.exit(f"error: no {alphabet.units_label} could be read from {path} with the alphabet {alphabet.name!r}")
    print(f"Read {tok.n_units:,} {alphabet.units_label} over {alphabet.name!r} from {path}"
          f" ({len(tok.skipped)} skipped, {sum(len(d.text) for d in tok.dropped)} characters dropped)")
    return tok.text


def _to_internal(values, alphabet, args):
    out = []
    for value in values:
        text = _tokenize(str(value), alphabet, args).text
        if text:
            out.append(text)
    return out


def _read_list(path: str, column: str | None = None, alphabet=None, args=None) -> list[str]:
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    if path.lower().endswith(".json"):
        data = json.loads(text)
        if column:
            values = data[column]
        else:
            values = data.get("corelist") or data.get("xmotifs") or []
        if alphabet is not None:
            stored = data.get("alphabet") if isinstance(data, dict) else None
            if stored is not None:
                # a saved session: its motifs are already internal strings of its own alphabet
                if stored.get("name") != alphabet.name:
                    sys.exit(f"error: {path} was read over {stored.get('name')!r}, not {alphabet.name!r}")
                return [str(x) for x in values if str(x).strip()]
            return _to_internal([x for x in values if str(x).strip()], alphabet, args)
        return [normalize_sequence(str(x)) for x in values if str(x).strip()]

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return []
    if "\t" in lines[0] or "," in lines[0]:
        sep = "\t" if "\t" in lines[0] else ","
        header = [h.strip() for h in lines[0].split(sep)]
        idx = header.index(column) if column and column in header else 0
        values = [row.split(sep)[idx] for row in lines[1:] if row.split(sep)[idx].strip()]
    else:
        values = [line.split()[0] for line in lines]
    if alphabet is not None:
        return _to_internal(values, alphabet, args)
    return [normalize_sequence(v) for v in values]


def _discover_xmotifs(sequence: str, min_xmotif_len: int, max_xmotif_len: int, min_count: int) -> list[str]:
    xmotifs = find_boundary(sequence, minlen=min_xmotif_len, maxlen=max_xmotif_len, at_least=min_count)
    return sorted(set(xmotifs), key=lambda x: (-len(x), x))


def cmd_score_exact(args: argparse.Namespace) -> None:
    alphabet = _alphabet(args)
    sequence = _read_sequence(args.input, alphabet, args)
    motifs = _read_list(args.motifs, args.motif_column, alphabet, args)
    xmotifs = _read_list(args.xmotifs, args.xmotif_column, alphabet, args) if args.xmotifs else []
    rows = score_exact_motifs(
        sequence,
        motifs,
        xmotifs=xmotifs,
        markov_order=args.markov_order,
        alpha=args.alpha,
        enrichment_threshold=args.enrichment_threshold,
        min_count=args.min_count,
        min_xmotif_type_support=args.min_xmotif_type_support,
        alphabet=alphabet,
    )
    write_rows_csv(args.output, rows, alphabet=alphabet)
    print(f"Wrote {len(rows)} exact motif score rows to {args.output}")


def cmd_rank_cores(args: argparse.Namespace) -> None:
    alphabet = _alphabet(args)
    sequence = _read_sequence(args.input, alphabet, args)
    if args.xmotifs:
        xmotifs = _read_list(args.xmotifs, args.xmotif_column, alphabet, args)
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
        alphabet=alphabet,
    )
    write_rows_csv(args.output, rows, alphabet=alphabet)
    print(f"Wrote {len(rows)} ranked core candidate rows to {args.output}")

    if args.legacy_cores_output:
        legacy = cores(sequence, xmotifs, minclen=args.candidate_min_len)
        write_rows_csv(args.legacy_cores_output, [{"core": c.upper()} for c in legacy], alphabet=alphabet)
        print(f"Wrote {len(legacy)} legacy core rows to {args.legacy_cores_output}")


def cmd_mutation_families(args: argparse.Namespace) -> None:
    alphabet = _alphabet(args)
    sequence = _read_sequence(args.input, alphabet, args)
    motifs = _read_list(args.motifs, args.motif_column, alphabet, args)
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
        alphabet=alphabet,
    )
    write_rows_csv(args.output, rows, alphabet=alphabet)
    print(f"Wrote {len(rows)} motif-radius test rows to {args.output}")

    if args.best_output:
        best = best_mutation_family_per_motif(rows)
        write_rows_csv(args.best_output, best, alphabet=alphabet)
        print(f"Wrote {len(best)} best-radius rows to {args.best_output}")


def cmd_gapped_motif(args: argparse.Namespace) -> None:
    alphabet = _alphabet(args)
    sequence = _read_sequence(args.input, alphabet, args)
    left, right = args.left, args.right
    if alphabet is not None:
        left, right = _tokenize(left, alphabet, args).text, _tokenize(right, alphabet, args).text
        if not left or not right:
            sys.exit("error: an anchor has no symbols of the alphabet")
    score = score_gapped_motif(
        sequence,
        left,
        right,
        min_gap=args.min_gap,
        max_gap=args.max_gap,
        markov_order=args.markov_order,
        alphabet=alphabet,
    )
    hits = find_gapped_motif_hits(
        sequence,
        left,
        right,
        min_gap=args.min_gap,
        max_gap=args.max_gap,
        alphabet=alphabet,
    )
    rows = [{**score, **hit} for hit in hits]
    if not rows:
        rows = [score]
    write_rows_csv(args.output, rows, alphabet=alphabet)
    print(f"Wrote {len(hits)} gapped motif hit row(s) to {args.output}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rna_lexis_stat_cli",
        description="RNA-Lexis-Stat workflow: Markov/FDR scoring, mutation families, and gapped motifs.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def add_alphabet_options(p: argparse.ArgumentParser) -> None:
        p.add_argument("--alphabet", default=None,
                       help="Read the text over an alphabet: nucleotide, codon, latin, or the path of an "
                            "alphabet file (JSON). Default: nucleotides, as always. Motifs, anchors and "
                            "queries are then read the same way (for codons, as nucleotides), and "
                            "results are shown as the alphabet's labels.")
        p.add_argument("--frame", type=int, default=0,
                       help="Characters to skip at the start of each sequence for a fixed-width alphabet "
                            "such as codon (the reading frame; default 0).")
        p.add_argument("--on-unknown", choices=["skip", "error"], default="skip",
                       help="What to do with a symbol that is not in the alphabet: skip it with a warning "
                            "(default), or stop with an error.")

    def add_common_score_options(p: argparse.ArgumentParser) -> None:
        p.add_argument("--markov-order", type=int, default=None,
                       help="Order of the Markov background. Default: 1 for nucleotides; with --alphabet, "
                            "order 1 when the text has at least one letter pair per transition "
                            "(length >= K**2), otherwise order 0.")
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
    mut.add_argument("--markov-order", type=int, default=None)
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
    gap.add_argument("--markov-order", type=int, default=None)
    gap.add_argument("--output", required=True)
    gap.set_defaults(func=cmd_gapped_motif)

    for p in (exact, rank, mut, gap):
        add_alphabet_options(p)
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    try:
        args.func(args)
    except UnknownSymbolError as exc:
        sys.exit(f"error: {exc}")
    except ValueError as exc:
        if getattr(args, "alphabet", None):
            sys.exit(f"error: {exc}")
        raise


if __name__ == "__main__":
    main()

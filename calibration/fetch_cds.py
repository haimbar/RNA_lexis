"""Fetch the coding sequence (CDS) of the canonical transcript of a list of human genes.

Public Ensembl REST data; only gene symbols are sent.  Writes data/human_cds.fasta
(one record per gene, header ">SYMBOL ENST... length_in_codons") and prints a summary.
A CDS is kept only if it is a whole number of codons, starts with ATG, ends in a stop
codon, has no stop codon inside and no ambiguity code; the stop codon is removed.

Usage: python3 fetch_cds.py
"""
import sys
import time
import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from rna_lexis.io import _fetch_ensembl_json

GENES = """
ACTB GAPDH TP53 MYC EGFR KRAS NRAS HRAS BRAF PIK3CA PTEN RB1 APC BRCA1 BRCA2 ATM CDKN2A CDK4 CDK6
CCND1 CCNE1 MDM2 BCL2 BAX CASP3 CASP8 CASP9 FAS TNF IL6 IL1B IFNG TGFB1 SMAD4 SMAD2 VEGFA FGF2
PDGFRA KIT ALK RET MET ERBB2 ERBB3 JAK2 STAT3 STAT1 NFKB1 RELA MAPK1 MAPK3 MAP2K1 AKT1 MTOR TSC1
TSC2 VHL HIF1A NOTCH1 WNT1 CTNNB1 GSK3B HDAC1 EZH2 DNMT1 DNMT3A TET2 IDH1 IDH2 SOX2 POU5F1 NANOG
KLF4 GATA1 RUNX1 FOXP3 PAX6 HOXA9 MYOD1 TTN DMD COL1A1 COL1A2 COL3A1 FBN1 ELN LMNA KRT5 KRT14 VIM
ACTA1 MYH7 MYBPC3 TPM1 TNNT2 HBB HBA1 INS GCG ALB TF APOE APOB LDLR PCSK9 CFTR SOD1 HTT APP PSEN1
MAPT SNCA PRKN LRRK2 FMR1 MECP2 SCN1A KCNQ1 RYR1 CACNA1C GRIN1 GABRA1 DRD2 HTR2A SLC6A4 BDNF NGF
CREB1 FOS JUN EGR1 ATF4 HSP90AA1 HSPA1A HSPB1 UBB UBC RPS6 RPL13A EEF1A1 TUBB TUBA1B CALM1 PPIA B2M
HLA-A CD4 CD8A CD19 MS4A1 TLR4 MYD88 CGAS STING1 DDX58 IFIH1 IRF3 ISG15 MX1 OAS1 TRIM5 SAMHD1 TP63
TP73 CHEK1 CHEK2 PARP1 RAD51 MLH1 MSH2 XRCC1 ERCC1 FANCA SPTAN1 AHNAK NEB MUC1 MUC5AC FLNA TNC
""".split()

STOPS = {"TAA", "TAG", "TGA"}
OUT = Path(__file__).resolve().parent / "data" / "human_cds.fasta"


def fetch(symbol):
    info = _fetch_ensembl_json(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?expand=0")
    tid = (info.get("canonical_transcript") or "").split(".")[0]
    if not tid:
        raise ValueError("no canonical transcript")
    seq = _fetch_ensembl_json(f"https://rest.ensembl.org/sequence/id/{tid}?type=cds")["seq"].upper()
    return tid, seq


def valid(seq):
    if len(seq) % 3 or not seq.startswith("ATG") or set(seq) - set("ACGT"):
        return False
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    return codons[-1] in STOPS and not any(c in STOPS for c in codons[:-1])


def main():
    kept, rejected = [], []
    for sym in GENES:
        try:
            tid, seq = fetch(sym)
        except Exception as exc:                       # network error, unknown symbol, ...
            rejected.append((sym, f"{type(exc).__name__}: {str(exc)[:60]}"))
            continue
        finally:
            time.sleep(0.12)                           # stay under Ensembl's rate limit
        if valid(seq):
            kept.append((sym, tid, seq[:-3]))
        else:
            rejected.append((sym, "CDS not a clean ATG...stop reading frame"))
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(f"; Ensembl REST (canonical transcripts, homo_sapiens), fetched {datetime.date.today()}\n")
        for sym, tid, seq in kept:
            fh.write(f">{sym} {tid} {len(seq) // 3}\n{seq}\n")
    n = sum(len(s) // 3 for _, _, s in kept)
    print(f"kept {len(kept)} genes, {n} codons; rejected {len(rejected)}")
    for r in rejected:
        print("  rejected", *r)


if __name__ == "__main__":
    main()

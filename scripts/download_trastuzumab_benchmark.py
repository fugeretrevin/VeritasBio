"""
scripts/download_trastuzumab_benchmark.py — Downloads and prepares the authentic Trastuzumab CDR-H3 dataset.

Dataset details:
    Target: Trastuzumab (Herceptin) therapeutic monoclonal antibody heavy chain (VH domain)
    Parent: Trastuzumab humanized Fab targeting human HER2 receptor (PDB: 1N8Z)
    CDR-H3 WT: WGGDGFYAMD (10 residues, positions 99-108 in VH domain)
    Library: Deep mutational scan of the 10-residue CDR-H3 loop
             (38,839 empirical variants: 11,300 binders, 27,539 non-binders).
    Reference: Mason, D. et al. (2021), "Optimization of therapeutic antibodies by predicting
               antigen specificity using machine learning", Nature Biomedical Engineering.
    Source: dahjan/DMS_opt (mHER_H3_AgPos.csv & mHER_H3_AgNeg.csv)
"""

from __future__ import annotations

import csv
import io
import json
import logging
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("download_trastuzumab_benchmark")

ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = Path(__file__).resolve().parent / "data"
WORKER_DATA_DIR = ROOT_DIR / "workers" / "data"

OUTPUT_CSV = DATA_DIR / "HER2_TRASTUZUMAB_Mason_2021.csv"
CACHE_JSON = WORKER_DATA_DIR / "HER2_TRASTUZUMAB.json"

POS_URL = "https://raw.githubusercontent.com/dahjan/DMS_opt/master/data/mHER_H3_AgPos.csv"
NEG_URL = "https://raw.githubusercontent.com/dahjan/DMS_opt/master/data/mHER_H3_AgNeg.csv"

# Trastuzumab Variable Heavy (VH) domain (120 residues)
TRASTUZUMAB_VH_SEQ = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
CDRH3_WT_SEQ = "WGGDGFYAMD"
CDRH3_START_POS = 99  # 1-indexed start position of WGGDGFYAMD in TRASTUZUMAB_VH_SEQ
CDRH3_POSITIONS = list(range(99, 109))  # [99, 100, 101, 102, 103, 104, 105, 106, 107, 108]


def cdr_to_mutation_str(cdr_seq: str) -> str:
    """Diffs a 10-residue CDR-H3 sequence against WT WGGDGFYAMD to return standard notation."""
    if len(cdr_seq) != 10:
        return cdr_seq
    if cdr_seq == CDRH3_WT_SEQ:
        return "WT"
    muts = []
    for i, (wt_c, mut_c) in enumerate(zip(CDRH3_WT_SEQ, cdr_seq)):
        if wt_c != mut_c:
            pos = CDRH3_START_POS + i
            muts.append(f"{wt_c}{pos}{mut_c}")
    return ":".join(muts) if muts else "WT"


def cdr_to_full_vh_sequence(cdr_seq: str) -> str:
    """Splices the 10-residue CDR-H3 sequence into the full 120-residue VH domain."""
    prefix = TRASTUZUMAB_VH_SEQ[:CDRH3_START_POS - 1]
    suffix = TRASTUZUMAB_VH_SEQ[CDRH3_START_POS - 1 + len(CDRH3_WT_SEQ):]
    return prefix + cdr_seq + suffix


def download_trastuzumab_dataset(force: bool = False) -> Path:
    """Downloads positive and negative Trastuzumab screens, standardizes to ProteinGym format, and caches metadata."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    WORKER_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if OUTPUT_CSV.exists() and not force and OUTPUT_CSV.stat().st_size > 1_000_000:
        logger.info("Trastuzumab dataset already exists at %s (%d bytes). Skipping download.", OUTPUT_CSV, OUTPUT_CSV.stat().st_size)
        write_metadata_cache()
        return OUTPUT_CSV

    logger.info("Downloading Trastuzumab AgPos from %s...", POS_URL)
    req_pos = urllib.request.Request(POS_URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req_pos, timeout=25.0) as resp:
        pos_content = resp.read().decode("utf-8")

    logger.info("Downloading Trastuzumab AgNeg from %s...", NEG_URL)
    req_neg = urllib.request.Request(NEG_URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req_neg, timeout=25.0) as resp:
        neg_content = resp.read().decode("utf-8")

    rows_written = 0
    seen_mutants = set()

    with open(OUTPUT_CSV, "w", encoding="utf-8", newline="") as out_f:
        writer = csv.writer(out_f)
        writer.writerow(["mutant", "mutated_sequence", "DMS_score", "DMS_score_bin", "stage"])

        # 1. Parse Binders (AgClass = 1)
        reader_pos = csv.DictReader(io.StringIO(pos_content))
        for r in reader_pos:
            aa_seq = r.get("AASeq", "").strip().upper()
            if len(aa_seq) != 10:
                continue
            mut_str = cdr_to_mutation_str(aa_seq)
            if mut_str in seen_mutants:
                continue
            seen_mutants.add(mut_str)

            try:
                frac = float(r.get("Fraction", 0.0))
            except ValueError:
                frac = 0.0

            # Scale fractional antigen binding to continuous score: 1.0 to 10.0 range
            score = max(0.1, round(frac * 1000.0, 4))
            full_seq = cdr_to_full_vh_sequence(aa_seq)
            writer.writerow([mut_str, full_seq, f"{score:.6f}", 1, "train"])
            rows_written += 1

        # 2. Parse Non-Binders (AgClass = 0)
        reader_neg = csv.DictReader(io.StringIO(neg_content))
        for r in reader_neg:
            aa_seq = r.get("AASeq", "").strip().upper()
            if len(aa_seq) != 10:
                continue
            mut_str = cdr_to_mutation_str(aa_seq)
            if mut_str in seen_mutants:
                continue
            seen_mutants.add(mut_str)

            full_seq = cdr_to_full_vh_sequence(aa_seq)
            writer.writerow([mut_str, full_seq, "0.000000", 0, "test"])
            rows_written += 1

    logger.info("Successfully wrote %d Trastuzumab variants to %s.", rows_written, OUTPUT_CSV)
    write_metadata_cache()
    return OUTPUT_CSV


def write_metadata_cache() -> None:
    """Writes the offline UniProt metadata cache for HER2_TRASTUZUMAB conforming to parse_uniprot_response schema."""
    features = [
        {
            "type": "Domain",
            "description": "Variable Heavy Domain (VH)",
            "location": {"start": {"value": 1}, "end": {"value": len(TRASTUZUMAB_VH_SEQ)}},
        }
    ]
    for pos in CDRH3_POSITIONS:
        features.append({
            "type": "Binding site",
            "description": "Complementarity-Determining Region 3 (CDR-H3) antigen contact",
            "location": {"start": {"value": pos}, "end": {"value": pos}},
        })
    for pos in [1, 2, 22, 92, 119, 120]:
        features.append({
            "type": "Active site",
            "description": "Invariant framework anchor / structural cysteine",
            "location": {"start": {"value": pos}, "end": {"value": pos}},
        })

    record = {
        "primaryAccession": "HER2_TRASTUZUMAB",
        "uniProtkbId": "TRASTUZUMAB_VH",
        "proteinDescription": {
            "recommendedName": {"fullName": {"value": "Trastuzumab (Herceptin) Fab Variable Heavy Domain"}}
        },
        "organism": {"scientificName": "Homo sapiens / Synthetic therapeutic antibody"},
        "sequence": {"value": TRASTUZUMAB_VH_SEQ, "length": len(TRASTUZUMAB_VH_SEQ)},
        "uniProtKBCrossReferences": [
            {"database": "PDB", "id": "1N8Z"}
        ],
        "features": features,
    }
    with open(CACHE_JSON, "w", encoding="utf-8") as f:
        json.dump(record, f, indent=2)
    logger.info("Cached HER2_TRASTUZUMAB metadata at %s.", CACHE_JSON)


if __name__ == "__main__":
    download_trastuzumab_dataset(force=("--force" in sys.argv))

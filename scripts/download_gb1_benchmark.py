"""
scripts/download_gb1_benchmark.py — Downloads and prepares the authentic GB1 fitness landscape.

Dataset details:
    Target: Protein G B1 domain (GB1), Streptococcus dysgalactiae / Streptococcus group G
    Parent: SPG1_STRSG wild type (UniProt P06654, 56 residues:
            MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE)
    Library: Combinatorially complete 4-site saturation at residues 39, 40, 41, 54
             (20^4 = 160,000 variants, 149,361 empirical measurements).
    Reference: Wu, Dai et al. (2016), "Adaptation in protein fitness landscapes is
               facilitated by indirect paths", eLife 5:e16965.
    Source: FLIP Benchmark (four_mutations_full_data.csv.zip)
"""

from __future__ import annotations

import csv
import io
import json
import logging
import sys
import urllib.request
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("download_gb1_benchmark")

ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = Path(__file__).resolve().parent / "data"
WORKER_DATA_DIR = ROOT_DIR / "workers" / "data"

OUTPUT_CSV = DATA_DIR / "GB1_STRSG_Wu_2016.csv"
CACHE_JSON = WORKER_DATA_DIR / "SPG1_STRSG.json"

GB1_ZIP_URL = "https://github.com/J-SNACKKB/FLIP/raw/main/splits/gb1/four_mutations_full_data.csv.zip"

# Wild-type 56-amino acid sequence of GB1 domain
GB1_WT_SEQ = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
GB1_POSITIONS = [39, 40, 41, 54]
GB1_WT_ALLELES = ["V", "D", "G", "V"]  # VDGV


def alleles_to_mutation_str(alleles: str) -> str:
    """Converts a 4-letter allele string like 'ADGV' to 'V39A' or 'WT'."""
    if len(alleles) != 4:
        return alleles
    if alleles == "VDGV":
        return "WT"
    muts = []
    for wt_c, mut_c, pos in zip(GB1_WT_ALLELES, alleles, GB1_POSITIONS):
        if wt_c != mut_c:
            muts.append(f"{wt_c}{pos}{mut_c}")
    return ":".join(muts) if muts else "WT"


def alleles_to_full_sequence(alleles: str) -> str:
    """Substitutes the 4 alleles into the 56-residue WT sequence."""
    seq_chars = list(GB1_WT_SEQ)
    for pos, mut_c in zip(GB1_POSITIONS, alleles):
        seq_chars[pos - 1] = mut_c
    return "".join(seq_chars)


def download_gb1_dataset(force: bool = False) -> Path:
    """Downloads FLIP zip, standardizes to ProteinGym format, and caches metadata."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    WORKER_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if OUTPUT_CSV.exists() and not force and OUTPUT_CSV.stat().st_size > 1_000_000:
        logger.info("GB1 dataset already exists at %s (%d bytes). Skipping download.", OUTPUT_CSV, OUTPUT_CSV.stat().st_size)
        write_metadata_cache()
        return OUTPUT_CSV

    logger.info("Downloading GB1 4-site combinatorial landscape from %s...", GB1_ZIP_URL)
    req = urllib.request.Request(GB1_ZIP_URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30.0) as resp:
        zip_bytes = resp.read()

    logger.info("Downloaded %d bytes. Unzipping and parsing...", len(zip_bytes))
    rows_written = 0
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as z:
        content = z.read("four_mutations_full_data.csv").decode("utf-8")
        reader = csv.DictReader(io.StringIO(content))

        with open(OUTPUT_CSV, "w", encoding="utf-8", newline="") as out_f:
            writer = csv.writer(out_f)
            writer.writerow(["mutant", "mutated_sequence", "DMS_score", "DMS_score_bin", "stage"])

            for r in reader:
                variants = r.get("Variants", "").strip().upper()
                raw_fit = r.get("Fitness", "").strip()
                if not variants or not raw_fit:
                    continue
                try:
                    fit = float(raw_fit)
                except ValueError:
                    continue

                mut_str = alleles_to_mutation_str(variants)
                full_seq = alleles_to_full_sequence(variants)
                score_bin = 1 if fit >= 1.0 else 0
                stage = r.get("sampled", "train")

                writer.writerow([mut_str, full_seq, f"{fit:.6f}", score_bin, stage])
                rows_written += 1

    logger.info("Successfully wrote %d GB1 variants to %s.", rows_written, OUTPUT_CSV)
    write_metadata_cache()
    return OUTPUT_CSV


def write_metadata_cache() -> None:
    """Writes the offline UniProt metadata cache for SPG1_STRSG conforming to parse_uniprot_response schema."""
    features = [
        {
            "type": "Domain",
            "description": "IgG-binding domain B1",
            "location": {"start": {"value": 1}, "end": {"value": len(GB1_WT_SEQ)}},
        }
    ]
    for pos in GB1_POSITIONS:
        features.append({
            "type": "Binding site",
            "description": "IgG Fc binding contact residue (Wu 2016)",
            "location": {"start": {"value": pos}, "end": {"value": pos}},
        })
    for pos in [1, len(GB1_WT_SEQ)]:
        features.append({
            "type": "Active site",
            "description": "Invariant framework terminal anchor",
            "location": {"start": {"value": pos}, "end": {"value": pos}},
        })

    record = {
        "primaryAccession": "SPG1_STRSG",
        "uniProtkbId": "SPG1_STRSG_GB1",
        "proteinDescription": {
            "recommendedName": {"fullName": {"value": "Immunoglobulin G-binding protein G (GB1 domain)"}}
        },
        "organism": {"scientificName": "Streptococcus dysgalactiae"},
        "sequence": {"value": GB1_WT_SEQ, "length": len(GB1_WT_SEQ)},
        "uniProtKBCrossReferences": [
            {"database": "PDB", "id": "2QMT"}
        ],
        "features": features,
    }
    with open(CACHE_JSON, "w", encoding="utf-8") as f:
        json.dump(record, f, indent=2)
    logger.info("Cached SPG1_STRSG metadata at %s.", CACHE_JSON)


if __name__ == "__main__":
    download_gb1_dataset(force=("--force" in sys.argv))

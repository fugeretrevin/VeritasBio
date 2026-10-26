"""
scripts/download_real_benchmark.py — Downloads and prepares the authentic ProteinGym GFP dataset.

Dataset details:
    Target: Green Fluorescent Protein (Aequorea victoria avGFP)
    Reference: Sarkisyan et al. (2016) Nature 533, 397-401
    Source: Harvard Marks Lab / ProteinGym Benchmark
    Total variants: ~51,715 empirical FACS-seq functional measurements
"""

from __future__ import annotations

import csv
import logging
import os
import sys
import urllib.request
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("download_benchmark")

DATA_DIR = Path(__file__).resolve().parent / "data"
OUTPUT_CSV = DATA_DIR / "GFP_AEQVI_Sarkisyan_2016.csv"

# Official HuggingFace ProteinGym mirror
PROTEINGYM_URL = (
    "https://huggingface.co/datasets/ICML2022/ProteinGym/resolve/main/"
    "ProteinGym_substitutions/GFP_AEQVI_Sarkisyan_2016.csv"
)


def download_proteingym_gfp(dest: Path = OUTPUT_CSV, force: bool = False) -> Path:
    """Downloads the official ProteinGym GFP substitution assay CSV if not present."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force and dest.stat().st_size > 1_000_000:
        logger.info("Real benchmark file already exists at %s (%d bytes). Skipping download.", dest, dest.stat().st_size)
        return dest

    logger.info("Downloading official ProteinGym GFP benchmark from HuggingFace...")
    logger.info("URL: %s", PROTEINGYM_URL)

    req = urllib.request.Request(
        PROTEINGYM_URL,
        headers={"User-Agent": "VeritasBio/1.0 (Computational Biology Benchmark)"},
    )

    with urllib.request.urlopen(req, timeout=60.0) as resp, open(dest, "wb") as f:
        total_size = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        chunk_size = 64 * 1024
        while True:
            chunk = resp.read(chunk_size)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            if total_size > 0:
                pct = (downloaded / total_size) * 100
                sys.stdout.write(f"\rDownloading: {downloaded / 1024:.1f} KB / {total_size / 1024:.1f} KB ({pct:.1f}%)")
                sys.stdout.flush()

    print()
    logger.info("Successfully downloaded %d bytes to %s", dest.stat().st_size, dest)
    return dest


def analyze_dataset(csv_path: Path = OUTPUT_CSV) -> dict:
    """Parses and computes summary statistics of the authentic DMS dataset."""
    logger.info("Analyzing %s...", csv_path)
    single_mutants = 0
    double_mutants = 0
    higher_mutants = 0
    scores = []

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            mutant = row["mutant"]
            try:
                score = float(row["DMS_score"])
                scores.append(score)
            except ValueError:
                continue

            mut_count = len(mutant.split(":"))
            if mut_count == 1:
                single_mutants += 1
            elif mut_count == 2:
                double_mutants += 1
            else:
                higher_mutants += 1

    total = len(scores)
    scores.sort()
    min_score = scores[0] if scores else 0.0
    max_score = scores[-1] if scores else 0.0
    median_score = scores[total // 2] if scores else 0.0
    top_1_pct_thresh = scores[int(total * 0.99)] if scores else 0.0

    stats = {
        "total_variants": total,
        "single_mutants": single_mutants,
        "double_mutants": double_mutants,
        "higher_order_mutants": higher_mutants,
        "min_score": min_score,
        "max_score": max_score,
        "median_score": median_score,
        "top_1_pct_threshold": top_1_pct_thresh,
    }

    logger.info("==================================================")
    logger.info("PROTEINGYM SARKISYAN 2016 REAL DATA SUMMARY")
    logger.info("==================================================")
    logger.info("Total Empirical Variants: %d", stats["total_variants"])
    logger.info("  - Single Substitutions:  %d", stats["single_mutants"])
    logger.info("  - Double Substitutions:  %d", stats["double_mutants"])
    logger.info("  - Higher-Order Mutants:  %d", stats["higher_order_mutants"])
    logger.info("DMS Score (log fluorescence) Range: [%.2f, %.2f]", stats["min_score"], stats["max_score"])
    logger.info("Median Score:              %.2f", stats["median_score"])
    logger.info("Top 1%% Score Threshold:   %.2f", stats["top_1_pct_threshold"])
    logger.info("==================================================")
    return stats


if __name__ == "__main__":
    download_proteingym_gfp()
    analyze_dataset()

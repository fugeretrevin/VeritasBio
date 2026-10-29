"""
scripts/download_trpb_benchmark.py — Downloads and prepares the authentic TrpB fitness landscape.

Dataset details:
    Target: Tryptophan synthase beta subunit (TrpB), Thermotoga maritima
    Parent: Tm9D8* — a laboratory-evolved TmTrpB carrying 10 substitutions, plus a C-terminal
            LEHHHHHH purification tag. This is the reference the assay was measured against,
            NOT the UniProt TmTrpB wild type (P50909).
    Library: combinatorially complete 4-site saturation at residues 183, 184, 227, 228
             (20^4 = 160,000 variants, every one measured).
    Reference: Johnston, Iqbal et al. (2024), "A combinatorially complete epistatic fitness
               landscape in an enzyme active site", bioRxiv 10.1101/2024.06.23.600144
    Source: https://huggingface.co/datasets/SaProtHub/Dataset-TrpB_fitness_landsacpe

What this script does and does NOT do:
    - The fitness values are copied verbatim from the published `label` column. No rescaling,
      smoothing, imputation or synthetic rows are introduced at any point.
    - The upstream file stores each variant as a full 397-residue sequence. This script
      diffs each sequence against the parent and rewrites it as ProteinGym-style mutation
      notation (e.g. 'V183A:F184I:V227K:S228G') so the Gatekeeper and Quant Engine can
      consume it exactly like the Sarkisyan avGFP benchmark. That is a lossless re-encoding.
    - `DMS_score_bin` is NOT present upstream. It is derived here by a single documented rule:
      1 if the variant's measured fitness is at least the parent's, else 0. See BIN_RULE.
    - The upstream train/valid/test assignment is preserved in a `stage` column so the
      published split is not lost.
"""

from __future__ import annotations

import csv
import json
import logging
import sys
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("download_trpb_benchmark")

ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = Path(__file__).resolve().parent / "data"
WORKER_DATA_DIR = ROOT_DIR / "workers" / "data"

RAW_CSV = DATA_DIR / "TrpB_Johnston_2024_raw.csv"
OUTPUT_CSV = DATA_DIR / "TRPB_THEMA_Johnston_2024.csv"

# Offline metadata cache consumed by workers.rag_worker.fetch_protein_metadata().
PARENT_ACCESSION = "TRPB_TM9D8"
PARENT_CACHE_JSON = WORKER_DATA_DIR / f"{PARENT_ACCESSION}.json"
UNIPROT_ACCESSION = "P50909"  # TmTrpB wild type — provides organism + structural annotations.
UNIPROT_CACHE_JSON = WORKER_DATA_DIR / f"{UNIPROT_ACCESSION}.json"

TRPB_URL = (
    "https://huggingface.co/datasets/SaProtHub/Dataset-TrpB_fitness_landsacpe/"
    "resolve/main/dataset.csv"
)
UNIPROT_URL = f"https://rest.uniprot.org/uniprotkb/{UNIPROT_ACCESSION}.json"

# The four saturated active-site residues (1-indexed; shared with UniProt P50909 numbering).
ACTIVE_SITE_POSITIONS: List[int] = [183, 184, 227, 228]

# Residues carried by the Tm9D8* parent at the four saturated sites (Johnston et al. 2024).
PARENT_RESIDUES: str = "VFVS"

# Tm9D8* was derived from wild-type TmTrpB by directed evolution; the published lineage
# carries exactly this many substitutions. Used as an integrity check, not as an input.
EXPECTED_PARENT_SUBSTITUTIONS: int = 10
EXPECTED_SEQUENCE_LENGTH: int = 397
PURIFICATION_TAG: str = "LEHHHHHH"

BIN_RULE = "DMS_score_bin = 1 if DMS_score >= parent (Tm9D8*) fitness, else 0"

_USER_AGENT = "PanoptesBio/1.0 (Computational Biology Benchmark)"


# -----------------------------------------------------------------------------
# Download
# -----------------------------------------------------------------------------
def _download(url: str, dest: Path, min_bytes: int, force: bool = False) -> Path:
    """Streams `url` to `dest`, skipping the fetch when a plausible file already exists."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force and dest.stat().st_size >= min_bytes:
        logger.info("Already present at %s (%d bytes). Skipping download.", dest, dest.stat().st_size)
        return dest

    logger.info("Downloading %s", url)
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})

    with urllib.request.urlopen(req, timeout=120.0) as resp, open(dest, "wb") as f:
        total_size = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        while True:
            chunk = resp.read(64 * 1024)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            if total_size > 0:
                pct = (downloaded / total_size) * 100
                sys.stdout.write(
                    f"\rDownloading: {downloaded / 1024:.1f} KB / {total_size / 1024:.1f} KB ({pct:.1f}%)"
                )
                sys.stdout.flush()

    if total_size > 0:
        print()
    logger.info("Downloaded %d bytes to %s", dest.stat().st_size, dest)
    return dest


def download_trpb_raw(dest: Path = RAW_CSV, force: bool = False) -> Path:
    """Downloads the published 160,000-variant TrpB landscape CSV from HuggingFace."""
    return _download(TRPB_URL, dest, min_bytes=60_000_000, force=force)


def download_uniprot_trpb(dest: Path = UNIPROT_CACHE_JSON, force: bool = False) -> Path:
    """Downloads the unmodified UniProt record for wild-type TmTrpB."""
    return _download(UNIPROT_URL, dest, min_bytes=10_000, force=force)


# -----------------------------------------------------------------------------
# Parsing & re-encoding
# -----------------------------------------------------------------------------
def _read_raw_rows(csv_path: Path) -> List[Tuple[str, float, str, str]]:
    """
    Reads (sequence, fitness, raw_fitness_text, stage) tuples from the upstream CSV.
    The raw text is carried through so the emitted benchmark reproduces the published
    value digit for digit rather than a float round-trip of it.
    """
    rows: List[Tuple[str, float, str, str]] = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        missing = {"protein", "label"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{csv_path} is missing expected column(s): {sorted(missing)}")
        for row in reader:
            seq = row["protein"].strip().upper()
            raw_score = (row["label"] or "").strip()
            try:
                score = float(raw_score)
            except (TypeError, ValueError):
                continue
            rows.append((seq, score, raw_score, (row.get("stage") or "").strip()))
    return rows


def find_variable_positions(sequences: List[str]) -> List[int]:
    """
    Recovers the saturated positions empirically rather than trusting a constant.
    Returns 1-indexed positions at which the library is not invariant.
    """
    if not sequences:
        raise ValueError("No sequences supplied")
    length = len(sequences[0])
    if any(len(s) != length for s in sequences):
        raise ValueError("Library contains sequences of differing length; alignment is ambiguous")

    variable: List[int] = []
    for idx in range(length):
        first = sequences[0][idx]
        for seq in sequences:
            if seq[idx] != first:
                variable.append(idx + 1)
                break
    return variable


def identify_parent_sequence(
    rows: List[Tuple[str, float, str, str]],
    positions: List[int],
) -> Tuple[str, float]:
    """
    Locates the Tm9D8* parent among the 160,000 members: the single variant carrying the
    published parent residues (PARENT_RESIDUES) at the four saturated sites.
    """
    matches = [
        (seq, score)
        for seq, score, _, _ in rows
        if "".join(seq[p - 1] for p in positions) == PARENT_RESIDUES
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one parent ({PARENT_RESIDUES} at {positions}); found {len(matches)}."
        )
    return matches[0]


def diff_against_uniprot(parent_seq: str, uniprot_seq: str) -> List[str]:
    """
    Lists the substitutions separating the Tm9D8* parent from wild-type TmTrpB, ignoring the
    C-terminal purification tag. Reported so the provenance of the reference is auditable.
    """
    core = parent_seq[: len(uniprot_seq)]
    return [
        f"{wt}{i}{mut}"
        for i, (wt, mut) in enumerate(zip(uniprot_seq, core), 1)
        if wt != mut
    ]


def encode_mutant(seq: str, parent_seq: str, positions: List[int]) -> str:
    """Re-encodes a library member as ProteinGym-style notation relative to the parent."""
    parts = [
        f"{parent_seq[p - 1]}{p}{seq[p - 1]}"
        for p in positions
        if seq[p - 1] != parent_seq[p - 1]
    ]
    return ":".join(parts)


def build_benchmark_csv(
    raw_path: Path = RAW_CSV,
    dest: Path = OUTPUT_CSV,
    uniprot_path: Path = UNIPROT_CACHE_JSON,
) -> Dict[str, object]:
    """
    Converts the upstream sequence-level landscape into a ProteinGym-style substitution CSV
    and writes the matching offline parent-metadata cache. Returns a provenance summary.
    """
    logger.info("Reading %s ...", raw_path)
    rows = _read_raw_rows(raw_path)
    if not rows:
        raise ValueError(f"No usable rows parsed from {raw_path}")
    logger.info("Parsed %d measured library members.", len(rows))

    sequences = [r[0] for r in rows]
    seq_len = len(sequences[0])
    if seq_len != EXPECTED_SEQUENCE_LENGTH:
        logger.warning(
            "Sequence length is %d, expected %d — verify the upstream file before using it.",
            seq_len, EXPECTED_SEQUENCE_LENGTH,
        )

    positions = find_variable_positions(sequences)
    logger.info("Saturated positions recovered from the data: %s", positions)
    if positions != ACTIVE_SITE_POSITIONS:
        raise ValueError(
            f"Saturated positions {positions} do not match the published active site "
            f"{ACTIVE_SITE_POSITIONS}. Refusing to emit a benchmark with uncertain numbering."
        )

    parent_seq, parent_score = identify_parent_sequence(rows, positions)
    logger.info("Parent Tm9D8* (%s at %s) measured fitness: %.9g", PARENT_RESIDUES, positions, parent_score)

    if not parent_seq.endswith(PURIFICATION_TAG):
        logger.warning("Parent sequence does not end in the expected %s tag.", PURIFICATION_TAG)

    # Integrity check: confirm the reference really is the evolved parent, not the WT enzyme.
    uniprot_raw = json.loads(uniprot_path.read_text(encoding="utf-8"))
    uniprot_seq = uniprot_raw["sequence"]["value"].strip().upper()
    parent_subs = diff_against_uniprot(parent_seq, uniprot_seq)
    logger.info(
        "Parent differs from UniProt %s at %d position(s): %s",
        UNIPROT_ACCESSION, len(parent_subs), ", ".join(parent_subs) or "none",
    )
    if len(parent_subs) != EXPECTED_PARENT_SUBSTITUTIONS:
        logger.warning(
            "Expected %d parent substitutions for Tm9D8*, found %d. The reference sequence is "
            "still taken verbatim from the dataset, but double-check the lineage.",
            EXPECTED_PARENT_SUBSTITUTIONS, len(parent_subs),
        )

    dest.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    above_parent = 0
    stage_counts: Counter[str] = Counter()

    with open(dest, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["mutant", "DMS_score", "DMS_score_bin", "stage"])
        for seq, score, raw_score, stage in rows:
            mutant = encode_mutant(seq, parent_seq, positions)
            if not mutant:
                continue  # the parent itself; recorded in the metadata cache instead
            is_above = 1 if score >= parent_score else 0
            above_parent += is_above
            stage_counts[stage] += 1
            writer.writerow([mutant, raw_score, is_above, stage])
            written += 1

    logger.info("Wrote %d variants to %s", written, dest)

    summary: Dict[str, object] = {
        "variants_written": written,
        "parent_score": parent_score,
        "at_or_above_parent": above_parent,
        "positions": positions,
        "parent_sequence": parent_seq,
        "parent_substitutions_vs_uniprot": parent_subs,
        "stage_counts": dict(stage_counts),
    }

    write_parent_metadata(parent_seq, parent_score, parent_subs, uniprot_raw)
    return summary


# -----------------------------------------------------------------------------
# Offline metadata cache for the RAG worker
# -----------------------------------------------------------------------------
def write_parent_metadata(
    parent_seq: str,
    parent_score: float,
    parent_subs: List[str],
    uniprot_raw: Dict,
    dest: Path = PARENT_CACHE_JSON,
) -> Path:
    """
    Writes a UniProt-shaped offline record for the Tm9D8* parent so that
    workers.rag_worker.fetch_protein_metadata('TRPB_TM9D8') resolves without a network call.

    This is a DERIVED record, and says so in its own `panoptesbio_provenance` block:
      - `sequence` is the parent taken verbatim from the Johnston et al. dataset.
      - organism, protein name and structural features are copied unchanged from UniProt
        P50909. Tm9D8* carries only substitutions (no insertions or deletions) across the
        389 residues UniProt covers, so those feature coordinates still apply; features that
        would extend past the UniProt sequence are dropped rather than extrapolated.
      - Nothing here is a measurement, and no measurement is altered anywhere in this file.
    """
    uniprot_len = len(uniprot_raw["sequence"]["value"])
    kept_features = [
        f
        for f in uniprot_raw.get("features", [])
        if (f.get("location", {}).get("end", {}).get("value") or 0) <= uniprot_len
    ]

    record = {
        "primaryAccession": PARENT_ACCESSION,
        "uniProtkbId": "TRPB_TM9D8_JOHNSTON2024",
        "proteinDescription": {
            "recommendedName": {
                "fullName": {
                    "value": "Tryptophan synthase beta chain, engineered parent Tm9D8*"
                }
            }
        },
        "organism": uniprot_raw.get("organism", {}),
        "sequence": {"value": parent_seq, "length": len(parent_seq)},
        "features": kept_features,
        "panoptesbio_provenance": {
            "record_type": "derived — not a UniProt entry",
            "sequence_source": (
                "Tm9D8* parent sequence, verbatim from "
                "huggingface.co/datasets/SaProtHub/Dataset-TrpB_fitness_landsacpe"
            ),
            "annotation_source": f"UniProt {UNIPROT_ACCESSION} (wild-type TmTrpB), copied unchanged",
            "substitutions_vs_uniprot_wild_type": parent_subs,
            "c_terminal_tag": PURIFICATION_TAG,
            "active_site_positions": ACTIVE_SITE_POSITIONS,
            "parent_measured_fitness": parent_score,
            "numbering": (
                f"1-indexed; positions 1-{uniprot_len} share UniProt {UNIPROT_ACCESSION} numbering, "
                f"positions {uniprot_len + 1}-{len(parent_seq)} are the {PURIFICATION_TAG} tag"
            ),
            "reference": (
                "Johnston, Iqbal et al. (2024) bioRxiv 10.1101/2024.06.23.600144 — "
                "A combinatorially complete epistatic fitness landscape in an enzyme active site"
            ),
        },
    }

    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(record, indent=2), encoding="utf-8")
    logger.info("Wrote derived parent metadata cache to %s", dest)
    return dest


# -----------------------------------------------------------------------------
# Reporting
# -----------------------------------------------------------------------------
def analyze_dataset(csv_path: Path = OUTPUT_CSV, parent_score: Optional[float] = None) -> Dict[str, float]:
    """Parses and reports summary statistics of the prepared TrpB benchmark."""
    logger.info("Analyzing %s...", csv_path)
    scores: List[float] = []
    order_counts: Counter[int] = Counter()

    with open(csv_path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                scores.append(float(row["DMS_score"]))
            except (TypeError, ValueError):
                continue
            order_counts[len(row["mutant"].split(":"))] += 1

    scores.sort()
    total = len(scores)
    stats = {
        "total_variants": total,
        "single_mutants": order_counts[1],
        "double_mutants": order_counts[2],
        "triple_mutants": order_counts[3],
        "quadruple_mutants": order_counts[4],
        "min_score": scores[0] if scores else 0.0,
        "max_score": scores[-1] if scores else 0.0,
        "median_score": scores[total // 2] if scores else 0.0,
        "top_1_pct_threshold": scores[int(total * 0.99)] if scores else 0.0,
    }

    logger.info("==================================================")
    logger.info("TRPB JOHNSTON 2024 REAL DATA SUMMARY")
    logger.info("==================================================")
    logger.info("Total Empirical Variants: %d", stats["total_variants"])
    logger.info("  - Single Substitutions:    %d", stats["single_mutants"])
    logger.info("  - Double Substitutions:    %d", stats["double_mutants"])
    logger.info("  - Triple Substitutions:    %d", stats["triple_mutants"])
    logger.info("  - Quadruple Substitutions: %d", stats["quadruple_mutants"])
    logger.info("Fitness Range:             [%.4f, %.4f]", stats["min_score"], stats["max_score"])
    logger.info("Median Fitness:            %.4f", stats["median_score"])
    logger.info("Top 1%% Fitness Threshold: %.4f", stats["top_1_pct_threshold"])
    if parent_score is not None:
        logger.info("Parent Tm9D8* Fitness:     %.4f", parent_score)
    logger.info("Bin rule:                  %s", BIN_RULE)
    logger.info("==================================================")
    return stats


if __name__ == "__main__":
    download_uniprot_trpb()
    download_trpb_raw()
    provenance = build_benchmark_csv()
    analyze_dataset(parent_score=float(provenance["parent_score"]))

"""
scripts/ingest_plate_data.py — Asynchronous Wet-Lab Plate Reader Data Ingestion CLI (Gap 2).

Purpose:
    Ingests physical 96-well and 384-well microplate reader experimental data (fluorescence,
    absorbance, luminescence, or binding Kd) into the PanoptesBio SQLite ledger.
    Connects offline physical biology to the online active learning engine.

Workflow:
    1. Parses robotic plate reader CSV (PanoptesBio anti-evaporation layout or generic plate reader export).
    2. Filters out perimeter hydration barriers, media blanks, and negative controls.
    3. Normalizes candidate signals against the Wild-Type reference parent well (B2).
    4. Applies fail-closed ledger authorization: `ledger_db.authorize_spend(...)`.
    5. Records empirical evidence batch atomically: `ledger_db.record_evidence(...)`.
    6. Returns top discovered hits and allows the next active learning round to train on wet-lab data.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Add repository root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from controller.gatekeeper import Variant
from ledger import database as ledger_db
from workers.rag_worker import fetch_protein_metadata

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("ingest_plate_data")


def parse_plate_csv(
    csv_path: Path,
    normalize_to_wt: bool = True,
) -> Dict[str, Any]:
    """
    Parses plate reader CSV files. Supports PanoptesBio robotic layout and generic plate reader exports.
    """
    if not csv_path.is_file():
        raise FileNotFoundError(f"Plate CSV file not found: {csv_path}")

    rows: List[Dict[str, str]] = []
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = [h.strip() for h in next(reader)]
        header_lower = [h.lower() for h in header]

        for raw_row in reader:
            if not raw_row or not any(raw_row):
                continue
            row_dict = {k: v.strip() for k, v in zip(header_lower, raw_row) if k}
            rows.append(row_dict)

    if not rows:
        raise ValueError(f"No data rows found in plate CSV: {csv_path}")

    # Detect score column
    score_col = None
    for candidate_col in ("score", "fitness", "reading", "value", "signal", "absorbance", "fluorescence", "od600", "rfu"):
        if candidate_col in header_lower:
            score_col = candidate_col
            break

    if not score_col:
        # Fallback: look for any column containing numerical data in row 0
        for k, val in rows[0].items():
            try:
                float(val)
                score_col = k
                break
            except ValueError:
                continue

    if not score_col:
        raise ValueError("Could not detect numerical fitness/reading column in plate CSV.")

    # Locate WT baseline control if present
    wt_reading: Optional[float] = None
    barrier_count = 0
    blank_count = 0
    candidate_records: List[Dict[str, Any]] = []

    for r in rows:
        sample_type = r.get("sample_type", "").lower()
        sample_name = r.get("sample_name", r.get("mutations", r.get("variant", ""))).strip()
        well_id = r.get("well", "").strip()

        # Check for evaporation barriers or blanks
        if sample_type == "evaporation_barrier" or "barrier" in sample_name.lower():
            barrier_count += 1
            continue
        if sample_type == "control_blank" or "blank" in sample_name.lower() or sample_type == "empty":
            blank_count += 1
            continue

        raw_score_str = r.get(score_col, "")
        try:
            raw_score = float(raw_score_str)
        except ValueError:
            continue

        # Detect WT parent control (Well B2 or sample_type control_wt)
        if sample_type == "control_wt" or "wt" in sample_name.lower() or well_id.upper() in ("B2", "B02"):
            if wt_reading is None:
                wt_reading = raw_score

        candidate_records.append({
            "well": well_id,
            "sample_name": sample_name,
            "sample_type": sample_type or "variant",
            "raw_score": raw_score,
            "variance": float(r["variance"]) if "variance" in r and r["variance"] else None,
        })

    # Normalize to WT parent baseline if requested and WT was found
    normalized_records: List[Dict[str, Any]] = []
    for c in candidate_records:
        score = c["raw_score"]
        if normalize_to_wt and wt_reading is not None and wt_reading > 0:
            norm_score = score / wt_reading
        else:
            norm_score = score

        # Filter out negative controls and blanks from final training evidence
        if c["sample_type"] not in ("control_negative", "control_blank"):
            normalized_records.append({
                "well": c["well"],
                "mutations": c["sample_name"],
                "raw_score": c["raw_score"],
                "normalized_score": round(norm_score, 4),
                "variance": c.get("variance"),
            })

    return {
        "score_column": score_col,
        "total_wells": len(rows),
        "barrier_wells_ignored": barrier_count,
        "blank_wells_ignored": blank_count,
        "wt_baseline_score": wt_reading,
        "valid_variants": normalized_records,
    }


def ingest_plate_csv(
    csv_path: Union[str, Path],
    db_path: Optional[Union[str, Path]] = None,
    target_accession: str = "P42212",
    normalize_to_wt: bool = True,
    nominal_cost_usd: float = 0.000100,
) -> Dict[str, Any]:
    """
    Programmatic entrypoint to ingest plate reader data into SQLite ledger.
    """
    csv_file = Path(csv_path)
    parsed = parse_plate_csv(csv_file, normalize_to_wt=normalize_to_wt)
    variants = parsed["valid_variants"]

    if not variants:
        raise ValueError(f"No candidate variants found in {csv_file}")

    # Fetch reference wild-type sequence from UniProt
    metadata = fetch_protein_metadata(target_accession)
    wt_seq = metadata.sequence

    # Construct sequence-score tuples for ledger
    evidence_tuples: List[Tuple[str, float]] = []
    annotated_variants: List[Dict[str, Any]] = []

    for v in variants:
        mut_str = v["mutations"]
        score = v["normalized_score"]
        try:
            if mut_str in ("WT", "WT_Control", f"{metadata.protein_name} WT"):
                seq = wt_seq
            elif ":" in mut_str or any(c.isdigit() for c in mut_str):
                variant_obj = Variant(notation=mut_str)
                seq = variant_obj.apply_to(wt_seq)
            elif len(mut_str) == len(wt_seq):
                seq = mut_str
            else:
                continue

            evidence_tuples.append((seq, score))
            v_annotated = dict(v)
            v_annotated["variant_sequence"] = seq
            annotated_variants.append(v_annotated)
        except Exception as parse_err:
            logger.debug("Skipping unparseable variant notation '%s': %s", mut_str, parse_err)

    if not evidence_tuples:
        raise ValueError("Could not convert any plate entries into valid protein sequences.")

    # Ledger Spend Authorization Handshake
    resolved_db = Path(db_path) if db_path else None
    budget = ledger_db.get_budget(1, db_path=resolved_db)
    if not budget:
        ledger_db.init_db(resolved_db)
        ledger_db.create_budget(total_allocation=5.00, db_path=resolved_db)

    job_id = f"job-wetlab-plate-{int(time.time())}"
    authorized = ledger_db.authorize_spend(
        job_id=job_id,
        target_model="wetlab_plate_reader",
        variant_count=len(evidence_tuples),
        estimated_cost=nominal_cost_usd,
        db_path=resolved_db,
    )

    if not authorized:
        raise PermissionError(f"Ledger rejected wet-lab job {job_id}: budget exhausted or illegal spend.")

    # Record Evidence & Mark Completed
    ledger_db.record_evidence(job_id=job_id, results=evidence_tuples, db_path=resolved_db)
    ledger_db.mark_job_completed(job_id=job_id, db_path=resolved_db)

    # Sort discovered variants by fitness
    annotated_variants.sort(key=lambda x: x["normalized_score"], reverse=True)
    top_hit = annotated_variants[0]

    return {
        "job_id": job_id,
        "csv_path": str(csv_file),
        "target_protein": metadata.protein_name,
        "target_accession": metadata.accession,
        "variants_ingested": len(evidence_tuples),
        "barrier_wells_ignored": parsed["barrier_wells_ignored"],
        "wt_baseline_score": parsed["wt_baseline_score"],
        "top_fitness_discovered": top_hit["normalized_score"],
        "top_variant": top_hit["mutations"],
        "top_variants_ranked": annotated_variants[:5],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="PanoptesBio Wet-Lab Plate Reader Data Ingestion CLI")
    parser.add_argument("--plate-csv", type=str, required=True, help="Path to plate reader CSV")
    parser.add_argument("--db-path", type=str, default=None, help="Path to SQLite ledger database")
    parser.add_argument("--target-accession", type=str, default="P42212", help="UniProt accession (P42212, TRPB_TM9D8)")
    parser.add_argument("--no-normalize-wt", action="store_true", help="Do not normalize scores against WT well")

    args = parser.parse_args()

    result = ingest_plate_csv(
        csv_path=Path(args.plate_csv),
        db_path=Path(args.db_path) if args.db_path else None,
        target_accession=args.target_accession,
        normalize_to_wt=not args.no_normalize_wt,
    )

    print("=" * 60)
    print("PANOPTESBIO WET-LAB PLATE DATA INGESTION SUCCESS")
    print("=" * 60)
    print(f"Job ID:                 {result['job_id']}")
    print(f"Target Protein:         {result['target_protein']} ({result['target_accession']})")
    print(f"Variants Ingested:      {result['variants_ingested']}")
    print(f"Barrier Wells Ignored:  {result['barrier_wells_ignored']}")
    print(f"WT Parent Baseline:     {result['wt_baseline_score']}")
    print(f"Champion Variant:       {result['top_variant']} (Score: {result['top_fitness_discovered']:.4f})")
    print("=" * 60)


if __name__ == "__main__":
    main()

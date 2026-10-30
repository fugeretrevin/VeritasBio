"""
scripts/ingest_plate_data.py — Asynchronous Wet-Lab Plate Reader Data Ingestion CLI (Gap 2 & SPEC-001 Seam 3).

Purpose:
    Ingests physical 96-well and 384-well microplate reader experimental data (fluorescence,
    absorbance, luminescence, or binding Kd) into the VeritasBio SQLite ledger.
    Connects offline physical biology to the online active learning engine.

Workflow:
    1. Parses robotic plate reader CSV (VeritasBio anti-evaporation layout or generic plate reader export).
    2. Filters out perimeter hydration barriers, media blanks, and negative controls.
    3. Applies Grubbs' test (alpha=0.05) and Huber / MAD robust outlier rejection across multi-well replicates.
    4. Normalizes candidate signals against the Wild-Type reference parent well (B2/B02).
    5. Applies fail-closed ledger authorization: `ledger_db.authorize_spend(...)`.
    6. Records empirical evidence batch atomically: `ledger_db.record_evidence(...)`.
    7. Returns top discovered hits and allows the next active learning round to train on wet-lab data.
"""

from __future__ import annotations

import argparse
import csv
import io
import logging
import math
import statistics
import sys
import time
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.optimize import curve_fit

# Add repository root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from controller.gatekeeper import Variant
from ledger import database as ledger_db
from workers.rag_worker import fetch_protein_metadata

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("ingest_plate_data")


def _grubbs_huber_filter_replicates(
    group: List[Dict[str, Any]],
    alpha: float = 0.05,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Filters contaminated or evaporated outlier wells from a multi-well replicate group
    using Grubbs' two-sided test (alpha=0.05) and Huber / MAD robust scale estimation.
    """
    if len(group) < 3:
        return list(group), []

    working = list(group)
    rejected: List[Dict[str, Any]] = []

    try:
        from scipy.stats import t as student_t
    except Exception:
        student_t = None

    while len(working) >= 3:
        vals = [float(item["raw_score"]) for item in working]
        n = len(vals)
        mean_val = sum(vals) / n
        var_val = sum((v - mean_val) ** 2 for v in vals) / max(n - 1, 1)
        std_val = math.sqrt(var_val)
        if std_val <= 1e-12:
            break

        max_idx = max(range(n), key=lambda i: abs(vals[i] - mean_val))
        g_stat = abs(vals[max_idx] - mean_val) / std_val

        if student_t is not None and n > 2:
            t_val = float(student_t.ppf(1.0 - alpha / (2.0 * n), n - 2))
            g_crit = ((n - 1) / math.sqrt(n)) * math.sqrt((t_val * t_val) / (n - 2 + t_val * t_val))
        else:
            g_crit = {3: 1.153, 4: 1.463, 5: 1.672, 6: 1.822}.get(n, 1.96)

        med = statistics.median(vals)
        mad = statistics.median([abs(v - med) for v in vals])
        robust_scale = max(1.4826 * mad, 1e-4 * max(abs(med), 1.0))
        mod_z = abs(vals[max_idx] - med) / robust_scale

        if g_stat > g_crit or mod_z > 4.0:
            outlier_item = dict(working.pop(max_idx))
            outlier_item["grubbs_g"] = round(g_stat, 4)
            outlier_item["modified_z"] = round(mod_z, 4)
            rejected.append(outlier_item)
        else:
            break

    return working, rejected


def parse_plate_csv(
    csv_path: Union[Path, str, io.StringIO],
    normalize_to_wt: bool = True,
    reject_replicate_outliers: bool = True,
) -> Dict[str, Any]:
    """
    Parses plate reader CSV files or raw CSV text. Supports VeritasBio robotic layout and generic plate reader exports.
    """
    close_stream = False
    if isinstance(csv_path, io.StringIO):
        reader_file = csv_path
    elif isinstance(csv_path, Path):
        if not csv_path.is_file():
            raise FileNotFoundError(f"Plate CSV file not found: {csv_path}")
        reader_file = open(csv_path, "r", encoding="utf-8-sig")
        close_stream = True
    elif isinstance(csv_path, str):
        if "\n" in csv_path or not Path(csv_path).is_file():
            reader_file = io.StringIO(csv_path)
        else:
            reader_file = open(Path(csv_path), "r", encoding="utf-8-sig")
            close_stream = True
    else:
        raise ValueError(f"Invalid csv_path type: {type(csv_path)}")

    try:
        reader = csv.reader(reader_file)
        header = [h.strip() for h in next(reader)]
        header_lower = [h.lower() for h in header]

        rows: List[Dict[str, str]] = []
        for raw_row in reader:
            if not raw_row or not any(raw_row):
                continue
            row_dict = {k: v.strip() for k, v in zip(header_lower, raw_row) if k}
            rows.append(row_dict)
    finally:
        if close_stream:
            reader_file.close()

    if not rows:
        raise ValueError("No data rows found in plate CSV.")

    # Detect score column
    score_col = None
    for candidate_col in ("score", "fitness", "reading", "value", "signal", "absorbance", "fluorescence", "od600", "rfu"):
        if candidate_col in header_lower:
            score_col = candidate_col
            break

    if not score_col:
        for k, val in rows[0].items():
            try:
                float(val)
                score_col = k
                break
            except ValueError:
                continue

    if not score_col:
        raise ValueError("Could not detect numerical fitness/reading column in plate CSV.")

    wt_reading: Optional[float] = None
    barrier_count = 0
    blank_count = 0
    candidate_records: List[Dict[str, Any]] = []

    for r in rows:
        sample_type = r.get("sample_type", "").lower()
        sample_name = r.get("sample_name", r.get("mutations", r.get("variant", ""))).strip()
        well_id = r.get("well", "").strip()

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

    outliers_rejected: List[Dict[str, Any]] = []
    filtered_candidates: List[Dict[str, Any]] = []

    if reject_replicate_outliers:
        by_sample: Dict[str, List[Dict[str, Any]]] = {}
        order_keys: List[str] = []
        for c in candidate_records:
            if c["sample_type"] in ("control_negative", "control_blank"):
                continue
            key = c["sample_name"] or f"__well_{c['well']}"
            if key not in by_sample:
                by_sample[key] = []
                order_keys.append(key)
            by_sample[key].append(c)

        for key in order_keys:
            grp = by_sample[key]
            clean_grp, rej_grp = _grubbs_huber_filter_replicates(grp)
            outliers_rejected.extend(rej_grp)
            if len(clean_grp) == 1:
                rec = dict(clean_grp[0])
                rec["standard_error"] = math.sqrt(rec["variance"]) if rec.get("variance") else 0.0
                filtered_candidates.append(rec)
            elif len(clean_grp) > 1:
                raw_vals = [float(item["raw_score"]) for item in clean_grp]
                mean_raw = sum(raw_vals) / len(raw_vals)
                norm_denom = wt_reading if (normalize_to_wt and wt_reading is not None and wt_reading > 0) else 1.0
                norm_vals = [v / norm_denom for v in raw_vals]
                mean_norm = sum(norm_vals) / len(norm_vals)
                emp_var = sum((nv - mean_norm) ** 2 for nv in norm_vals) / max(len(norm_vals) - 1, 1)
                std_err = math.sqrt(emp_var / len(norm_vals))
                rep_rec = dict(clean_grp[0])
                rep_rec["raw_score"] = round(mean_raw, 4)
                rep_rec["variance"] = round(emp_var, 6)
                rep_rec["standard_error"] = round(std_err, 6)
                rep_rec["replicate_count"] = len(clean_grp)
                filtered_candidates.append(rep_rec)
    else:
        filtered_candidates = [
            c for c in candidate_records if c["sample_type"] not in ("control_negative", "control_blank")
        ]

    normalized_records: List[Dict[str, Any]] = []
    for c in filtered_candidates:
        score = c["raw_score"]
        if normalize_to_wt and wt_reading is not None and wt_reading > 0:
            norm_score = score / wt_reading
        else:
            norm_score = score

        if c["sample_type"] not in ("control_negative", "control_blank"):
            normalized_records.append({
                "well": c["well"],
                "mutations": c["sample_name"],
                "raw_score": c["raw_score"],
                "normalized_score": round(norm_score, 4),
                "variance": c.get("variance"),
                "standard_error": c.get("standard_error"),
                "replicate_count": c.get("replicate_count", 1),
            })

    return {
        "score_column": score_col,
        "total_wells": len(rows),
        "barrier_wells_ignored": barrier_count,
        "blank_wells_ignored": blank_count,
        "wt_baseline_score": wt_reading,
        "outliers_rejected": outliers_rejected,
        "valid_variants": normalized_records,
    }


def ingest_plate_csv(
    csv_path: Optional[Union[str, Path, io.StringIO]] = None,
    db_path: Optional[Union[str, Path]] = None,
    target_accession: str = "P42212",
    normalize_to_wt: bool = True,
    nominal_cost_usd: float = 0.000100,
    operator_uid: Optional[str] = None,
    operator_email: Optional[str] = None,
    csv_input: Optional[Union[str, Path, io.StringIO]] = None,
) -> Dict[str, Any]:
    """
    Programmatic entrypoint to ingest plate reader data into SQLite ledger.
    """
    target_input = csv_path if csv_path is not None else csv_input
    if target_input is None:
        raise ValueError("Must provide csv_path or csv_input")

    parsed = parse_plate_csv(target_input, normalize_to_wt=normalize_to_wt)
    variants = parsed["valid_variants"]

    if not variants:
        raise ValueError("No candidate variants found in plate input.")

    metadata = fetch_protein_metadata(target_accession)
    wt_seq = metadata.sequence

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

    resolved_db = Path(db_path) if db_path else None
    budget = ledger_db.get_budget(1, db_path=resolved_db)
    if not budget:
        ledger_db.init_db(resolved_db)
        ledger_db.create_budget(total_allocation=5.00, db_path=resolved_db)

    job_id = f"job-wetlab-plate-{int(time.time())}-{uuid.uuid4().hex[:6]}"
    authorized = ledger_db.authorize_spend(
        job_id=job_id,
        target_model="wetlab_plate_reader",
        variant_count=len(evidence_tuples),
        estimated_cost=nominal_cost_usd,
        db_path=resolved_db,
        operator_uid=operator_uid,
        operator_email=operator_email,
    )

    if not authorized:
        raise PermissionError(f"Ledger rejected wet-lab job {job_id}: budget exhausted or illegal spend.")

    ledger_db.record_evidence(job_id=job_id, results=evidence_tuples, db_path=resolved_db)
    ledger_db.mark_job_completed(job_id=job_id, db_path=resolved_db)

    annotated_variants.sort(key=lambda x: x["normalized_score"], reverse=True)
    top_hit = annotated_variants[0]

    audit_blocks = ledger_db.list_audit_blocks(db_path=resolved_db)
    latest_block = audit_blocks[-1] if audit_blocks else {}

    csv_desc = str(csv_input) if isinstance(csv_input, (str, Path)) and "\n" not in str(csv_input) else "in-memory / uploaded microplate"

    return {
        "job_id": job_id,
        "csv_path": csv_desc,
        "target_protein": metadata.protein_name,
        "target_accession": metadata.accession,
        "variants_ingested": len(evidence_tuples),
        "barrier_wells_ignored": parsed["barrier_wells_ignored"],
        "outliers_rejected": parsed.get("outliers_rejected", []),
        "wt_baseline_score": parsed["wt_baseline_score"],
        "top_fitness_discovered": top_hit["normalized_score"],
        "top_variant": top_hit["mutations"],
        "top_variants_ranked": annotated_variants[:10],
        "block_hash": latest_block.get("block_hash"),
        "block_index": latest_block.get("block_index"),
        "timestamp": latest_block.get("timestamp"),
    }


def fit_dose_response_4pl(
    concentrations: Sequence[float],
    responses: Sequence[float],
) -> Dict[str, float]:
    """
    Fits a 4-Parameter Logistic (Hill) dose-response curve (SPEC-002 Seam 4):
      y(c) = bottom + (top - bottom) / (1 + (c / ic50)^hill_slope)
    Returns IC50/EC50/Kd, Hill slope, asymptotes, R^2, and propagated potency variance
    for the Heteroscedastic Gaussian Process.
    """
    c_arr = np.asarray(concentrations, dtype=np.float64)
    y_arr = np.asarray(responses, dtype=np.float64)
    if len(c_arr) < 4 or len(c_arr) != len(y_arr):
        raise ValueError("fit_dose_response_4pl requires at least 4 paired concentration-response points")

    valid = c_arr > 0
    c_pos = c_arr[valid]
    y_pos = y_arr[valid]

    bottom_init = float(np.min(y_pos))
    top_init = float(np.max(y_pos))
    mid_y = 0.5 * (top_init + bottom_init)
    mid_idx = int(np.argmin(np.abs(y_pos - mid_y)))
    ic50_init = max(float(c_pos[mid_idx]), 1e-6)
    hill_init = 1.0 if y_pos[0] >= y_pos[-1] else -1.0

    def _4pl(c: np.ndarray, bottom: float, top: float, ic50: float, hill: float) -> np.ndarray:
        ic50_safe = np.maximum(ic50, 1e-12)
        return bottom + (top - bottom) / (1.0 + np.power(c / ic50_safe, hill))

    try:
        popt, pcov = curve_fit(
            _4pl,
            c_pos,
            y_pos,
            p0=[bottom_init, top_init, ic50_init, hill_init],
            bounds=(
                [-np.inf, -np.inf, 1e-9, -5.0],
                [np.inf, np.inf, np.inf, 5.0],
            ),
            maxfev=5000,
        )
        bottom_fit, top_fit, ic50_fit, hill_fit = [float(x) for x in popt]
        ic50_var = float(pcov[2, 2]) if (pcov is not None and np.isfinite(pcov[2, 2]) and pcov[2, 2] > 0) else (0.05 * ic50_fit) ** 2
    except Exception:
        bottom_fit, top_fit, ic50_fit, hill_fit = bottom_init, top_init, ic50_init, hill_init
        ic50_var = (0.10 * ic50_fit) ** 2

    y_pred = _4pl(c_pos, bottom_fit, top_fit, ic50_fit, hill_fit)
    ss_res = float(np.sum((y_pos - y_pred) ** 2))
    ss_tot = float(np.sum((y_pos - np.mean(y_pos)) ** 2))
    r_sq = 1.0 - (ss_res / max(ss_tot, 1e-12))

    rel_se = math.sqrt(max(ic50_var, 1e-12)) / max(ic50_fit, 1e-9)
    log10_se = rel_se / math.log(10.0)
    potency_var = max(rel_se ** 2, 1e-5)

    return {
        "ic50": round(ic50_fit, 6),
        "hill_slope": round(abs(hill_fit), 4),
        "top_asymptote": round(top_fit, 4),
        "bottom_asymptote": round(bottom_fit, 4),
        "r_squared": round(r_sq, 5),
        "log10_ic50_se": round(log10_se, 6),
        "potency_variance": round(potency_var, 6),
    }


def fit_michaelis_menten_kinetics(
    timepoints_sec: Sequence[float],
    absorbances: Sequence[float],
    extinction_coeff_m_cm: float = 6220.0,
    pathlength_cm: float = 0.5,
) -> Dict[str, float]:
    """
    Extracts initial catalytic reaction velocity v0 (Delta Abs / sec and uM / sec)
    from time-series absorbance progression curves while automatically trimming
    late-stage substrate-depletion plateau points (SPEC-002 Seam 4).
    """
    t_arr = np.asarray(timepoints_sec, dtype=np.float64)
    a_arr = np.asarray(absorbances, dtype=np.float64)
    if len(t_arr) < 3 or len(t_arr) != len(a_arr):
        raise ValueError("fit_michaelis_menten_kinetics requires at least 3 paired time-absorbance points")

    # Identify linear initial rate regime by trimming late points where instantaneous slope drops < 50% of initial slope
    dt0 = max(float(t_arr[1] - t_arr[0]), 1e-6)
    init_slope = float(a_arr[1] - a_arr[0]) / dt0
    k_end = len(t_arr)
    for k in range(2, len(t_arr)):
        dt_k = max(float(t_arr[k] - t_arr[k - 1]), 1e-6)
        inst_slope = float(a_arr[k] - a_arr[k - 1]) / dt_k
        if abs(init_slope) > 1e-9 and abs(inst_slope) < 0.45 * abs(init_slope):
            k_end = k
            break

    t_lin = t_arr[:k_end]
    a_lin = a_arr[:k_end]
    n_lin = len(t_lin)

    t_mean = float(np.mean(t_lin))
    a_mean = float(np.mean(a_lin))
    denom = float(np.sum((t_lin - t_mean) ** 2))
    slope = float(np.sum((t_lin - t_mean) * (a_lin - a_mean))) / max(denom, 1e-12)
    intercept = a_mean - slope * t_mean

    a_pred = slope * t_lin + intercept
    ss_res = float(np.sum((a_lin - a_pred) ** 2))
    ss_tot = float(np.sum((a_lin - a_mean) ** 2))
    r_sq = 1.0 - (ss_res / max(ss_tot, 1e-12))
    slope_var = (ss_res / max(n_lin - 2, 1)) / max(denom, 1e-12) if n_lin > 2 else 1e-6

    # Beer-Lambert conversion: dC/dt (uM/s) = (dAbs/dt) / (extinction_coeff * pathlength) * 1e6
    v0_um_per_sec = (slope / max(extinction_coeff_m_cm * pathlength_cm, 1e-6)) * 1e6

    return {
        "initial_slope_abs_per_sec": round(slope, 6),
        "v0_um_per_sec": round(v0_um_per_sec, 6),
        "r_squared": round(r_sq, 5),
        "linear_points_used": n_lin,
        "v0_variance": round(max(slope_var, 1e-8), 8),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="VeritasBio Wet-Lab Plate Reader Data Ingestion CLI")
    parser.add_argument("--plate-csv", type=str, required=True, help="Path to plate reader CSV")
    parser.add_argument("--db-path", type=str, default=None, help="Path to SQLite ledger database")
    parser.add_argument("--target-accession", type=str, default="P42212", help="UniProt accession (P42212, TRPB_TM9D8)")
    parser.add_argument("--no-normalize-wt", action="store_true", help="Do not normalize scores against WT well")

    args = parser.parse_args()

    result = ingest_plate_csv(
        csv_input=Path(args.plate_csv),
        db_path=Path(args.db_path) if args.db_path else None,
        target_accession=args.target_accession,
        normalize_to_wt=not args.no_normalize_wt,
    )

    print("=" * 60)
    print("VERITASBIO WET-LAB PLATE DATA INGESTION SUCCESS")
    print("=" * 60)
    print(f"Job ID:                 {result['job_id']}")
    print(f"Target Protein:         {result['target_protein']} ({result['target_accession']})")
    print(f"Variants Ingested:      {result['variants_ingested']}")
    print(f"Barrier Wells Ignored:  {result['barrier_wells_ignored']}")
    print(f"WT Parent Baseline:     {result['wt_baseline_score']}")
    print(f"Champion Variant:       {result['top_variant']} (Score: {result['top_fitness_discovered']:.4f})")
    if result.get("block_hash"):
        print(f"Audit Block Hash:       {result['block_hash']}")
    print("=" * 60)


if __name__ == "__main__":
    main()

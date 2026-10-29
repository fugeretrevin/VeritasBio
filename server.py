"""
server.py — FastAPI High-Performance Backend Bridge for PanoptesBio / VeritasBio.

Exposes REST and streaming APIs for:
  - Ledger & Budget Telemetry (Micro-USD double-entry accounting)
  - FDA 21 CFR Part 11 Cryptographic Merkle Hash Chain & Electronic Signatures
  - Firebase Authentication Token Verification
  - 3D Protein Structures & Contact Maps (PDB)
  - Active Learning Campaign Orchestration (TrpB, GB1, Trastuzumab, Scorpion Antivenom, avGFP)
  - Autonomous De Novo Binder Diffusion & Inverse Folding Studio (RFdiffusion + ProteinMPNN)
  - Wet-Lab Robotics & Synthesis Exports (Twist Bioscience CSV, Opentrons OT-2, Benchling LIMS)
"""

from __future__ import annotations

import json
import logging
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException, Query, Request, status
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from controller.auth import (
    AuthenticatedUser,
    create_part11_electronic_signature,
    verify_firebase_id_token,
)
from controller.firebase_config import get_firebase_config
from controller.gatekeeper import (
    BatchEvaluationRequest,
    DeNovoDiffusionRequest,
    Gatekeeper,
    SingleMutation,
    Variant,
)
from controller.orchestrator import Orchestrator, OrchestratorConfig
from ledger import database as ledger_db
from workers.diffusion_worker import DeNovoBinderRequest, generate_budgeted_de_novo_binders
from workers.mock_workers import MockLLMWorker
from workers.rag_worker import fetch_protein_metadata
from scripts.ingest_plate_data import ingest_plate_csv
import numpy as np
from workers.math_worker import compute_multi_objective_pareto_front, compute_biological_prior
from workers.lims_worker import BenchlingClient
from workers.wetlab_worker import (
    generate_cryptographic_audit_certificate,
    generate_opentrons_protocol,
    generate_plate_layout,
    generate_twist_order_csv,
)

logger = logging.getLogger("panoptesbio.server")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")

app = FastAPI(
    title="PanoptesBio Cloud Economics Engine",
    description="Deterministic Cloud Economics & Biologics Evolution Backend API",
    version="3.2.0",
)

# Enable CORS for local Vite development server
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

DATA_DIR = Path(__file__).resolve().parent / "workers" / "data"
DEFAULT_DB_PATH = Path(__file__).resolve().parent / "ledger" / "ledger.db"

# Ensure ledger is initialized
ledger_db.init_db(DEFAULT_DB_PATH)
if ledger_db.get_budget(1, DEFAULT_DB_PATH) is None:
    ledger_db.create_budget(5.00, budget_id=1, db_path=DEFAULT_DB_PATH)


# -----------------------------------------------------------------------------
# Request & Response Schemas
# -----------------------------------------------------------------------------
class VerifyTokenRequest(BaseModel):
    id_token: str


class SignActionRequest(BaseModel):
    id_token: Optional[str] = "panoptes_demo_token"
    action: str = "Approval of Modal A10G Compute Spend"
    meaning: str = "Authorized compute budget spend under 21 CFR §11.200"


class LaunchCampaignRequest(BaseModel):
    target_accession: str = "TRPB_TM9D8"
    total_budget_usd: float = 1.00
    max_iterations: int = 4
    batch_size: int = 15
    target_model: str = "mock_esm"
    enable_pocket_hopping: bool = True
    auth_token: Optional[str] = "panoptes_demo_token"


class ParetoFrontierRequest(BaseModel):
    target_protein: str = "avGFP"
    wt_sequence: Optional[str] = None
    variants: Optional[List[str]] = None
    weights: Optional[Dict[str, float]] = None
    active_site_positions: Optional[List[int]] = None
    auth_token: Optional[str] = "panoptes_demo_token"


class DeNovoRequestPayload(BaseModel):
    target_scaffold: str = "SCFV_9004G_CN2"
    num_designs: int = 3
    binder_length: int = 65
    budget_usd: float = 0.50
    auth_token: Optional[str] = "panoptes_demo_token"


# -----------------------------------------------------------------------------
# API Endpoints
# -----------------------------------------------------------------------------
@app.get("/api/health")
def get_health() -> Dict[str, Any]:
    """Returns engine health, budget status, and Firebase configuration status."""
    rem = ledger_db.get_remaining_budget(1, DEFAULT_DB_PATH)
    budget = ledger_db.get_budget(1, DEFAULT_DB_PATH)
    fb_cfg = get_firebase_config()
    return {
        "status": "online",
        "engine": "PanoptesBio / VeritasBio v3.2",
        "database": "sqlite_immutable_merkle",
        "budget": {
            "total_allocation_usd": budget["total_allocation"] if budget else 5.0,
            "cumulative_spent_usd": budget["cumulative_spent"] if budget else 0.0,
            "remaining_usd": rem,
        },
        "firebase": {
            "projectId": fb_cfg.project_id,
            "authDomain": fb_cfg.auth_domain,
            "storageBucket": fb_cfg.storage_bucket,
        },
    }


@app.get("/api/targets")
def list_targets() -> List[Dict[str, Any]]:
    """Returns the list of validated empirical benchmarks and therapeutic targets."""
    return [
        {
            "id": "TRPB_TM9D8",
            "name": "Tryptophan Synthase Beta (TrpB)",
            "benchmark": "Johnston et al. Science 2024 (160K Combinatorial DMS)",
            "pdb_id": "8EH0",
            "description": "Catalytic communication loop engineering gating substrate access to the PLP cofactor.",
            "hotspots": [183, 184, 227, 228],
        },
        {
            "id": "SCFV_9004G_CN2",
            "name": "Scorpion Antivenom scFv 9004G : Cn2 Neurotoxin",
            "benchmark": "Centruroides noxius Neutralizing Campaign",
            "pdb_id": "2YBR",
            "description": "Humanized antibody neutralizing scorpion neurotoxin with canonical disulfide preservation.",
            "hotspots": [98, 99, 100, 101, 102],
        },
        {
            "id": "HER2_TRASTUZUMAB",
            "name": "Trastuzumab Antibody : HER2 Epitope",
            "benchmark": "Therapeutic Biologics Paratope Landscape (38K DMS)",
            "pdb_id": "1N8Z",
            "description": "CDR-H3 paratope optimization for oncology receptor affinity.",
            "hotspots": [99, 100, 101, 102],
        },
        {
            "id": "SPG1_STRSG",
            "name": "Protein G (GB1)",
            "benchmark": "Wu et al. eLife 2016 (150K Epistatic DMS)",
            "pdb_id": "2QMT",
            "description": "Thermodynamic folding stability and core beta-sheet packing.",
            "hotspots": [39, 40, 41, 54],
        },
        {
            "id": "P42212",
            "name": "Aequorea victoria GFP (avGFP)",
            "benchmark": "Sarkisyan et al. Nature 2016 (54K Nature DMS)",
            "pdb_id": "1EMA",
            "description": "Internal chromophore triad maturation and pi-pi stacking.",
            "hotspots": [64, 65, 68, 203],
        },
    ]


@app.get("/api/structure/{target_id}")
def get_structure_pdb(target_id: str) -> Dict[str, Any]:
    """Returns PDB crystallographic coordinates for the requested target or de novo design."""
    clean_id = target_id.upper().strip()

    pdb_filename = "8EH0.pdb"
    if "2YBR" in clean_id or "SCORPION" in clean_id or "CN2" in clean_id:
        pdb_filename = "2YBR.pdb"
    elif "1N8Z" in clean_id or "HER2" in clean_id or "TRASTUZUMAB" in clean_id:
        pdb_filename = "1N8Z.pdb"
    elif "1EMA" in clean_id or "GFP" in clean_id or "P42212" in clean_id:
        pdb_filename = "1EMA.pdb"
    elif "2QMT" in clean_id or "GB1" in clean_id or "SPG1" in clean_id:
        pdb_filename = "2QMT.pdb"

    pdb_path = DATA_DIR / pdb_filename
    if not pdb_path.is_file():
        raise HTTPException(status_code=404, detail=f"Structure PDB {pdb_filename} not found")

    content = pdb_path.read_text(encoding="utf-8")
    return {
        "target_id": target_id,
        "pdb_filename": pdb_filename,
        "pdb_content": content,
    }


@app.get("/api/ledger")
def get_ledger_state() -> Dict[str, Any]:
    """Returns the complete double-entry ledger state and cryptographic Merkle audit blocks."""
    budget = ledger_db.get_budget(1, db_path=DEFAULT_DB_PATH)
    jobs = ledger_db.list_jobs(db_path=DEFAULT_DB_PATH)
    evidence = ledger_db.get_evidence(db_path=DEFAULT_DB_PATH)
    audit_blocks = ledger_db.list_audit_blocks(db_path=DEFAULT_DB_PATH)
    is_valid, count, msg = ledger_db.verify_audit_integrity(db_path=DEFAULT_DB_PATH)

    # Sort evidence by score descending
    sorted_ev = sorted(evidence, key=lambda e: float(e.get("computed_score", 0.0)), reverse=True)[:50]

    return {
        "budget": budget,
        "jobs_count": len(jobs),
        "jobs": list(reversed(jobs))[:30],
        "top_variants": sorted_ev[:15],
        "audit_integrity": {
            "is_valid": is_valid,
            "block_count": count,
            "verification_message": msg,
        },
        "audit_blocks": list(reversed(audit_blocks))[:30],
    }


@app.post("/api/auth/verify")
def verify_token(req: VerifyTokenRequest) -> Dict[str, Any]:
    """Verifies a Firebase ID token and returns authenticated operator identity and RBAC role."""
    try:
        user = verify_firebase_id_token(req.id_token)
        return {"authenticated": True, "user": user.model_dump()}
    except Exception as exc:
        raise HTTPException(status_code=401, detail=f"Authentication failed: {exc}") from exc


@app.post("/api/auth/sign")
def sign_part11_action(req: SignActionRequest) -> Dict[str, Any]:
    """Generates an FDA 21 CFR §11.200 HMAC-SHA256 compliant electronic signature."""
    try:
        user = verify_firebase_id_token(req.id_token)
        sig = create_part11_electronic_signature(user, req.action, req.meaning)
        return {"success": True, "signature": sig}
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/diffusion/generate")
def generate_diffusion_binders(req: DeNovoRequestPayload) -> Dict[str, Any]:
    """Generates budgeted de novo binders docked to the target epitope via RFdiffusion + ProteinMPNN."""
    try:
        # Authenticate operator if token provided
        op_uid = None
        op_email = None
        if req.auth_token:
            user = verify_firebase_id_token(req.auth_token)
            op_uid = user.uid
            op_email = user.email

        binder_req = DeNovoBinderRequest(
            target_scaffold=req.target_scaffold,
            num_designs=req.num_designs,
            binder_length=req.binder_length,
            budget_usd=req.budget_usd,
            budget_id=1,
            operator_uid=op_uid,
            operator_email=op_email,
        )

        result = generate_budgeted_de_novo_binders(binder_req, db_path=DEFAULT_DB_PATH)
        return {
            "success": result.success,
            "job_id": result.job_id,
            "designs_count": len(result.designs),
            "estimated_cost_usd": result.estimated_cost_usd,
            "designs": [d.model_dump() for d in result.designs],
            "audit_block_hash": result.audit_block_hash,
            "message": result.message,
        }
    except Exception as exc:
        logger.exception("Diffusion generation error: %s", exc)
        raise HTTPException(status_code=500, detail=str(exc)) from exc



@app.get("/api/campaign/stream")
def stream_campaign(
    target_accession: str = Query("TRPB_TM9D8"),
    total_budget_usd: float = Query(1.00),
    max_iterations: int = Query(4),
    batch_size: int = Query(15),
    target_model: str = Query("mock_esm"),
    enable_pocket_hopping: bool = Query(True),
    auth_token: Optional[str] = Query("panoptes_demo_token"),
):
    """
    Server-Sent Events (SSE) streaming endpoint for active learning campaigns.
    Yields live round-by-round OrchestratorEvents to the React frontend in real-time.
    """
    def event_generator():
        try:
            config = OrchestratorConfig(
                target_accession=target_accession,
                total_budget_usd=total_budget_usd,
                max_iterations=max_iterations,
                batch_size=batch_size,
                target_model=target_model,
                use_mock_llm=True,
                enable_live_llm=False,
                enable_regret_stopping=True,
                enable_pocket_hopping=enable_pocket_hopping,
                db_path=str(DEFAULT_DB_PATH),
            )
            orch = Orchestrator(config)
            for event in orch.run_stream():
                payload = {
                    "iteration": event.iteration,
                    "step": event.step,
                    "message": event.message,
                    "remaining_budget_usd": event.remaining_budget_usd,
                    "cumulative_spent_usd": event.cumulative_spent_usd,
                    "top_variants": event.top_variants if event.top_variants else None,
                    "decision": event.decision if event.decision else None,
                }
                yield f"data: {json.dumps(payload)}\n\n"

            # Final complete event with top discovered evidence
            evidence = ledger_db.get_evidence(db_path=DEFAULT_DB_PATH)
            sorted_ev = sorted(evidence, key=lambda e: float(e.get("computed_score", 0.0)), reverse=True)[:15]
            final_payload = {
                "step": "complete",
                "message": "Optimization campaign completed.",
                "top_variants": sorted_ev,
            }
            yield f"data: {json.dumps(final_payload)}\n\n"
        except Exception as exc:
            logger.exception("Streaming campaign failed: %s", exc)
            err_payload = {"step": "error", "message": str(exc)}
            yield f"data: {json.dumps(err_payload)}\n\n"

    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no",
        },
    )

@app.post("/api/campaign/run")
def run_campaign(req: LaunchCampaignRequest) -> Dict[str, Any]:
    """Executes an active learning optimization run with Bayesian regret stopping & pocket-hopping."""
    try:
        config = OrchestratorConfig(
            target_accession=req.target_accession,
            total_budget_usd=req.total_budget_usd,
            max_iterations=req.max_iterations,
            batch_size=req.batch_size,
            target_model=req.target_model,
            use_mock_llm=True,
            enable_live_llm=False,
            enable_regret_stopping=True,
            enable_pocket_hopping=req.enable_pocket_hopping,
            db_path=str(DEFAULT_DB_PATH),
        )

        orch = Orchestrator(config)
        events = []
        for event in orch.run_stream():
            events.append({
                "iteration": event.iteration,
                "step": event.step,
                "message": event.message,
                "remaining_budget_usd": event.remaining_budget_usd,
                "cumulative_spent_usd": event.cumulative_spent_usd,
            })

        evidence = ledger_db.get_evidence(db_path=DEFAULT_DB_PATH)
        sorted_ev = sorted(evidence, key=lambda e: float(e.get("computed_score", 0.0)), reverse=True)[:10]

        return {
            "success": True,
            "iterations_completed": events[-1]["iteration"] if events else 0,
            "total_events": len(events),
            "events": events[-15:],  # Return latest 15 events
            "top_variants": sorted_ev,
        }
    except Exception as exc:
        logger.exception("Campaign run failed: %s", exc)
        raise HTTPException(status_code=500, detail=str(exc)) from exc


@app.post("/api/wetlab/export")
def export_wetlab_manifests(
    target_name: str = "PanoptesBio_Leads",
    plate_format: str = "96_well",
) -> Dict[str, Any]:
    """Generates physical biology translation files: Twist CSV, Opentrons OT-2 script, FDA certificate."""
    evidence = ledger_db.get_evidence(db_path=DEFAULT_DB_PATH)
    sorted_ev = sorted(evidence, key=lambda e: float(e.get("computed_score", 0.0)), reverse=True)[:56]

    twist_csv = generate_twist_order_csv(sorted_ev, protein_name=target_name)
    var_names = [e.get("variant_sequence", "")[:12] for e in sorted_ev]
    opentrons_py = generate_opentrons_protocol(var_names, target_protein=target_name)
    fda_cert = generate_cryptographic_audit_certificate(db_path=DEFAULT_DB_PATH, target_protein=target_name)

    # Generate plate layout and Benchling LIMS manifest
    plate_layout = generate_plate_layout(sorted_ev, plate_format=plate_format, anti_edge_barrier=True)
    benchling_client = BenchlingClient(mock_mode=True)
    benchling_csv = benchling_client.export_lims_manifest(
        plate_name=f"{target_name}_PLT01",
        plate_map=plate_layout,
    )

    return {
        "twist_csv": twist_csv,
        "opentrons_script": opentrons_py,
        "fda_certificate_json": json.dumps(fda_cert, indent=2),
        "benchling_manifest": benchling_csv,
        "plate_csv": plate_layout.get("csv_layout", ""),
    }


@app.post("/api/campaign/ingest-plate")
async def ingest_plate_endpoint(request: Request) -> Dict[str, Any]:
    """
    Ingests physical 96-well / 384-well microplate reader CSV experimental data
    into the VeritasBio SQLite ledger and Merkle audit trail.
    Supports:
      - Multipart form-data with file upload ('file')
      - Form fields ('csv_content', 'target_accession', 'normalize_to_wt', 'auth_token')
      - JSON body ({"csv_content": "...", "target_accession": "...", "normalize_to_wt": True})
    """
    try:
        content_type = request.headers.get("content-type", "")
        csv_text = ""
        target_accession = "P42212"
        normalize_to_wt = True
        auth_token = "panoptes_demo_token"

        if "application/json" in content_type:
            body = await request.json()
            csv_text = body.get("csv_content", "")
            target_accession = body.get("target_accession", "P42212")
            normalize_to_wt = bool(body.get("normalize_to_wt", True))
            auth_token = body.get("auth_token", "panoptes_demo_token")
        elif "multipart/form-data" in content_type or "application/x-www-form-urlencoded" in content_type:
            form = await request.form()
            if "file" in form and hasattr(form["file"], "read"):
                uploaded = form["file"]
                raw_bytes = await uploaded.read()
                csv_text = raw_bytes.decode("utf-8-sig", errors="replace")
            elif "csv_content" in form:
                csv_text = str(form["csv_content"])
            target_accession = str(form.get("target_accession", "P42212"))
            norm_val = form.get("normalize_to_wt", "true")
            normalize_to_wt = str(norm_val).lower() in ("true", "1", "yes")
            auth_token = str(form.get("auth_token", "panoptes_demo_token"))
        else:
            raw_bytes = await request.body()
            csv_text = raw_bytes.decode("utf-8-sig", errors="replace")

        if not csv_text.strip():
            raise HTTPException(status_code=400, detail="Microplate CSV content or file is empty.")

        # Authenticate operator if token provided
        user = verify_firebase_id_token(auth_token) if auth_token else None
        operator_uid = user.uid if user else "usr_panoptes_operator"
        operator_email = user.email if user else "scientist@panoptesbio.com"

        res = ingest_plate_csv(
            csv_input=csv_text,
            db_path=DEFAULT_DB_PATH,
            target_accession=target_accession,
            normalize_to_wt=normalize_to_wt,
            operator_uid=operator_uid,
            operator_email=operator_email,
        )
        return {"success": True, **res}
    except HTTPException:
        raise
    except Exception as exc:
        logger.exception("Microplate data ingestion failed: %s", exc)
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/campaign/pareto-frontier")
def get_pareto_frontier(request: ParetoFrontierRequest) -> Dict[str, Any]:
    """
    Computes true multi-objective Pareto trade-off surface across 4 dimensions:
      1. Predicted functional fitness (to maximize)
      2. Zero-shot conformational stability prior (to maximize)
      3. Wet-lab synthesizability & expression score (to maximize)
      4. Biologics developability & CMC safety score (to maximize)
    Returns non-dominated Rank 0 leads and pre-formatted chart data for interactive UI exploration.
    """
    target_clean = request.target_protein.upper().strip()
    wt_seq = request.wt_sequence
    active_positions: Optional[Set[int]] = set(request.active_site_positions) if request.active_site_positions else None

    if not wt_seq:
        if "SCFV" in target_clean or "2YBR" in target_clean or "SCORPION" in target_clean or "CN2" in target_clean:
            scaffold_path = DATA_DIR / "SCFV_9004G_CN2.json"
            if scaffold_path.is_file():
                meta = json.loads(scaffold_path.read_text(encoding="utf-8"))
                wt_seq = meta.get("sequence", {}).get("value")
        elif "TRPB" in target_clean or "8EH0" in target_clean:
            meta = fetch_protein_metadata("TRPB_TM9D8")
            wt_seq = meta.get("sequence")
        elif "P42212" in target_clean or "GFP" in target_clean or "1EMA" in target_clean:
            meta = fetch_protein_metadata("P42212")
            wt_seq = meta.get("sequence")
        elif "GB1" in target_clean or "2QMT" in target_clean or "SPG1" in target_clean:
            meta = fetch_protein_metadata("SPG1_STRSG")
            wt_seq = meta.get("sequence")
        elif "HER2" in target_clean or "1N8Z" in target_clean or "TRASTUZUMAB" in target_clean:
            meta = fetch_protein_metadata("HER2_TRASTUZUMAB")
            wt_seq = meta.get("sequence")

    if not wt_seq:
        meta = fetch_protein_metadata("P42212")
        wt_seq = meta.get("sequence")

    if not wt_seq:
        raise HTTPException(status_code=400, detail=f"Could not resolve wild-type sequence for {request.target_protein}")

    candidates: List[Variant] = []
    evidence = ledger_db.get_evidence(db_path=DEFAULT_DB_PATH)
    score_map: Dict[str, float] = {}
    for ev in evidence:
        s = ev.get("variant_sequence")
        sc = ev.get("computed_score")
        if s and sc is not None:
            score_map[s] = float(sc)

    if request.variants:
        for notat in request.variants:
            try:
                v = Variant(notation=notat)
                v.validate_against_wt(wt_seq)
                candidates.append(v)
            except Exception as e:
                logger.warning("Skipping invalid variant notation %s: %s", notat, e)
    else:
        for ev in evidence:
            s = ev.get("variant_sequence")
            if s and len(s) == len(wt_seq) and len(candidates) < 30:
                try:
                    candidates.append(Variant.from_sequence(s, wt_seq))
                except Exception:
                    pass

    if not candidates:
        hotspot_pos = [64, 65, 68]
        if "SCFV" in target_clean:
            hotspot_pos = [98, 99, 100]
        elif "TRPB" in target_clean:
            hotspot_pos = [183, 184, 227]
        for p in hotspot_pos:
            if p <= len(wt_seq):
                orig_aa = wt_seq[p - 1]
                for mut_aa in ["A", "V", "L", "S", "T"]:
                    if mut_aa != orig_aa and len(candidates) < 15:
                        candidates.append(Variant(
                            notation=f"{orig_aa}{p}{mut_aa}",
                            mutations=[SingleMutation(wt_aa=orig_aa, position=p, mut_aa=mut_aa)],
                        ))

    predicted_fitness_list: List[float] = []
    stability_priors_list: List[float] = []

    for c in candidates:
        seq = c.apply_to(wt_seq)
        if seq in score_map:
            fit = score_map[seq]
        else:
            fit = compute_biological_prior(c)
        predicted_fitness_list.append(fit)
        stability_priors_list.append(compute_biological_prior(c) * 0.5)

    fitness_arr = np.array(predicted_fitness_list, dtype=np.float64)
    stab_arr = np.array(stability_priors_list, dtype=np.float64)

    result = compute_multi_objective_pareto_front(
        candidates=candidates,
        wt_sequence=wt_seq,
        predicted_fitness=fitness_arr,
        stability_priors=stab_arr,
        weights=request.weights,
        active_site_positions=active_positions,
    )

    chart_data = [
        {
            "label": p.notation,
            "x_fitness": round(p.fitness, 3),
            "y_developability": round(p.developability_score, 3),
            "synthesizability": round(p.synthesizability_score, 3),
            "stability": round(p.stability_score, 3),
            "rank": p.rank,
            "is_frontier": p.rank == 0,
            "weighted_utility": round(p.weighted_utility, 3),
            "liabilities_count": len(p.liabilities.get("clinical_liabilities", [])),
            "clinical_liabilities": p.liabilities.get("clinical_liabilities", []),
        }
        for p in result.all_points
    ]

    return {
        "target_protein": request.target_protein,
        "points": [p.model_dump() for p in result.all_points],
        "frontier_points": [p.model_dump() for p in result.frontier_points],
        "rank_0_count": result.rank_0_count,
        "dominated_count": result.dominated_count,
        "hypervolume": result.hypervolume,
        "weights_used": result.weights_used,
        "chart_data": chart_data,
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("server:app", host="127.0.0.1", port=8000, reload=True)

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
import secrets
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from fastapi import FastAPI, HTTPException, Query, Request, status
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from controller.auth import (
    AuthenticatedUser,
    authenticate_user_account,
    create_part11_electronic_signature,
    register_user_account,
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

# CORS allowlist. `allow_origins=["*"]` together with `allow_credentials=True`
# makes the middleware echo back whichever Origin the request carries, so any
# website a logged-in operator visited could issue credentialed calls against this
# API. Origins come from VERITAS_CORS_ORIGINS (comma-separated); the default is the
# local Vite dev server only.
_DEFAULT_CORS_ORIGINS = "http://localhost:5173,http://127.0.0.1:5173"
CORS_ORIGINS = [
    origin.strip()
    for origin in os.environ.get("VERITAS_CORS_ORIGINS", _DEFAULT_CORS_ORIGINS).split(",")
    if origin.strip()
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["Authorization", "Content-Type"],
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


class RegisterUserRequest(BaseModel):
    # No `role` field: registration is unauthenticated, so the role must not be
    # client-influenced. register_user_account() mints the lowest-privilege role.
    email: str
    password: str
    display_name: str = ""
    organization_id: str = "panoptesbio"


class LoginUserRequest(BaseModel):
    email: str
    password: str


class SignActionRequest(BaseModel):
    id_token: Optional[str] = None
    action: str = "Approval of Modal A10G Compute Spend"
    meaning: str = "Authorized compute budget spend under 21 CFR §11.200"


class LaunchCampaignRequest(BaseModel):
    target_accession: str = "TRPB_TM9D8"
    total_budget_usd: float = Field(default=1.00, gt=0, le=1000.0)
    max_iterations: int = Field(default=4, ge=1, le=100)
    batch_size: int = Field(default=15, ge=1, le=500)
    target_model: str = "mock_esm"
    enable_pocket_hopping: bool = True
    auth_token: Optional[str] = None


class ParetoFrontierRequest(BaseModel):
    target_protein: str = "avGFP"
    wt_sequence: Optional[str] = None
    variants: Optional[List[str]] = None
    weights: Optional[Dict[str, float]] = None
    active_site_positions: Optional[List[int]] = None
    auth_token: Optional[str] = None


class DeNovoRequestPayload(BaseModel):
    target_scaffold: str = "SCFV_9004G_CN2"
    num_designs: int = Field(default=3, ge=1, le=100)
    binder_length: int = Field(default=65, ge=10, le=2000)
    budget_usd: float = Field(default=0.50, gt=0, le=1000.0)
    auth_token: Optional[str] = None


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
@app.get("/api/campaign/pdb/{target_id}")
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


@app.post("/api/auth/register")
def register_account(req: RegisterUserRequest) -> Dict[str, Any]:
    """Creates a new operator account and issues an authenticated session token."""
    try:
        user, token = register_user_account(
            email=req.email,
            password=req.password,
            display_name=req.display_name,
            organization_id=req.organization_id,
        )
        return {
            "authenticated": True,
            "user": user.model_dump(),
            "id_token": token,
        }
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/auth/login")
def login_account(req: LoginUserRequest) -> Dict[str, Any]:
    """Authenticates an existing operator account and issues a session token."""
    try:
        user, token = authenticate_user_account(
            email=req.email,
            password=req.password,
        )
        return {
            "authenticated": True,
            "user": user.model_dump(),
            "id_token": token,
        }
    except Exception as exc:
        raise HTTPException(status_code=401, detail=str(exc)) from exc


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



class PrepareCampaignStreamRequest(BaseModel):
    """Parameters for a campaign the caller is asking to have authorized."""

    target_accession: str = Field(default="TRPB_TM9D8", max_length=64)
    total_budget_usd: float = Field(default=1.00, gt=0, le=1000.0)
    max_iterations: int = Field(default=4, ge=1, le=100)
    batch_size: int = Field(default=15, ge=1, le=500)
    target_model: str = Field(default="mock_esm", max_length=64)
    enable_pocket_hopping: bool = Field(default=True)
    auth_token: Optional[str] = Field(default=None)


# Campaigns authorized by POST and awaiting their SSE replay connection.
# Single-use: the ticket is consumed when the stream opens, so a replayed or
# guessed URL starts nothing.
_PENDING_CAMPAIGN_STREAMS: Dict[str, Dict[str, Any]] = {}
_MAX_PENDING_STREAMS = 32


@app.post("/api/campaign/prepare-stream")
def prepare_campaign_stream(request: PrepareCampaignStreamRequest) -> Dict[str, Any]:
    """
    Authorize a campaign and return a single-use ticket for its SSE stream.

    Spending has to happen behind a POST. The streaming endpoint used to take the
    budget, batch size and model as GET query parameters and start a real campaign
    that moved real money, which any page could trigger with an <img> tag or a
    link - CORS does not stop the request being sent, only the response being read.
    The GET now carries nothing but an opaque ticket minted here.
    """
    if len(_PENDING_CAMPAIGN_STREAMS) >= _MAX_PENDING_STREAMS:
        # Drop the oldest rather than growing without bound.
        oldest = min(_PENDING_CAMPAIGN_STREAMS, key=lambda k: _PENDING_CAMPAIGN_STREAMS[k]["created"])
        _PENDING_CAMPAIGN_STREAMS.pop(oldest, None)

    ticket = secrets.token_urlsafe(32)
    _PENDING_CAMPAIGN_STREAMS[ticket] = {
        "config": request.model_dump(exclude={"auth_token"}),
        "created": time.time(),
    }
    return {"stream_ticket": ticket, "stream_url": f"/api/campaign/stream/{ticket}"}


@app.get("/api/campaign/stream/{ticket}")
def stream_campaign(ticket: str):
    """
    Server-Sent Events (SSE) replay of a campaign authorized by POST.

    Carries no spend parameters: the ticket is the only input, and it is consumed
    on connect so it cannot be reused.
    """
    pending = _PENDING_CAMPAIGN_STREAMS.pop(ticket, None)
    if pending is None:
        raise HTTPException(
            status_code=404,
            detail=(
                "Unknown or already-used stream ticket. POST to "
                "/api/campaign/prepare-stream first."
            ),
        )
    params = pending["config"]

    def event_generator():
        try:
            config = OrchestratorConfig(
                target_accession=params["target_accession"],
                total_budget_usd=params["total_budget_usd"],
                max_iterations=params["max_iterations"],
                batch_size=params["batch_size"],
                target_model=params["target_model"],
                use_mock_llm=True,
                enable_live_llm=False,
                enable_regret_stopping=True,
                enable_pocket_hopping=params["enable_pocket_hopping"],
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
        auth_token = None

        if "application/json" in content_type:
            body = await request.json()
            csv_text = body.get("csv_content", "")
            target_accession = body.get("target_accession", "P42212")
            normalize_to_wt = bool(body.get("normalize_to_wt", True))
            auth_token = body.get("auth_token")
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
            auth_token = form.get("auth_token")
            auth_token = str(auth_token) if auth_token else None
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
            wt_seq = meta.sequence
        elif "P42212" in target_clean or "GFP" in target_clean or "1EMA" in target_clean:
            meta = fetch_protein_metadata("P42212")
            wt_seq = meta.sequence
        elif "GB1" in target_clean or "2QMT" in target_clean or "SPG1" in target_clean:
            meta = fetch_protein_metadata("SPG1_STRSG")
            wt_seq = meta.sequence
        elif "HER2" in target_clean or "1N8Z" in target_clean or "TRASTUZUMAB" in target_clean:
            meta = fetch_protein_metadata("HER2_TRASTUZUMAB")
            wt_seq = meta.sequence

    if not wt_seq:
        meta = fetch_protein_metadata("P42212")
        wt_seq = meta.sequence

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




class RelaxStructureRequest(BaseModel):
    target_protein: str = Field(default="P42212")
    mutations: List[str] = Field(..., description="List of mutations, e.g. ['A206W'] or ['S65T']")
    max_iterations: int = Field(default=200, ge=1, le=10000)
    harmonic_k: float = Field(default=50.0, gt=0, le=10000.0)


@app.post("/api/campaign/relax-structure")
async def relax_structure_endpoint(request: RelaxStructureRequest):
    """
    Track C: Fast CPU L-BFGS sidechain rotamer & local energy relaxation.
    Resolves false-negative steric clashes in candidate mutant protein structures.
    """
    from workers.relaxation_worker import relax_variant_structure

    pdb_filename = "1EMA.pdb"
    target_upper = request.target_protein.upper()
    if "TRPB" in target_upper or "8EH0" in target_upper:
        pdb_filename = "8EH0.pdb"
    elif "SCFV" in target_upper or "2YBR" in target_upper or "CN2" in target_upper:
        pdb_filename = "2YBR.pdb"
    elif "TRASTUZUMAB" in target_upper or "HER2" in target_upper or "1N8Z" in target_upper:
        pdb_filename = "1N8Z.pdb"
    elif "GB1" in target_upper or "2QMT" in target_upper or "SPG1" in target_upper:
        pdb_filename = "2QMT.pdb"

    pdb_path = DATA_DIR / pdb_filename
    if not pdb_path.is_file():
        pdb_path = DATA_DIR / "1EMA.pdb"

    try:
        res = relax_variant_structure(
            pdb_path=pdb_path,
            mutations=request.mutations,
            max_iterations=request.max_iterations,
            harmonic_k=request.harmonic_k,
        )
        return {
            "success": True,
            "target_protein": request.target_protein,
            "pdb_filename": pdb_filename,
            "result": res.model_dump(),
        }
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Relaxation failed: {str(exc)}")


class CloudLabDispatchRequest(BaseModel):
    provider: str = Field(default="emerald_cloud_lab", description="emerald_cloud_lab or strateos")
    target_accession: str = Field(default="P42212")
    variants: List[str] = Field(..., min_length=1)
    fantasy_scores: Optional[List[float]] = Field(default=None)
    assay_type: str = Field(default="fluorescence")
    job_id: Optional[str] = Field(default=None)
    nominal_cost_usd: float = Field(default=0.000500)
    # db_path deliberately removed: the ledger location is server-side state.
    # A client-supplied path let a caller point the audit chain and budget at a
    # fresh database, which combined with on-demand budget creation made the
    # spend cap unenforceable - request a new path, get a new $10 allocation.


@app.post("/api/cloudlab/dispatch")
async def cloudlab_dispatch_endpoint(request: CloudLabDispatchRequest):
    """
    SPEC-001 Seam 5: Dispatches robotic synthesis & assay protocol to Emerald Cloud Lab (SLL)
    or Strateos (Autoprotocol), stages Kriging-believer fantasy observations, and transitions
    the ledger job to IN_FLIGHT_SYNTHESIS.
    """
    import uuid
    from controller.gatekeeper import Variant
    from workers.cloudlab_worker import EmeraldCloudLabClient, StrateosClient
    from workers.rag_worker import fetch_protein_metadata

    target_db = DEFAULT_DB_PATH
    ledger_db.init_db(target_db)
    if not ledger_db.get_budget(1, db_path=target_db):
        # Fail closed rather than minting a budget on demand. Auto-creating a $10
        # allocation here meant the cap could be reset by any caller.
        raise HTTPException(
            status_code=409,
            detail=(
                "No budget is provisioned on this server. An operator must allocate "
                "one before compute or wet-lab spend can be authorized."
            ),
        )

    job_id = request.job_id or f"job-cloudlab-{int(time.time())}-{uuid.uuid4().hex[:6]}"
    existing_job = ledger_db.get_job(job_id, db_path=target_db)
    if existing_job is None:
        authorized = ledger_db.authorize_spend(
            job_id=job_id,
            target_model=f"cloudlab_{request.provider.lower()}",
            variant_count=len(request.variants),
            estimated_cost=request.nominal_cost_usd,
            db_path=target_db,
        )
        if not authorized:
            raise HTTPException(status_code=402, detail="Ledger rejected cloud-lab spend authorization")

    meta = fetch_protein_metadata(request.target_accession)
    wt_seq = meta.sequence

    fantasy_list: List[Tuple[str, float]] = []
    for idx, var_str in enumerate(request.variants):
        f_score = (
            float(request.fantasy_scores[idx])
            if request.fantasy_scores and idx < len(request.fantasy_scores)
            else 1.25
        )
        try:
            seq = wt_seq if var_str == "WT" else Variant(notation=var_str).apply_to(wt_seq)
        except Exception:
            seq = var_str
        fantasy_list.append((seq, f_score))

    cur_job = ledger_db.get_job(job_id, db_path=target_db)
    if cur_job and cur_job["status"] == "APPROVED":
        ledger_db.stage_speculative_evidence(job_id=job_id, fantasy_results=fantasy_list, db_path=target_db)

    if request.provider.lower() in ("strateos", "autoprotocol"):
        dispatch_info = StrateosClient().dispatch_protocol(campaign_id=job_id, variants=request.variants)
    else:
        dispatch_info = EmeraldCloudLabClient().dispatch_protocol(
            campaign_id=job_id,
            variants=request.variants,
            assay_type=request.assay_type,
        )

    ledger_db.mark_job_in_flight(
        job_id=job_id,
        provider=dispatch_info["provider"],
        external_run_id=dispatch_info["external_run_id"],
        simulated=bool(dispatch_info.get("simulated", True)),
        db_path=target_db,
    )

    return {
        "success": True,
        "job_id": job_id,
        "status": "IN_FLIGHT_SYNTHESIS",
        "provider": dispatch_info["provider"],
        "external_run_id": dispatch_info["external_run_id"],
        "protocol_format": dispatch_info["protocol_format"],
        "protocol_preview": dispatch_info["protocol_preview"],
    }


@app.post("/api/cloudlab/webhook")
async def cloudlab_webhook_endpoint(request: Request):
    """
    SPEC-001 Seam 5: HMAC-SHA256 authenticated webhook callback from remote cloud laboratories.
    Cleanses replicate outliers via Grubbs'/Huber filtering and atomically reconciles
    provisional Kriging-believer rows into verified Evidence.
    """
    from controller.gatekeeper import Variant
    from scripts.ingest_plate_data import parse_plate_csv
    from workers.cloudlab_worker import verify_cloudlab_webhook_signature
    from workers.rag_worker import fetch_protein_metadata

    raw_bytes = await request.body()
    sig_header = request.headers.get("X-CloudLab-Signature", "")
    if not verify_cloudlab_webhook_signature(raw_bytes, sig_header):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid or missing X-CloudLab-Signature HMAC header",
        )

    try:
        payload = json.loads(raw_bytes.decode("utf-8"))
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Malformed webhook JSON: {exc}")

    job_id = payload.get("job_id")
    csv_data = payload.get("csv_data")
    target_accession = payload.get("target_accession", "P42212")
    # The ledger location is server-side state and is never taken from the payload,
    # signed or not: a caller able to choose the database could write "physical
    # evidence" into a chain of their own and have it accepted as reconciled.
    target_db = DEFAULT_DB_PATH

    if not job_id or not csv_data:
        raise HTTPException(status_code=400, detail="Webhook payload requires job_id and csv_data")

    parsed = parse_plate_csv(csv_data, normalize_to_wt=True, reject_replicate_outliers=True)
    meta = fetch_protein_metadata(target_accession)
    wt_seq = meta.sequence

    observed_tuples: List[Tuple[str, float]] = []
    for v in parsed["valid_variants"]:
        mut_str = v["mutations"]
        score = float(v["normalized_score"])
        try:
            if mut_str in ("WT", "WT_Control"):
                seq = wt_seq
            else:
                seq = Variant(notation=mut_str).apply_to(wt_seq)
        except Exception:
            seq = mut_str
        observed_tuples.append((seq, score))

    reconciled_count = ledger_db.commit_physical_evidence_and_reconcile(
        job_id=job_id,
        observed_results=observed_tuples,
        db_path=target_db,
    )
    job_now = ledger_db.get_job(job_id, db_path=target_db)

    return {
        "success": True,
        "job_id": job_id,
        "job_status": job_now["status"] if job_now else "COMPLETED",
        "reconciled_count": reconciled_count,
        "outliers_rejected": parsed.get("outliers_rejected", []),
        "valid_variants": parsed["valid_variants"],
    }


class CurveSeriesItem(BaseModel):
    mutations: str = Field(..., description="Variant notation e.g. 'WT' or 'S65T:Y66H'")
    x_values: List[float] = Field(..., description="Concentrations (nM) for 4PL or timepoints (sec) for kinetics")
    y_values: List[float] = Field(..., description="Measured responses or absorbances")


class CurveIngestionRequest(BaseModel):
    curve_type: str = Field(default="dose_response_4pl", description="dose_response_4pl or michaelis_menten")
    target_accession: str = Field(default="P42212")
    series: List[CurveSeriesItem] = Field(..., min_length=1)
    nominal_cost_usd: float = Field(default=0.000200)
    # db_path deliberately removed: the ledger location is server-side state.
    # A client-supplied path let a caller point the audit chain and budget at a
    # fresh database, which combined with on-demand budget creation made the
    # spend cap unenforceable - request a new path, get a new $10 allocation.


@app.post("/api/campaign/ingest-curves")
async def ingest_curves_endpoint(request: CurveIngestionRequest):
    """
    SPEC-002 Seam 4: Fits 4PL dose-response (IC50/EC50/Kd) or Michaelis-Menten initial velocity (v0)
    curves across candidate dilution/time series, normalizes potency against WT, and records
    variance-weighted evidence into the SQLite ledger.
    """
    import uuid
    from controller.gatekeeper import Variant
    from scripts.ingest_plate_data import fit_dose_response_4pl, fit_michaelis_menten_kinetics
    from workers.rag_worker import fetch_protein_metadata

    target_db = DEFAULT_DB_PATH
    ledger_db.init_db(target_db)
    if not ledger_db.get_budget(1, db_path=target_db):
        # Fail closed rather than minting a budget on demand. Auto-creating a $10
        # allocation here meant the cap could be reset by any caller.
        raise HTTPException(
            status_code=409,
            detail=(
                "No budget is provisioned on this server. An operator must allocate "
                "one before compute or wet-lab spend can be authorized."
            ),
        )

    meta = fetch_protein_metadata(request.target_accession)
    wt_seq = meta.sequence

    fitted_items: List[Dict[str, Any]] = []
    wt_ref_metric: Optional[float] = None

    for item in request.series:
        if request.curve_type.lower() in ("michaelis_menten", "kinetics"):
            fit_res = fit_michaelis_menten_kinetics(item.x_values, item.y_values)
            raw_metric = fit_res["v0_um_per_sec"]
            variance = fit_res["v0_variance"]
        else:
            fit_res = fit_dose_response_4pl(item.x_values, item.y_values)
            raw_metric = fit_res["ic50"]
            variance = fit_res["potency_variance"]

        if item.mutations.upper() in ("WT", "WT_CONTROL") and wt_ref_metric is None:
            wt_ref_metric = raw_metric

        fitted_items.append({
            "mutations": item.mutations,
            "raw_metric": raw_metric,
            "variance": variance,
            "fit_details": fit_res,
        })

    evidence_tuples: List[Tuple[str, float]] = []
    for f_item in fitted_items:
        mut_str = f_item["mutations"]
        if request.curve_type.lower() in ("michaelis_menten", "kinetics"):
            norm_score = (f_item["raw_metric"] / wt_ref_metric) if (wt_ref_metric and wt_ref_metric > 0) else f_item["raw_metric"]
        else:
            # Lower IC50/Kd = higher potency: fold improvement = IC50_WT / IC50_variant
            norm_score = (wt_ref_metric / max(f_item["raw_metric"], 1e-9)) if (wt_ref_metric and wt_ref_metric > 0) else (100.0 / max(f_item["raw_metric"], 1e-6))

        f_item["normalized_score"] = round(float(norm_score), 4)
        try:
            seq = wt_seq if mut_str.upper() in ("WT", "WT_CONTROL") else Variant(notation=mut_str).apply_to(wt_seq)
        except Exception:
            seq = mut_str
        evidence_tuples.append((seq, float(norm_score)))

    job_id = f"job-curves-{int(time.time())}-{uuid.uuid4().hex[:6]}"
    if not ledger_db.authorize_spend(
        job_id=job_id,
        target_model=f"wetlab_{request.curve_type.lower()}",
        variant_count=len(evidence_tuples),
        estimated_cost=request.nominal_cost_usd,
        db_path=target_db,
    ):
        raise HTTPException(status_code=402, detail="Ledger rejected curve ingestion spend")

    ledger_db.record_evidence(job_id=job_id, results=evidence_tuples, db_path=target_db)
    ledger_db.mark_job_completed(job_id=job_id, db_path=target_db)

    return {
        "success": True,
        "job_id": job_id,
        "curve_type": request.curve_type,
        "variants_ingested": len(evidence_tuples),
        "fitted_curves": fitted_items,
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("server:app", host="127.0.0.1", port=8008, reload=True)

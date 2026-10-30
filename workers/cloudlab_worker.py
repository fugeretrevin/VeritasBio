"""
workers/cloudlab_worker.py — Autonomous Cloud-Lab Adapter (SPEC-001 Seam 5).

Purpose:
    Programmatic client adapters for Emerald Cloud Lab (Symbolic Lab Language, SLL)
    and Strateos (Autoprotocol JSON) remote robotic wet-lab execution, plus
    HMAC-SHA256 webhook authentication for automated microplate reader callbacks.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import logging
import os
import time
import uuid
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# No committed default. The previous literal was public in the repository, so any
# reader could forge an X-CloudLab-Signature header and post fabricated assay
# results that the ledger would accept as reconciled physical evidence. When unset,
# signing and verification both refuse rather than falling back to a known key.
CLOUDLAB_WEBHOOK_SECRET: bytes = os.environ.get("VERITAS_CLOUDLAB_WEBHOOK_SECRET", "").encode("utf-8")


class WebhookSecretNotConfigured(RuntimeError):
    """Raised when webhook HMAC operations are attempted with no secret configured."""


def compute_cloudlab_webhook_signature(
    payload_bytes: bytes,
    secret: Optional[bytes] = None,
) -> str:
    """Computes HMAC-SHA256 hex digest for cloud-lab webhook payloads."""
    key = secret if secret is not None else CLOUDLAB_WEBHOOK_SECRET
    if not key:
        raise WebhookSecretNotConfigured(
            "VERITAS_CLOUDLAB_WEBHOOK_SECRET is not set; refusing to compute a webhook "
            "signature with an empty key."
        )
    return hmac.new(key, payload_bytes, hashlib.sha256).hexdigest()


def verify_cloudlab_webhook_signature(
    payload_bytes: bytes,
    signature_header: str,
    secret: Optional[bytes] = None,
) -> bool:
    """Constant-time verification of X-CloudLab-Signature HMAC header."""
    if not signature_header or not isinstance(signature_header, str):
        return False
    clean_sig = signature_header.strip()
    if clean_sig.startswith("sha256="):
        clean_sig = clean_sig[len("sha256="):]
    try:
        expected = compute_cloudlab_webhook_signature(payload_bytes, secret=secret)
    except WebhookSecretNotConfigured:
        # Fail closed: with no secret configured every signature is unverifiable,
        # so no webhook can be trusted.
        logger.error(
            "Rejecting cloud-lab webhook: VERITAS_CLOUDLAB_WEBHOOK_SECRET is not configured."
        )
        return False
    return hmac.compare_digest(clean_sig, expected)


class EmeraldCloudLabClient:
    """
    Generates Emerald Cloud Lab Symbolic Lab Language (SLL) protocols locally.

    Protocol generation is real; dispatch is simulated. This module contains no
    HTTP client, so nothing is ever submitted to Emerald and `external_run_id` is
    a locally minted identifier, not an ECL run. `simulated` is therefore always
    True, regardless of whether ECL_API_TOKEN is set - a token previously flipped
    `mock_mode` to False while still submitting nothing, which made the return
    value claim a live run that had not happened.
    """

    def __init__(self, api_token: Optional[str] = None) -> None:
        self.api_token = api_token or os.environ.get("ECL_API_TOKEN", "")
        self.has_credentials = bool(self.api_token)
        # No transport exists in this module; dispatch cannot be live.
        self.simulated = True
        self.mock_mode = True

    def build_sll_protocol(
        self,
        campaign_id: str,
        variants: List[str],
        assay_type: str = "fluorescence",
        well_volume_ul: float = 100.0,
    ) -> str:
        """Generates an Emerald Cloud Lab Symbolic Lab Language (SLL) notebook expression."""
        readout_fn = (
            "ExperimentFluorescenceIntensity"
            if assay_type.lower() != "absorbance"
            else "ExperimentAbsorbanceIntensity"
        )
        variant_str = ", ".join(f'"{v}"' for v in variants)
        return (
            f"(* VeritasBio Autonomous SLL Protocol — Campaign {campaign_id} *)\n"
            f"Module[{{assayPlate, transferPrimitives, readoutProtocol}},\n"
            f"  assayPlate = Model[Container, Plate, \"96-well Corning Black Flat Bottom\"];\n"
            f"  transferPrimitives = {{\n"
            f"    (* Perimeter Hydration Barrier: Rows A/H, Cols 1/12 *)\n"
            f"    Transfer[Source -> Model[Sample, \"Milli-Q Water\"], Destination -> assayPlate, "
            f"DestinationWell -> {{\"A1\", \"A12\", \"H1\", \"H12\"}}, Amount -> {well_volume_ul} Microliter],\n"
            f"    (* Wild-Type Reference Baseline in Well B02 *)\n"
            f"    Transfer[Source -> Model[Sample, \"WT_Parent_Control\"], Destination -> assayPlate, "
            f"DestinationWell -> \"B2\", Amount -> {well_volume_ul} Microliter],\n"
            f"    (* Candidate Variants: {{{variant_str}}} *)\n"
            f"    ExperimentRoboticLiquidHandling[Variants -> {{{variant_str}}}, "
            f"AssayVolume -> {well_volume_ul} Microliter, Replicates -> 3]\n"
            f"  }};\n"
            f"  ExperimentIncubate[assayPlate, Temperature -> 37 Celsius, Time -> 4 Hour];\n"
            f"  readoutProtocol = {readout_fn}[\n"
            f"    SamplesIn -> assayPlate,\n"
            f"    ExcitationWavelength -> 488 Nanometer,\n"
            f"    EmissionWavelength -> 510 Nanometer\n"
            f"  ]\n"
            f"]"
        )

    def dispatch_protocol(
        self,
        campaign_id: str,
        variants: List[str],
        assay_type: str = "fluorescence",
    ) -> Dict[str, Any]:
        sll_script = self.build_sll_protocol(
            campaign_id=campaign_id,
            variants=variants,
            assay_type=assay_type,
        )
        run_id = f"ecl-sim-{int(time.time())}-{uuid.uuid4().hex[:6]}"
        return {
            "provider": "emerald_cloud_lab",
            "external_run_id": run_id,
            "external_run_id_is_simulated": True,
            "simulated": self.simulated,
            "mock_mode": self.mock_mode,
            "has_credentials": self.has_credentials,
            "protocol_format": "SLL",
            "protocol_preview": sll_script,
        }


class StrateosClient:
    """
    Generates Strateos Autoprotocol JSON locally.

    Protocol generation is real; dispatch is simulated. This module contains no
    HTTP client, so nothing is submitted to Strateos and `external_run_id` is a
    locally minted identifier rather than a Strateos run. `simulated` is always
    True regardless of STRATEOS_API_TOKEN.
    """

    def __init__(self, api_token: Optional[str] = None, organization_id: Optional[str] = None) -> None:
        self.api_token = api_token or os.environ.get("STRATEOS_API_TOKEN", "")
        self.organization_id = organization_id or os.environ.get("STRATEOS_ORG_ID", "veritasbio-cro")
        self.has_credentials = bool(self.api_token)
        # No transport exists in this module; dispatch cannot be live.
        self.simulated = True
        self.mock_mode = True

    def build_autoprotocol(
        self,
        campaign_id: str,
        variants: List[str],
        well_volume_ul: float = 100.0,
    ) -> Dict[str, Any]:
        """Generates a Strateos-compliant Autoprotocol JSON payload."""
        inner_rows = ["B", "C", "D", "E", "F", "G"]
        inner_cols = list(range(2, 12))
        available_wells = [f"{r}{c}" for r in inner_rows for c in inner_cols if not (r == "B" and c == 2)]

        transfers: List[Dict[str, Any]] = [
            {
                "from": "wt_stock/0",
                "to": "assay_plate/B2",
                "volume": f"{well_volume_ul}:microliter",
                "sample_name": "WT",
                "sample_type": "control_wt",
            }
        ]
        well_cursor = 0
        for var_str in variants:
            for rep in range(3):
                if well_cursor >= len(available_wells):
                    break
                target_well = available_wells[well_cursor]
                well_cursor += 1
                transfers.append({
                    "from": f"cDNA_expression_stock/{var_str}",
                    "to": f"assay_plate/{target_well}",
                    "volume": f"{well_volume_ul}:microliter",
                    "sample_name": var_str,
                    "replicate": rep + 1,
                })

        return {
            "campaign_id": campaign_id,
            "refs": {
                "assay_plate": {
                    "new": "96-flat-black",
                    "store": {"where": "cold_4"},
                }
            },
            "instructions": [
                {
                    "op": "pipette",
                    "groups": [{"transfer": transfers}],
                },
                {
                    "op": "incubate",
                    "object": "assay_plate",
                    "where": "warm_37",
                    "duration": "4:hour",
                    "shaking": True,
                },
                {
                    "op": "fluorescence",
                    "object": "assay_plate",
                    "excitation": "488:nanometer",
                    "emission": "510:nanometer",
                    "dataref": f"readout_{campaign_id}",
                },
            ],
        }

    def dispatch_protocol(
        self,
        campaign_id: str,
        variants: List[str],
    ) -> Dict[str, Any]:
        autoprotocol = self.build_autoprotocol(campaign_id=campaign_id, variants=variants)
        run_id = f"stx-sim-{int(time.time())}-{uuid.uuid4().hex[:6]}"
        return {
            "provider": "strateos",
            "external_run_id": run_id,
            "external_run_id_is_simulated": True,
            "simulated": self.simulated,
            "mock_mode": self.mock_mode,
            "has_credentials": self.has_credentials,
            "protocol_format": "Autoprotocol",
            "protocol_preview": json.dumps(autoprotocol, indent=2),
            "autoprotocol": autoprotocol,
        }


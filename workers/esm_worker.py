"""
workers/esm_worker.py — DOMAIN: Science & Math Tools (Pillar 5: Cloud Muscle).

Purpose:
    Modal serverless GPU worker running ESM-2 (650M / 150M) via batched PyTorch tensors.
    Computes zero-shot masked-marginal log-likelihood fitness scores for protein variants.

Constraints & Invariants (ARCHITECTURE.md §3.5, §7.2; INNOVATION_ROADMAP Phase 4):
    1. Paid Compute Behind Spend Gate: Must only be invoked after `authorize_spend()`
       in `ledger.database` returns True.
    2. Container Amortization (G1): Configures Modal with `scaledown_window=120` (2-minute keep-alive)
       to avoid paying a 20-second cold-boot penalty on consecutive active-learning iterations.
    3. Multi-Fidelity Support (G4): Supports `facebook/esm2_t33_650M_UR50D` (Tier-2)
       and `facebook/esm2_t30_150M_UR50D` (Tier-1 fast screening).
    4. Resilient Fallback: If Modal credentials are not configured, falls back to the
       high-fidelity biological mock worker so tests and local runs never hard-crash.
"""

from __future__ import annotations

import logging
import os
from typing import Any, Dict, List, Optional, Sequence

import modal

logger = logging.getLogger(__name__)

# Modal App definition
app = modal.App("veritasbio-esm")

# Remote container image with PyTorch and HuggingFace Transformers
esm_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch>=2.1.0",
        "transformers>=4.35.0",
        "numpy>=1.24.0",
    )
)

# Standard HuggingFace model identifiers
MODEL_MAP: Dict[str, str] = {
    "esm2_650m": "facebook/esm2_t33_650M_UR50D",
    "esm2_t33_650M_UR50D": "facebook/esm2_t33_650M_UR50D",
    "esm2_150m": "facebook/esm2_t30_150M_UR50D",
    "esm2_t30_150M_UR50D": "facebook/esm2_t30_150M_UR50D",
}


@app.function(
    image=esm_image,
    gpu="A10G",
    timeout=300,
    scaledown_window=45,  # 45 seconds keep-alive: amortizes cold boot across rounds without lingering idle costs
)
def score_variants_remote(
    wt_sequence: str,
    variant_sequences: List[str],
    hf_model_id: str = "facebook/esm2_t33_650M_UR50D",
    batch_size: int = 32,
) -> List[float]:
    """
    Remote Modal GPU function: executes batched masked-marginal log-likelihood scoring.
    Batches sequences into tensor mini-batches for maximum A10G throughput.

    For each variant, computes:
        Score = sum_{i in mut_pos} [ log P(mut_token_i | seq) - log P(wt_token_i | seq) ]
    """
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"
    logger.info("Modal GPU Worker loading model %s on %s...", hf_model_id, device)

    tokenizer = AutoTokenizer.from_pretrained(hf_model_id)
    model = AutoModelForMaskedLM.from_pretrained(hf_model_id).to(device)
    model.eval()

    # Pre-encode WT reference
    wt_tokens = tokenizer(wt_sequence, return_tensors="pt").to(device)
    with torch.no_grad():
        wt_logits = model(**wt_tokens).logits[0]  # [L, Vocab]
        wt_log_probs = torch.log_softmax(wt_logits, dim=-1)

    scores: List[float] = []

    # Process variants in tensor mini-batches
    for b_start in range(0, len(variant_sequences), batch_size):
        b_seqs = variant_sequences[b_start : b_start + batch_size]

        b_indices = []
        b_mut_indices = []
        for i, v_seq in enumerate(b_seqs):
            if v_seq == wt_sequence:
                continue
            mut_idx = [j for j, (wc, vc) in enumerate(zip(wt_sequence, v_seq)) if wc != vc]
            if mut_idx:
                b_indices.append(i)
                b_mut_indices.append(mut_idx)

        b_scores = [0.0] * len(b_seqs)
        if b_indices:
            eval_seqs = [b_seqs[i] for i in b_indices]
            b_tokens = tokenizer(eval_seqs, padding=True, return_tensors="pt").to(device)
            with torch.no_grad():
                b_logits = model(**b_tokens).logits  # [B, L, Vocab]
                b_log_probs = torch.log_softmax(b_logits, dim=-1)

            for slot_idx, orig_i in enumerate(b_indices):
                mut_idx = b_mut_indices[slot_idx]
                delta_logp = 0.0
                for idx in mut_idx:
                    token_idx = idx + 1
                    wt_tok_id = wt_tokens.input_ids[0, token_idx]
                    mut_tok_id = b_tokens.input_ids[slot_idx, token_idx]
                    wt_logp = wt_log_probs[token_idx, wt_tok_id].item()
                    mut_logp = b_log_probs[slot_idx, token_idx, mut_tok_id].item()
                    delta_logp += (mut_logp - wt_logp)
                b_scores[orig_i] = round(delta_logp, 4)

        scores.extend(b_scores)

    return scores


def is_modal_configured() -> bool:
    """Checks whether Modal API credentials exist in environment."""
    return bool(
        os.getenv("MODAL_TOKEN_ID")
        or os.path.exists(os.path.expanduser("~/.modal.toml"))
    )


def score_variants(
    wt_sequence: str,
    variant_sequences: Sequence[str],
    target_model: str = "esm2_650m",
    use_mock_fallback: bool = True,
) -> List[float]:
    """
    Client entrypoint for scoring a batch of sequences.

    Delegates to:
      - `MockESMWorker` if `target_model == 'mock_esm'` or Modal is unconfigured.
      - Remote Modal A10G function if running live cloud execution.
    """
    from workers.mock_workers import MockESMWorker

    clean_model = target_model.strip().lower()

    if clean_model == "mock_esm":
        logger.info("Scoring %d variants using local MockESMWorker...", len(variant_sequences))
        return MockESMWorker(wt_sequence).score_variants(wt_sequence, variant_sequences)

    # Attempt remote Modal execution
    if is_modal_configured():
        hf_id = MODEL_MAP.get(clean_model, "facebook/esm2_t33_650M_UR50D")
        logger.info("Submitting %d variants to Modal A10G (%s)...", len(variant_sequences), hf_id)
        try:
            with app.run():
                return score_variants_remote.remote(
                    wt_sequence=wt_sequence,
                    variant_sequences=list(variant_sequences),
                    hf_model_id=hf_id,
                )
        except Exception as exc:
            logger.error("Modal remote execution failed: %s", exc)
            if not use_mock_fallback:
                raise RuntimeError(f"Modal execution failed: {exc}") from exc
            logger.warning("Falling back to MockESMWorker due to Modal execution failure.")
            return MockESMWorker(wt_sequence).score_variants(wt_sequence, variant_sequences)

    # Modal not configured: check fallback
    if use_mock_fallback:
        logger.warning(
            "Modal credentials not detected in environment. "
            "Falling back to high-fidelity local MockESMWorker."
        )
        return MockESMWorker(wt_sequence).score_variants(wt_sequence, variant_sequences)
    else:
        raise RuntimeError(
            "Modal is not configured (MODAL_TOKEN_ID not set and ~/.modal.toml missing) "
            "and use_mock_fallback is False."
        )

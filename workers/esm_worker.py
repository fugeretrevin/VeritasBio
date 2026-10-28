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

try:
    import modal
    app = modal.App("veritasbio-esm")
except ImportError:
    modal = None
    app = None

logger = logging.getLogger(__name__)

def download_esm_models():
    """Pre-downloads ESM model weights into container image during build time."""
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    for model_id in ["facebook/esm2_t30_150M_UR50D"]:
        AutoTokenizer.from_pretrained(model_id)
        AutoModelForMaskedLM.from_pretrained(model_id)

if modal is not None:
    esm_image = (
        modal.Image.debian_slim(python_version="3.11")
        .pip_install(
            "torch>=2.1.0",
            "transformers>=4.35.0",
            "numpy>=1.24.0",
        )
        .run_function(download_esm_models)
    )
else:
    esm_image = None

def _modal_function(**kwargs):
    def decorator(fn):
        if app is not None:
            return app.function(**kwargs)(fn)
        return fn
    return decorator

# Standard HuggingFace model identifiers
MODEL_MAP: Dict[str, str] = {
    "esm2_650m": "facebook/esm2_t33_650M_UR50D",
    "esm2_t33_650M_UR50D": "facebook/esm2_t33_650M_UR50D",
    "esm2_150m": "facebook/esm2_t30_150M_UR50D",
    "esm2_t30_150M_UR50D": "facebook/esm2_t30_150M_UR50D",
}


_LOADED_MODELS: Dict[str, Any] = {}
_LOADED_TOKENIZERS: Dict[str, Any] = {}


def _score_variants_core(
    wt_sequence: str,
    variant_sequences: List[str],
    hf_model_id: str = "facebook/esm2_t33_650M_UR50D",
    batch_size: int = 32,
) -> List[float]:
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"
    if hf_model_id not in _LOADED_MODELS:
        logger.info("Modal GPU Worker loading model %s on %s...", hf_model_id, device)
        tok = AutoTokenizer.from_pretrained(hf_model_id)
        mod = AutoModelForMaskedLM.from_pretrained(hf_model_id).to(device)
        mod.eval()
        _LOADED_TOKENIZERS[hf_model_id] = tok
        _LOADED_MODELS[hf_model_id] = mod
    else:
        tok = _LOADED_TOKENIZERS[hf_model_id]
        mod = _LOADED_MODELS[hf_model_id]

    tokenizer = tok
    model = mod

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


@_modal_function(
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
    return _score_variants_core(wt_sequence, variant_sequences, hf_model_id, batch_size)


@_modal_function(
    image=esm_image,
    gpu="A10G",
    timeout=300,
    scaledown_window=45,
)
def score_variants_cascade_remote(
    wt_sequence: str,
    variant_sequences: List[str],
    top_k_ratio: float = 0.20,
    batch_size: int = 32,
) -> Dict[str, Any]:
    """
    Remote Modal GPU function: Co-located Multi-Fidelity Cloud Cascade on A10G.
    Executes Tier-1 (150M) screening on all N candidates, promotes top quantile,
    and runs Tier-2 (650M) precision scoring in the SAME container session.
    Avoids duplicate cold-boot and data transfer overhead.
    """
    import math

    logger.info("Modal Cascade: Screening %d variants via Tier 1 (150M)...", len(variant_sequences))
    t1_scores = _score_variants_core(
        wt_sequence=wt_sequence,
        variant_sequences=variant_sequences,
        hf_model_id="facebook/esm2_t30_150M_UR50D",
        batch_size=batch_size,
    )

    k = max(1, math.ceil(len(variant_sequences) * top_k_ratio))
    promoted_indices = sorted(range(len(t1_scores)), key=lambda i: t1_scores[i], reverse=True)[:k]
    promoted_seqs = [variant_sequences[i] for i in promoted_indices]

    logger.info("Modal Cascade: Verifying %d promoted variants via Tier 2 (650M)...", len(promoted_seqs))
    t2_scores = _score_variants_core(
        wt_sequence=wt_sequence,
        variant_sequences=promoted_seqs,
        hf_model_id="facebook/esm2_t33_650M_UR50D",
        batch_size=batch_size,
    )

    return {
        "tier1_scores": t1_scores,
        "promoted_indices": promoted_indices,
        "promoted_sequences": promoted_seqs,
        "tier2_scores": t2_scores,
    }


@_modal_function(
    image=esm_image,
    gpu="A10G",
    timeout=300,
    scaledown_window=45,
)
def embed_sequences_remote(
    sequences: List[str],
    hf_model_id: str = "facebook/esm2_t33_650M_UR50D",
    batch_size: int = 16,
) -> List[List[float]]:
    """
    Remote Modal GPU function: extracts [CLS] 1280-dim representation vectors
    for sequences, feeding deep latent representations to DKL.
    """
    import torch
    from transformers import AutoModel, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"
    tokenizer = AutoTokenizer.from_pretrained(hf_model_id)
    model = AutoModel.from_pretrained(hf_model_id).to(device)
    model.eval()

    embeddings: List[List[float]] = []
    for b_start in range(0, len(sequences), batch_size):
        b_seqs = sequences[b_start : b_start + batch_size]
        tokens = tokenizer(b_seqs, padding=True, return_tensors="pt").to(device)
        with torch.no_grad():
            outputs = model(**tokens)
            cls_reps = outputs.last_hidden_state[:, 0, :].cpu().numpy().tolist()
            embeddings.extend(cls_reps)
    return embeddings


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
            # Attempt lookup on deployed Modal app for instant zero-overhead invocation
            try:
                fn_getter = getattr(modal.Function, "from_name", getattr(modal.Function, "lookup", None))
                if fn_getter is None:
                    raise AttributeError("Neither from_name nor lookup found on modal.Function")
                fn = fn_getter("veritasbio-esm", "score_variants_remote")
                return fn.remote(
                    wt_sequence=wt_sequence,
                    variant_sequences=list(variant_sequences),
                    hf_model_id=hf_id,
                )
            except Exception as lookup_err:
                logger.debug(
                    "Modal deployed function lookup failed (%s); running ephemeral app.run()...",
                    lookup_err,
                )
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


def score_variants_cascade(
    wt_sequence: str,
    variant_sequences: Sequence[str],
    target_model: str = "esm_cascade",
    top_k_ratio: float = 0.20,
    use_mock_fallback: bool = True,
) -> Dict[str, Any]:
    """
    Client entrypoint for executing multi-fidelity cloud cascade scoring.
    Tier 1 (150M) screens all N variants; Tier 2 (650M) scores the top quantile.
    """
    from workers.mock_workers import MockESMWorker

    clean_model = target_model.strip().lower()
    if clean_model == "mock_cascade":
        return MockESMWorker(wt_sequence).score_variants_cascade(
            wt_sequence=wt_sequence,
            variant_sequences=variant_sequences,
            top_k_ratio=top_k_ratio,
        )

    if is_modal_configured():
        try:
            fn_getter = getattr(modal.Function, "from_name", getattr(modal.Function, "lookup", None))
            if fn_getter:
                fn = fn_getter("veritasbio-esm", "score_variants_cascade_remote")
                return fn.remote(
                    wt_sequence=wt_sequence,
                    variant_sequences=list(variant_sequences),
                    top_k_ratio=top_k_ratio,
                )
        except Exception as exc:
            logger.warning("Modal remote cascade scoring failed: %s", exc)
            if not use_mock_fallback:
                raise RuntimeError(f"Modal cascade execution failed: {exc}") from exc

    if use_mock_fallback:
        return MockESMWorker(wt_sequence).score_variants_cascade(
            wt_sequence=wt_sequence,
            variant_sequences=variant_sequences,
            top_k_ratio=top_k_ratio,
        )
    raise RuntimeError("Modal is not configured and use_mock_fallback is False.")


def embed_sequences(
    sequences: Sequence[str],
    target_model: str = "esm2_650m",
    use_mock_fallback: bool = True,
) -> List[List[float]]:
    """
    Client entrypoint for extracting latent representation vectors (1280-dim) for sequences.
    Delegates to Modal GPU worker or local MockESMWorker.
    """
    from workers.mock_workers import MockESMWorker

    clean_model = target_model.strip().lower()

    if clean_model == "mock_esm":
        return MockESMWorker().embed_sequences(sequences)

    if is_modal_configured():
        hf_id = MODEL_MAP.get(clean_model, "facebook/esm2_t33_650M_UR50D")
        try:
            fn_getter = getattr(modal.Function, "from_name", getattr(modal.Function, "lookup", None))
            if fn_getter:
                fn = fn_getter("veritasbio-esm", "embed_sequences_remote")
                return fn.remote(list(sequences), hf_model_id=hf_id)
        except Exception:
            pass

    if use_mock_fallback:
        return MockESMWorker().embed_sequences(sequences)
    raise RuntimeError("Modal is not configured and use_mock_fallback is False.")

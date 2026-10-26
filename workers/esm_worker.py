"""
workers/esm_worker.py — DOMAIN: Science & Math Tools (Pillar 5: Cloud Muscle).

Purpose:
    Modal @function running on serverless A10G GPUs. Boots, loads the ESM
    foundation model, scores variants via batched PyTorch tensors, spins down.

Constraints (ARCHITECTURE.md §7.2):
    - Paid compute. Must only be invoked after authorize_spend() returns True.

Status: scaffold only — not implemented yet.
"""

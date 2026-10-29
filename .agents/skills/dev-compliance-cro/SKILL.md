---
name: dev-compliance-cro
description: >-
  Use this skill when developing wet-lab translation tools, Twist Bioscience cDNA synthesis CSVs,
  Opentrons OT-2 robotics pipetting scripts, Benchling LIMS sync, or FDA 21 CFR Part 11 cryptographic
  audit trails in workers/wetlab_worker.py, workers/lims_worker.py, or ledger/database.py.
---

# Wet-Lab CRO & Compliance Skill

This skill governs the translation of in silico discoveries into actionable physical biology and regulatory compliance.

## Core Invariants
1. **Physical Synthesizability First:**
   - Always run high-expression codon optimization (removing forbidden restriction sites: EcoRI, BamHI, HindIII, XhoI, NotI).
   - Enforce Gibson assembly flanking homology overhangs for expression vectors.
2. **Evaporation Barrier Isolation:**
   - 96-well and 384-well robotic plate maps must support anti-edge barrier isolation (36 perimeter wells reserved for sterile PBS/water) to eliminate meniscus evaporative bias.
3. **FDA 21 CFR Part 11 / GAMP 5 Immutability:**
   - Every state change in the SQLite ledger is backed by an append-only Merkle block with SHA-256 hash linking and HMAC digital signatures.
   - Verification must confirm zero breaks in the block chain (`verify_audit_integrity`).
4. **LIMS Interoperability:**
   - Benchling integration must format standard cDNA registration, 96-well virtual plate positions, and exportable CSV manifests.

## Key APIs
- `codon_optimize_sequence(protein_seq, organism="ecoli")`: High-expression cDNA translation.
- `generate_twist_order_csv(top_variants, price_per_bp=0.09)`: Twist/IDT order portal CSV.
- `generate_plate_layout(variants, plate_format="96_well", anti_edge_barrier=True)`: Automated plate mapping.
- `generate_opentrons_protocol(variant_names, target_protein)`: Valid OT-2 Python protocol script.
- `generate_cryptographic_audit_certificate(...)`: FDA 21 CFR Part 11 compliance certificate.
- `ledger_db.verify_audit_integrity(db_path)`: Complete cryptographic validation of all blocks.

## Verification Checklist
```powershell
python -m unittest tests/test_enterprise_compliance.py
python -m unittest tests/test_wetlab_worker.py
```
Ensure:
- Twist CSV has valid column headers and valid cDNA sequences.
- Opentrons script syntax passes Python compile checks.
- Audit ledger verification returns `(True, count, message)`.

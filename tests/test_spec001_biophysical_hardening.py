"""
tests/test_spec001_biophysical_hardening.py — TDD Verification Suite for SPEC-001.

Covers the public behavioral seams of SPEC-001:
  1. Seam 1: Anisotropic / Gaussian Network Model allosteric hinge detection (`compute_anm_hinge_residues`).
  2. Seam 2: Two-stage elastic mainchain relaxation (`relax_variant_structure` with `allow_backbone_accommodation=True`).
  3. Seam 3: Robust replicate outlier rejection (Grubbs' test & Huber M-estimation) in `parse_plate_csv`.
  4. Seam 4: Asynchronous Kriging-believer ledger staging & physical reconciliation (`stage_speculative_evidence`, `commit_physical_evidence_and_reconcile`).
  5. Seam 5: Autonomous Cloud-Lab SLL/Strateos dispatch & HMAC webhook ingestion (`workers/cloudlab_worker.py`, `server.py`).
"""

from __future__ import annotations

import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PDB_1EMA = PROJECT_ROOT / "workers" / "data" / "1EMA.pdb"
PDB_8EH0 = PROJECT_ROOT / "workers" / "data" / "8EH0.pdb"


class TestSlice1AllostericHinges(unittest.TestCase):
    """Vertical Slice 1: Normal-mode allosteric hinge detection on real crystal structures."""

    def test_compute_anm_hinge_residues_returns_valid_distal_hinges_on_trpb_and_gfp(self) -> None:
        from workers.spatial_filter import compute_anm_hinge_residues, get_spatial_contact_map

        hinges_trpb = compute_anm_hinge_residues(pdb_path=PDB_8EH0, cutoff_distance=12.0, top_k=10)
        cmap_trpb = get_spatial_contact_map(PDB_8EH0)

        self.assertEqual(len(hinges_trpb), 10)
        for res_id in hinges_trpb:
            self.assertIn(res_id, cmap_trpb.ca_coordinates)
            # Must not mutate the locked 4.2 A PLP cofactor catalytic core
            d_cof = cmap_trpb.get_cofactor_distance(res_id)
            if d_cof is not None:
                self.assertGreaterEqual(d_cof, 4.2)

        hinges_gfp = compute_anm_hinge_residues(pdb_path=PDB_1EMA, cutoff_distance=12.0, top_k=8)
        cmap_gfp = get_spatial_contact_map(PDB_1EMA)
        self.assertEqual(len(hinges_gfp), 8)
        for res_id in hinges_gfp:
            self.assertIn(res_id, cmap_gfp.ca_coordinates)
            self.assertNotIn(res_id, {65, 66, 67})


class TestSlice2ElasticBackboneRelaxation(unittest.TestCase):
    """Vertical Slice 2: Two-stage elastic mainchain accommodation (k_CA=500 -> k_CA=5.0)."""

    def test_two_stage_backbone_relaxation_accommodates_crowded_mutant_within_rmsd_cap(self) -> None:
        from workers.relaxation_worker import relax_variant_structure

        # Rigid Stage-1 only (allow_backbone_accommodation=False)
        res_rigid = relax_variant_structure(
            pdb_path=PDB_1EMA,
            mutations=["F64W:V68W"],
            allow_backbone_accommodation=False,
            strain_threshold=20.0,
        )
        # Two-stage with elastic backbone accommodation (allow_backbone_accommodation=True)
        res_elastic = relax_variant_structure(
            pdb_path=PDB_1EMA,
            mutations=["F64W:V68W"],
            allow_backbone_accommodation=True,
            strain_threshold=20.0,
        )

        self.assertTrue(res_elastic.backbone_stage2_triggered)
        self.assertFalse(res_rigid.backbone_stage2_triggered)
        self.assertLess(res_elastic.final_energy_kcal, res_rigid.final_energy_kcal)
        self.assertGreater(res_elastic.rmsd_ca_angstrom, res_rigid.rmsd_ca_angstrom)
        self.assertLessEqual(res_elastic.rmsd_ca_angstrom, 1.20)
        self.assertFalse(res_elastic.denaturation_rejected)


class TestSlice3RobustPlateOutlierRejection(unittest.TestCase):
    """Vertical Slice 3: Grubbs' test (alpha=0.05) & Huber M-estimation on multi-well replicates."""

    def test_parse_plate_csv_rejects_spiked_replicate_outlier_and_computes_variance(self) -> None:
        from scripts.ingest_plate_data import parse_plate_csv

        csv_content = "\n".join([
            "well,sample_name,sample_type,reading",
            "A01,Perimeter_Water,evaporation_barrier,10.0",
            "B02,WT,control_wt,1000.0",
            "B03,S65T:Y66H,variant,2480.0",
            "B04,S65T:Y66H,variant,2520.0",
            "B05,S65T:Y66H,variant,2500.0",
            "B06,S65T:Y66H,variant,2510.0",
            "H12,S65T:Y66H,variant,25000.0",  # 10x contaminated outlier well
            "C02,F64L,variant,1500.0",
            "C03,F64L,variant,1530.0",
            "C04,F64L,variant,1470.0",
        ])

        parsed = parse_plate_csv(csv_content, normalize_to_wt=True, reject_replicate_outliers=True)
        self.assertEqual(parsed["wt_baseline_score"], 1000.0)
        self.assertGreaterEqual(len(parsed["outliers_rejected"]), 1)
        rejected_wells = {r["well"] for r in parsed["outliers_rejected"]}
        self.assertIn("H12", rejected_wells)

        # Aggregated replicate record for S65T:Y66H should be ~2.5025 (not corrupted to ~7.0 by 25000.0)
        by_mut = {v["mutations"]: v for v in parsed["valid_variants"]}
        self.assertIn("S65T:Y66H", by_mut)
        s65t_rec = by_mut["S65T:Y66H"]
        self.assertAlmostEqual(s65t_rec["normalized_score"], 2.5025, delta=0.05)
        self.assertIsNotNone(s65t_rec["variance"])
        self.assertGreater(s65t_rec["variance"], 0.0)
        self.assertLess(s65t_rec["variance"], 0.01)
        self.assertIsNotNone(s65t_rec.get("standard_error"))


class TestSlice4AsyncKrigingLedger(unittest.TestCase):
    """Vertical Slice 4: Asynchronous Kriging-believer staging & physical reconciliation with FDA Part 11 Merkle chain."""

    def test_speculative_staging_and_physical_reconciliation_preserves_merkle_audit_chain(self) -> None:
        import tempfile
        from ledger import database as ledger_db

        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / "test_async_ledger.db"
            ledger_db.init_db(db_path)
            ledger_db.create_budget(total_allocation=5.00, db_path=db_path)

            # Fail-closed spend authorization required first
            job_id = "job-async-kriging-001"
            approved = ledger_db.authorize_spend(
                job_id=job_id,
                target_model="esm2_650m_kriging_believer",
                variant_count=2,
                estimated_cost=0.002500,
                db_path=db_path,
            )
            self.assertTrue(approved)

            # Stage provisional Kriging fantasy observations
            staged_count = ledger_db.stage_speculative_evidence(
                job_id=job_id,
                fantasy_results=[("MSKGEELFTG", 1.85), ("MSKGEELFTH", 2.10)],
                db_path=db_path,
            )
            self.assertEqual(staged_count, 2)

            job_after_stage = ledger_db.get_job(job_id, db_path=db_path)
            self.assertEqual(job_after_stage["status"], "STAGED_FOR_WETLAB")

            active_spec = ledger_db.get_speculative_evidence(job_id=job_id, active_only=True, db_path=db_path)
            self.assertEqual(len(active_spec), 2)

            # Transition to IN_FLIGHT_SYNTHESIS
            ledger_db.mark_job_in_flight(job_id=job_id, provider="strateos", external_run_id="run-991", db_path=db_path)
            job_in_flight = ledger_db.get_job(job_id, db_path=db_path)
            self.assertEqual(job_in_flight["status"], "IN_FLIGHT_SYNTHESIS")

            # Reconcile with physical plate reader observations
            reconciled = ledger_db.commit_physical_evidence_and_reconcile(
                job_id=job_id,
                observed_results=[("MSKGEELFTG", 1.92), ("MSKGEELFTH", 2.34)],
                db_path=db_path,
            )
            self.assertEqual(reconciled, 2)

            # Provisional rows are retired; permanent Evidence rows are committed; job is COMPLETED
            self.assertEqual(len(ledger_db.get_speculative_evidence(job_id=job_id, active_only=True, db_path=db_path)), 0)
            perm_ev = ledger_db.get_evidence(job_id=job_id, db_path=db_path)
            self.assertEqual(len(perm_ev), 2)
            self.assertAlmostEqual(perm_ev[1]["computed_score"], 2.34, places=4)
            self.assertEqual(ledger_db.get_job(job_id, db_path=db_path)["status"], "COMPLETED")

            # Verify FDA 21 CFR Part 11 Merkle hash chain is 100% intact
            is_valid, block_count, msg = ledger_db.verify_audit_integrity(db_path=db_path)
            self.assertTrue(is_valid, msg)
            self.assertGreaterEqual(block_count, 6)


class TestSlice5CloudLabWebhook(unittest.TestCase):
    """Vertical Slice 5: Autonomous Cloud-Lab SLL/Strateos dispatch & HMAC-authenticated webhook callback."""

    def test_cloudlab_dispatch_and_hmac_webhook_reconciles_in_flight_job(self) -> None:
        import json
        import tempfile
        from fastapi.testclient import TestClient
        from ledger import database as ledger_db
        from server import app
        from workers.cloudlab_worker import (
            EmeraldCloudLabClient,
            StrateosClient,
            compute_cloudlab_webhook_signature,
        )

        # 1. Verify SLL and Strateos protocol generation directly
        ecl = EmeraldCloudLabClient()
        sll_script = ecl.build_sll_protocol(
            campaign_id="camp-001",
            variants=["S65T", "S65T:Y66H"],
            assay_type="fluorescence",
        )
        self.assertIn("ExperimentRoboticLiquidHandling", sll_script)
        self.assertIn("ExperimentFluorescenceIntensity", sll_script)

        strateos = StrateosClient()
        autoprotocol = strateos.build_autoprotocol(
            campaign_id="camp-001",
            variants=["S65T", "S65T:Y66H"],
        )
        self.assertIn("instructions", autoprotocol)

        # 2. Verify REST dispatch -> IN_FLIGHT_SYNTHESIS -> HMAC webhook reconciliation
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / "test_cloudlab.db"
            ledger_db.init_db(db_path)
            ledger_db.create_budget(total_allocation=5.00, db_path=db_path)

            # The ledger path is server-side state: routes no longer accept a
            # db_path from the client, so point the server at this temp ledger.
            import server as server_mod
            saved_db = server_mod.DEFAULT_DB_PATH
            server_mod.DEFAULT_DB_PATH = db_path
            self.addCleanup(setattr, server_mod, "DEFAULT_DB_PATH", saved_db)

            client = TestClient(app)
            dispatch_res = client.post(
                "/api/cloudlab/dispatch",
                json={
                    "provider": "emerald_cloud_lab",
                    "target_accession": "P42212",
                    "variants": ["S65T", "F64L:S65T"],
                    "fantasy_scores": [1.75, 2.20],
                },
            )
            self.assertEqual(dispatch_res.status_code, 200, dispatch_res.text)
            dispatch_data = dispatch_res.json()
            job_id = dispatch_data["job_id"]
            self.assertEqual(dispatch_data["status"], "IN_FLIGHT_SYNTHESIS")
            self.assertIn("ExperimentRoboticLiquidHandling", dispatch_data["protocol_preview"])

            # Unsigned / forged webhook must be rejected (HTTP 401/403)
            webhook_payload = {
                "job_id": job_id,
                "target_accession": "P42212",
                "csv_data": "\n".join([
                    "well,sample_name,sample_type,reading",
                    "B02,WT,control_wt,1000.0",
                    "B03,S65T,variant,1850.0",
                    "B04,S65T,variant,1870.0",
                    "B05,S65T,variant,1860.0",
                    "H12,S65T,variant,19000.0",
                    "C02,F64L:S65T,variant,2420.0",
                    "C03,F64L:S65T,variant,2440.0",
                    "C04,F64L:S65T,variant,2430.0",
                ]),
            }
            raw_bytes = json.dumps(webhook_payload, sort_keys=True).encode("utf-8")

            bad_res = client.post(
                "/api/cloudlab/webhook",
                content=raw_bytes,
                headers={"Content-Type": "application/json", "X-CloudLab-Signature": "forged-sig"},
            )
            self.assertIn(bad_res.status_code, (401, 403))

            # Valid HMAC signature succeeds, drops outlier H12, reconciles ledger to COMPLETED
            valid_sig = compute_cloudlab_webhook_signature(raw_bytes)
            good_res = client.post(
                "/api/cloudlab/webhook",
                content=raw_bytes,
                headers={"Content-Type": "application/json", "X-CloudLab-Signature": valid_sig},
            )
            self.assertEqual(good_res.status_code, 200, good_res.text)
            good_data = good_res.json()
            self.assertEqual(good_data["job_status"], "COMPLETED")
            self.assertGreaterEqual(len(good_data["outliers_rejected"]), 1)

            is_valid, _, msg = ledger_db.verify_audit_integrity(db_path=db_path)
            self.assertTrue(is_valid, msg)


if __name__ == "__main__":
    unittest.main()

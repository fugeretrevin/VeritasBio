"""
tests/test_enterprise_compliance.py — Unit Tests for Gap 4 Enterprise Compliance & Benchling LIMS.

Validates:
    1. FDA 21 CFR Part 11 cryptographic hash chain:
       - Genesis block initialization on init_db().
       - Monotonic block numbering and SHA-256 parent hash chaining.
       - HMAC-SHA256 digital signature authenticity.
       - Complete lifecycle events (Budget, Job, Evidence, Completion, LIMS Sync).
       - Schema-level immutability triggers (UPDATE and DELETE prohibited).
       - Tamper detection in verify_audit_integrity() on data or signature modification.
    2. Benchling LIMS Integration:
       - DNA sequence registration with cDNA bases, annotations, and mock storage.
       - Virtual plate mapping with anti-edge evaporation layouts and container coordinates.
       - Standardized Bulk CSV manifest generation.
       - Audit block creation for LIMS events.
    3. Multi-Tenant PostgreSQL Schema DDL:
       - File existence and presence of RLS policies and immutability triggers.
"""

from __future__ import annotations

import csv
import io
import os
import sqlite3
import tempfile
import unittest
from pathlib import Path

from ledger import database as ledger_db
from workers.lims_worker import BenchlingClient, DEFAULT_PLATE_FORMAT


class TestEnterpriseCompliance(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.db_path = Path(self.temp_dir.name) / "test_compliance.db"
        ledger_db.init_db(self.db_path)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def test_genesis_block_creation(self) -> None:
        """Verifies that init_db atomically creates the FDA 21 CFR Part 11 Genesis Block."""
        blocks = ledger_db.list_audit_blocks(self.db_path)
        self.assertEqual(len(blocks), 1)

        genesis = blocks[0]
        self.assertEqual(genesis["block_index"], 0)
        self.assertEqual(genesis["action_type"], "GENESIS")
        self.assertEqual(genesis["prev_block_hash"], ledger_db.GENESIS_PREV_HASH)
        self.assertTrue(len(genesis["block_hash"]) == 64)
        self.assertTrue(len(genesis["signature"]) == 64)

        is_valid, count, msg = ledger_db.verify_audit_integrity(self.db_path)
        self.assertTrue(is_valid, msg)
        self.assertEqual(count, 1)

    def test_audit_chain_tracks_complete_lifecycle(self) -> None:
        """Verifies that budget, job authorization, evidence, and completion are hash-chained."""
        # 1. Budget creation
        ledger_db.create_budget(100.0, db_path=self.db_path)

        # 2. Spend authorization
        approved = ledger_db.authorize_spend(
            job_id="job_comp_001",
            target_model="esm2_650m",
            variant_count=10,
            estimated_cost=2.50,
            db_path=self.db_path,
        )
        self.assertTrue(approved)

        # 3. Evidence recording
        n_ev = ledger_db.record_evidence(
            job_id="job_comp_001",
            results=[("S65T", 2.45), ("Y66H", 1.80)],
            db_path=self.db_path,
        )
        self.assertEqual(n_ev, 2)

        # 4. Job completion
        ledger_db.mark_job_completed("job_comp_001", db_path=self.db_path)

        # Verify blocks: 0: GENESIS, 1: BUDGET_CREATED, 2: SPEND_AUTHORIZED, 3: EVIDENCE, 4: COMPLETED
        blocks = ledger_db.list_audit_blocks(self.db_path)
        self.assertEqual(len(blocks), 5)

        action_types = [b["action_type"] for b in blocks]
        self.assertEqual(
            action_types,
            [
                "GENESIS",
                "BUDGET_CREATED",
                "SPEND_AUTHORIZED",
                "EVIDENCE_BATCH_RECORDED",
                "JOB_COMPLETED",
            ],
        )

        # Cryptographic chain verification
        is_valid, count, msg = ledger_db.verify_audit_integrity(self.db_path)
        self.assertTrue(is_valid, msg)
        self.assertEqual(count, 5)

    def test_schema_level_immutability_triggers(self) -> None:
        """Verifies that SQLite triggers prohibit UPDATE or DELETE on AuditLedgerBlocks."""
        with ledger_db._connect(self.db_path) as conn:
            # Attempting to delete the genesis block must be aborted by trigger
            with self.assertRaises(sqlite3.DatabaseError) as ctx_del:
                conn.execute("DELETE FROM AuditLedgerBlocks WHERE block_index = 0;")
            self.assertIn("AuditLedgerBlocks cannot be deleted", str(ctx_del.exception))

            # Attempting to update the genesis block must be aborted by trigger
            with self.assertRaises(sqlite3.DatabaseError) as ctx_upd:
                conn.execute(
                    "UPDATE AuditLedgerBlocks SET payload_json = '{\"tampered\": true}' WHERE block_index = 0;"
                )
            self.assertIn("AuditLedgerBlocks are cryptographically immutable", str(ctx_upd.exception))

    def test_verify_audit_integrity_detects_hash_tampering(self) -> None:
        """Verifies that if block hash or payload is modified, verify_audit_integrity flags it."""
        # Create a budget and job
        ledger_db.create_budget(50.0, db_path=self.db_path)
        ledger_db.authorize_spend("job_tamper_test", "mock_esm", 5, 0.50, db_path=self.db_path)

        # Drop the immutability trigger temporarily to simulate an adversary with raw disk access
        with ledger_db._connect(self.db_path) as conn:
            conn.execute("DROP TRIGGER IF EXISTS trg_audit_no_update;")
            # Tamper with block 1 payload
            conn.execute(
                "UPDATE AuditLedgerBlocks SET payload_json = '{\"tampered\": 9999}' WHERE block_index = 1;"
            )

        # Now verify_audit_integrity must fail closed and detect the exact tampered block!
        is_valid, bad_idx, err_msg = ledger_db.verify_audit_integrity(self.db_path)
        self.assertFalse(is_valid)
        self.assertEqual(bad_idx, 1)
        self.assertIn("Data tampering detected at block 1", err_msg)

    def test_verify_audit_integrity_detects_signature_forgery(self) -> None:
        """Verifies that if an attacker modifies the HMAC digital signature, it is rejected."""
        ledger_db.create_budget(50.0, db_path=self.db_path)

        # Drop trigger to simulate raw disk write
        with ledger_db._connect(self.db_path) as conn:
            conn.execute("DROP TRIGGER IF EXISTS trg_audit_no_update;")
            conn.execute(
                "UPDATE AuditLedgerBlocks SET signature = '0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef' WHERE block_index = 1;"
            )

        is_valid, bad_idx, err_msg = ledger_db.verify_audit_integrity(self.db_path)
        self.assertFalse(is_valid)
        self.assertEqual(bad_idx, 1)
        self.assertIn("Digital signature forgery detected at block 1", err_msg)

    def test_benchling_client_mock_registration(self) -> None:
        """Verifies BenchlingClient registers cDNA constructs and records audit blocks."""
        client = BenchlingClient(mock_mode=True)

        res = client.register_dna_sequence(
            name="avGFP_S65T_Lead",
            bases="ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCCGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTT",
            annotations=[{"name": "Chromophore_S65T", "type": "mutation", "start": 65, "end": 65}],
            db_path=self.db_path,
        )

        self.assertTrue(res["id"].startswith("seq_"))
        self.assertEqual(res["name"], "avGFP_S65T_Lead")
        self.assertTrue(res["isMock"])

        # Check that audit block was created
        blocks = ledger_db.list_audit_blocks(self.db_path)
        last_block = blocks[-1]
        self.assertEqual(last_block["action_type"], "LIMS_SYNC")
        self.assertIn("avGFP_S65T_Lead", last_block["payload_json"])

        # Integrity check
        is_valid, _, msg = ledger_db.verify_audit_integrity(self.db_path)
        self.assertTrue(is_valid, msg)

    def test_benchling_client_virtual_plate_and_manifest_export(self) -> None:
        """Verifies virtual plate layout mapping and bulk CSV manifest generation."""
        client = BenchlingClient(mock_mode=True)

        # Mock plate layout with 4 wells
        plate_map = {
            "wells": {
                "B2": {"sample_name": "WT_Parent", "role": "wildtype", "volume_uL": 20.0, "predicted_fitness": 1.0, "developability_status": "PASS"},
                "B3": {"sample_name": "Neg_Control", "role": "negative_control", "volume_uL": 20.0, "predicted_fitness": 0.0, "developability_status": "PASS"},
                "B4": {"sample_name": "S65T:Y66H", "role": "experimental", "volume_uL": 20.0, "predicted_fitness": 2.45, "developability_status": "PASS"},
                "B5": {"sample_name": "F64L:S65T", "role": "experimental", "volume_uL": 20.0, "predicted_fitness": 2.89, "developability_status": "PASS"},
            }
        }

        plate_res = client.create_virtual_plate(
            plate_name="PanoptesBio_Run01_PlateA",
            plate_format=DEFAULT_PLATE_FORMAT,
            plate_map=plate_map,
            db_path=self.db_path,
        )

        self.assertTrue(plate_res["id"].startswith("plt_"))
        self.assertEqual(plate_res["containerCount"], 4)

        # Export manifest
        manifest_path = Path(self.temp_dir.name) / "benchling_manifest.csv"
        csv_str = client.export_lims_manifest(
            plate_name="PanoptesBio_Run01_PlateA",
            plate_map=plate_map,
            output_path=manifest_path,
        )

        self.assertTrue(manifest_path.exists())
        self.assertIn("Plate_Name,Well_Coordinate,Container_Barcode", csv_str)
        self.assertIn("PanoptesBio_Run01_PlateA,B4,CNT-PanoptesBio_Run01_PlateA-B4,S65T:Y66H,experimental,20.0,50.0,2.4500,PASS", csv_str)

    def test_postgres_schema_ddl_integrity(self) -> None:
        """Verifies ledger/postgres_schema.sql exists and contains valid DDL structures."""
        schema_path = Path(__file__).resolve().parent.parent / "ledger" / "postgres_schema.sql"
        self.assertTrue(schema_path.exists(), "postgres_schema.sql not found")

        content = schema_path.read_text(encoding="utf-8")
        self.assertIn("CREATE TABLE IF NOT EXISTS organizations", content)
        self.assertIn("CREATE TABLE IF NOT EXISTS budgets", content)
        self.assertIn("CREATE TABLE IF NOT EXISTS jobs", content)
        self.assertIn("CREATE TABLE IF NOT EXISTS evidence", content)
        self.assertIn("CREATE TABLE IF NOT EXISTS audit_ledger_blocks", content)
        self.assertIn("ENABLE ROW LEVEL SECURITY", content)
        self.assertIn("CREATE POLICY tenant_isolation", content)
        self.assertIn("trg_abort_immutable_operation", content)


if __name__ == "__main__":
    unittest.main()

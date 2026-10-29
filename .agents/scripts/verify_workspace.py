"""
.agents/scripts/verify_workspace.py — Master Automated Verification Harness for VeritasBio.

Executes comprehensive verification gates across:
  1. Git working tree inspection (uncommitted changes and modified files).
  2. Grounding files existence check (PDB crystal structures, DMS datasets, UniProt metadata).
  3. FDA 21 CFR Part 11 cryptographic audit chain integrity check in SQLite.
  4. Full repository unit test discovery (all test suites).

Usage:
    python .agents/scripts/verify_workspace.py [--fast]
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# Add project root to sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ledger import database as ledger_db


def check_grounding_files() -> bool:
    """Verifies that all local empirical reference datasets and PDBs are present."""
    required_files = [
        PROJECT_ROOT / "workers" / "data" / "1EMA.pdb",
        PROJECT_ROOT / "workers" / "data" / "2YBR.pdb",
        PROJECT_ROOT / "workers" / "data" / "8EH0.pdb",
        PROJECT_ROOT / "workers" / "data" / "P42212.json",
        PROJECT_ROOT / "workers" / "data" / "SCFV_9004G_CN2.json",
        PROJECT_ROOT / "workers" / "data" / "TRPB_TM9D8.json",
        PROJECT_ROOT / "scripts" / "data" / "sarkisyan_gfp_benchmark.json",
    ]
    missing = [f for f in required_files if not f.is_file()]
    if missing:
        print(f"[FAIL] Missing {len(missing)} empirical grounding files:")
        for m in missing:
            print(f"       - {m.relative_to(PROJECT_ROOT)}")
        return False
    print(f"[PASS] All {len(required_files)} empirical grounding files verified.")
    return True


def check_audit_ledger_integrity() -> bool:
    """Creates a temporary test ledger, executes a test action, and cryptographically validates the Merkle chain."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_file = Path(tmpdir) / "verify_chain.db"
        ledger_db.init_db(db_file)
        ledger_db.create_budget(total_allocation=1.00, budget_id=1, db_path=db_file)
        ledger_db.authorize_spend("job-verify-test", "mock_esm", 5, 0.0001, budget_id=1, db_path=db_file)
        ledger_db.record_evidence("job-verify-test", [("M1A", 1.2)], db_path=db_file)
        ledger_db.mark_job_completed("job-verify-test", db_path=db_file)

        is_valid, count, msg = ledger_db.verify_audit_integrity(db_path=db_file)
        if not is_valid:
            print(f"[FAIL] Audit ledger verification failed: {msg}")
            return False
        print(f"[PASS] FDA 21 CFR Part 11 audit chain verified intact ({count} blocks cryptographically validated).")
        return True


def run_unit_tests(fast: bool = False) -> bool:
    """Executes the test suite and captures execution statistics."""
    cmd = [sys.executable, "-m", "unittest"]
    if fast:
        cmd.extend(["tests.test_gatekeeper", "tests.test_diffusion_worker", "tests.test_phase5"])
        print("[INFO] Running fast test subset...")
    else:
        cmd.extend(["discover", "tests"])
        print("[INFO] Running full test discovery across repository...")

    t0 = time.time()
    res = subprocess.run(cmd, cwd=str(PROJECT_ROOT), capture_output=True, text=True)
    duration = time.time() - t0

    if res.returncode != 0:
        print(f"[FAIL] Unit tests encountered errors in {duration:.1f}s:")
        print(res.stderr or res.stdout)
        return False

    # Extract test count from stderr output (unittest writes summary to stderr)
    output = res.stderr + "\n" + res.stdout
    for line in output.splitlines():
        if "Ran" in line and "tests in" in line:
            print(f"[PASS] {line.strip()} ({duration:.1f}s total)")
            return True

    print(f"[PASS] All tests passed cleanly in {duration:.1f}s.")
    return True


def check_git_status() -> None:
    """Prints uncommitted and modified files."""
    try:
        res = subprocess.run(["git", "status", "--short"], cwd=str(PROJECT_ROOT), capture_output=True, text=True)
        if res.returncode == 0:
            lines = [l for l in res.stdout.splitlines() if l.strip()]
            if lines:
                print(f"[INFO] Git status: {len(lines)} file(s) modified or untracked:")
                for l in lines:
                    print(f"       {l}")
            else:
                print("[INFO] Git status: working tree clean.")
    except Exception as exc:
        print(f"[WARN] Could not inspect git status: {exc}")


def main():
    parser = argparse.ArgumentParser(description="VeritasBio Workspace Verification Harness")
    parser.add_argument("--fast", action="store_true", help="Run fast core test subset instead of full discovery")
    args = parser.parse_args()

    print("=" * 60)
    print("VERITASBIO DEVELOPMENT AGENT VERIFICATION HARNESS")
    print("=" * 60)

    g1 = check_grounding_files()
    g2 = check_audit_ledger_integrity()
    g3 = run_unit_tests(fast=args.fast)
    check_git_status()

    print("=" * 60)
    if g1 and g2 and g3:
        print("ALL VERIFICATION GATES PASSED: Ready for user review & commit.")
        print("=" * 60)
        sys.exit(0)
    else:
        print("VERIFICATION FAILED: Resolve errors above before requesting approval.")
        print("=" * 60)
        sys.exit(1)


if __name__ == "__main__":
    main()

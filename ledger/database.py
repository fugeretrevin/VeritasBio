"""
ledger/database.py — DOMAIN: Immutable Memory.

SQLite CRUD wrapper for the VeritasBio ledger (ARCHITECTURE.md §5, §6).

Design rules:
    * authorize_spend() is the ONLY function that can approve a job or move money
      (ARCHITECTURE.md §7.2). It is fail-closed: it returns True only when every
      check passes and the deduction + APPROVED job record were committed
      atomically. Every other outcome returns False and it never raises.
    * The ledger is append-only. There are no public update/delete functions for
      Budget spend, Jobs, or Evidence; the schema triggers enforce the same rule.
    * Deterministic: no randomness, no LLM involvement, no network calls.

Money units:
    The public API speaks float USD. Internally (and in SQLite) money is stored as
    INTEGER micro-USD (1 USD = 1,000,000) for exact arithmetic. Conversions are
    conservative: costs round UP, allocations round DOWN, so rounding can never
    let a job overspend.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import logging
import math
import os
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple, Union

logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
LEDGER_DIR: Path = Path(__file__).resolve().parent
SCHEMA_PATH: Path = LEDGER_DIR / "schema.sql"
DEFAULT_DB_PATH: Path = LEDGER_DIR / "veritasbio_ledger.db"
DEFAULT_BUDGET_ID: int = 1
_BUSY_TIMEOUT_SECONDS: float = 5.0

MICRO_USD_PER_USD: int = 1_000_000
# Sanity ceiling well inside SQLite's signed 64-bit INTEGER range (~$4.6 trillion).
_MAX_MICRO_USD: int = 2**62

STATUS_PENDING = "PENDING"
STATUS_APPROVED = "APPROVED"
STATUS_REJECTED = "REJECTED"
STATUS_COMPLETED = "COMPLETED"
STATUS_FAILED = "FAILED"
STATUS_STAGED_FOR_WETLAB = "STAGED_FOR_WETLAB"
STATUS_IN_FLIGHT_SYNTHESIS = "IN_FLIGHT_SYNTHESIS"
JOB_STATUSES = frozenset(
    {
        STATUS_PENDING,
        STATUS_APPROVED,
        STATUS_REJECTED,
        STATUS_COMPLETED,
        STATUS_FAILED,
        STATUS_STAGED_FOR_WETLAB,
        STATUS_IN_FLIGHT_SYNTHESIS,
    }
)

HMAC_AUDIT_SECRET: bytes = os.environ.get(
    "VERITAS_AUDIT_SECRET", "veritasbio-fda-21cfr11-genesis-secret-key"
).encode("utf-8")
GENESIS_PREV_HASH: str = "0" * 64

DbPath = Union[str, Path]


class LedgerError(Exception):
    """Raised by non-authorization ledger operations that cannot be completed."""


class _SpendRejected(Exception):
    """Internal: forces a rollback inside authorize_spend()."""


# -----------------------------------------------------------------------------
# Connection / transaction helpers
# -----------------------------------------------------------------------------
@contextmanager
def _connect(db_path: Optional[DbPath] = None) -> Iterator[sqlite3.Connection]:
    """Open a connection with FK enforcement and explicit transaction control."""
    conn = sqlite3.connect(
        str(db_path or DEFAULT_DB_PATH),
        timeout=_BUSY_TIMEOUT_SECONDS,
        isolation_level=None,  # autocommit; transactions are opened explicitly
    )
    try:
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON;")
        conn.execute("PRAGMA journal_mode = WAL;")
        conn.execute("PRAGMA busy_timeout = 5000;")
        yield conn
    finally:
        conn.close()


@contextmanager
def _transaction(conn: sqlite3.Connection) -> Iterator[None]:
    """BEGIN IMMEDIATE takes the write lock up front, so check-then-deduct is atomic."""
    conn.execute("BEGIN IMMEDIATE;")
    try:
        yield
    except BaseException:
        if conn.in_transaction:
            conn.execute("ROLLBACK;")
        raise
    else:
        conn.execute("COMMIT;")


# -----------------------------------------------------------------------------
# Validation & money helpers
# -----------------------------------------------------------------------------
def _is_non_empty_str(value: Any) -> bool:
    return isinstance(value, str) and value.strip() != ""


def _is_positive_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


def _is_finite_number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
    )


def _is_positive_finite_number(value: Any) -> bool:
    return _is_finite_number(value) and value > 0


def _usd_to_micro(value_usd: float, *, round_up: bool) -> int:
    """
    Convert float USD to integer micro-USD.

    round_up=True  for costs       (never undercharge)
    round_up=False for allocations (never over-grant)
    Raises ValueError if the result exceeds the ledger's sanity ceiling.
    """
    rounding = ROUND_CEILING if round_up else ROUND_FLOOR
    micro = int((Decimal(str(value_usd)) * MICRO_USD_PER_USD).to_integral_value(rounding=rounding))
    if micro > _MAX_MICRO_USD:
        raise ValueError(f"amount ${value_usd} exceeds ledger maximum")
    return micro


def _micro_to_usd(micro: int) -> float:
    return micro / MICRO_USD_PER_USD


def _job_row_to_dict(row: sqlite3.Row) -> Dict[str, Any]:
    job = dict(row)
    job["estimated_cost"] = _micro_to_usd(job["estimated_cost"])
    return job


def _compute_block_hash(
    block_index: int,
    timestamp: str,
    action_type: str,
    payload_json: str,
    prev_block_hash: str,
) -> str:
    raw = f"{block_index}|{timestamp}|{action_type}|{payload_json}|{prev_block_hash}".encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _compute_block_signature(block_hash: str) -> str:
    return hmac.new(HMAC_AUDIT_SECRET, block_hash.encode("utf-8"), hashlib.sha256).hexdigest()


def _append_audit_block_tx(
    conn: sqlite3.Connection,
    action_type: str,
    payload: Dict[str, Any],
) -> Dict[str, Any]:
    """Appends an immutable, cryptographic hash-linked block to AuditLedgerBlocks within the active transaction."""
    row = conn.execute(
        "SELECT block_index, block_hash FROM AuditLedgerBlocks ORDER BY block_index DESC LIMIT 1;"
    ).fetchone()

    if row is None:
        block_index = 0
        prev_block_hash = GENESIS_PREV_HASH
    else:
        block_index = int(row["block_index"]) + 1
        prev_block_hash = str(row["block_hash"])

    timestamp = datetime.now(timezone.utc).isoformat()
    payload_json = json.dumps(payload, sort_keys=True)
    block_hash = _compute_block_hash(block_index, timestamp, action_type, payload_json, prev_block_hash)
    signature = _compute_block_signature(block_hash)

    conn.execute(
        "INSERT INTO AuditLedgerBlocks (block_index, timestamp, action_type, payload_json, prev_block_hash, block_hash, signature) "
        "VALUES (?, ?, ?, ?, ?, ?, ?);",
        (block_index, timestamp, action_type, payload_json, prev_block_hash, block_hash, signature),
    )

    return {
        "block_index": block_index,
        "timestamp": timestamp,
        "action_type": action_type,
        "prev_block_hash": prev_block_hash,
        "block_hash": block_hash,
        "signature": signature,
    }


# -----------------------------------------------------------------------------
# Initialization
# -----------------------------------------------------------------------------
def _migrate_evidence_provenance_tx(conn: sqlite3.Connection) -> None:
    """
    Add the Evidence provenance columns to a ledger created before they existed.

    `CREATE TABLE IF NOT EXISTS` will not alter an existing table, so a ledger
    written by an earlier version keeps the old 5-column Evidence shape. Rows
    already present get the conservative defaults ('unknown' / is_mock = 1):
    their true origin was never recorded, so they must not read as verified.
    Idempotent — safe to run on every init.
    """
    existing = {row["name"] for row in conn.execute("PRAGMA table_info(Evidence);")}
    if not existing:
        return
    if "compute_engine" not in existing:
        conn.execute(
            "ALTER TABLE Evidence ADD COLUMN compute_engine TEXT NOT NULL DEFAULT 'unknown';"
        )
        logger.info("Migrated Evidence: added compute_engine column")
    if "is_mock" not in existing:
        conn.execute("ALTER TABLE Evidence ADD COLUMN is_mock INTEGER NOT NULL DEFAULT 1;")
        logger.info("Migrated Evidence: added is_mock column")


def _migrate_audit_action_types_tx(conn: sqlite3.Connection) -> None:
    """
    Widen the AuditLedgerBlocks.action_type CHECK constraint on an older ledger.

    SQLite cannot ALTER a CHECK constraint, so a ledger created before a new block
    type existed will reject it. The table is rebuilt with the current constraint and
    every row copied verbatim - including block_hash, prev_block_hash and signature -
    so the hash chain is preserved exactly and `verify_audit_integrity` still passes.
    Nothing is recomputed; this changes only the constraint, never the data.

    Idempotent: returns immediately when the constraint already permits the type.
    """
    row = conn.execute(
        "SELECT sql FROM sqlite_master WHERE type = 'table' AND name = 'AuditLedgerBlocks';"
    ).fetchone()
    if row is None or not row["sql"]:
        return  # fresh database; schema.sql will create it with the current constraint
    if "STRATEGY_PLANNED" in row["sql"]:
        return  # already current

    before = conn.execute("SELECT COUNT(*) AS n FROM AuditLedgerBlocks;").fetchone()["n"]
    logger.info("Migrating AuditLedgerBlocks constraint (%d existing blocks preserved)", before)

    conn.executescript(
        """
        PRAGMA legacy_alter_table = OFF;
        ALTER TABLE AuditLedgerBlocks RENAME TO _AuditLedgerBlocks_old;
        """
    )
    # schema.sql (executed right after this in init_db) recreates the table with the
    # current constraint; copy the old rows back and drop the temporary table.
    conn.executescript(SCHEMA_PATH.read_text(encoding="utf-8"))
    # Copy by introspected column intersection rather than a hardcoded list, so the
    # migration cannot drift from the schema (the audit timestamp column is named
    # `timestamp`, and guessing it wrong fails the whole init).
    old_cols = [r["name"] for r in conn.execute("PRAGMA table_info(_AuditLedgerBlocks_old);")]
    new_cols = [r["name"] for r in conn.execute("PRAGMA table_info(AuditLedgerBlocks);")]
    shared = [c for c in new_cols if c in old_cols]
    required = {"block_index", "action_type", "payload_json", "prev_block_hash", "block_hash"}
    missing = required - set(shared)
    if missing:
        raise LedgerError(
            f"Audit migration cannot proceed: columns {sorted(missing)} absent from the "
            f"existing table. Refusing to rebuild and risk the hash chain."
        )
    col_list = ", ".join(shared)
    conn.execute(
        f"INSERT INTO AuditLedgerBlocks ({col_list}) "
        f"SELECT {col_list} FROM _AuditLedgerBlocks_old ORDER BY block_index;"
    )
    conn.execute("DROP TABLE _AuditLedgerBlocks_old;")
    after = conn.execute("SELECT COUNT(*) AS n FROM AuditLedgerBlocks;").fetchone()["n"]
    if after != before:
        raise LedgerError(
            f"Audit migration lost blocks: {before} before, {after} after. Refusing to continue."
        )
    logger.info("AuditLedgerBlocks constraint migrated; %d blocks intact", after)


def init_db(db_path: Optional[DbPath] = None) -> Path:
    """Create the ledger file (if needed) and apply schema.sql. Idempotent."""
    path = Path(db_path) if db_path else DEFAULT_DB_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    schema_sql = SCHEMA_PATH.read_text(encoding="utf-8")
    with _connect(path) as conn:
        # Widen the audit action_type constraint BEFORE applying schema.sql: the
        # migration rebuilds the table and re-applies the schema itself, so running it
        # first keeps a single code path for both fresh and pre-existing ledgers.
        _migrate_audit_action_types_tx(conn)
        conn.executescript(schema_sql)
        _migrate_evidence_provenance_tx(conn)
        with _transaction(conn):
            # Ensure FDA 21 CFR Part 11 Genesis Block is present
            row = conn.execute("SELECT COUNT(*) as cnt FROM AuditLedgerBlocks;").fetchone()
            if row and row["cnt"] == 0:
                _append_audit_block_tx(
                    conn,
                    action_type="GENESIS",
                    payload={"system": "VeritasBio", "standard": "FDA 21 CFR Part 11", "version": "1.0"},
                )
    logger.info("Ledger initialized at %s", path)
    return path


# -----------------------------------------------------------------------------
# Budget
# -----------------------------------------------------------------------------
def create_budget(
    total_allocation: float,
    budget_id: int = DEFAULT_BUDGET_ID,
    db_path: Optional[DbPath] = None,
) -> None:
    """Create a budget row (USD) with zero spend. Raises if it already exists."""
    if not _is_positive_finite_number(total_allocation):
        raise ValueError("total_allocation must be a positive, finite number")
    if not _is_positive_int(budget_id):
        raise ValueError("budget_id must be a positive integer")
    allocation_micro = _usd_to_micro(total_allocation, round_up=False)
    if allocation_micro <= 0:
        raise ValueError("total_allocation is below the ledger's smallest unit ($0.000001)")
    with _connect(db_path) as conn, _transaction(conn):
        try:
            conn.execute(
                "INSERT INTO Budget (id, total_allocation, cumulative_spent) VALUES (?, ?, 0);",
                (budget_id, allocation_micro),
            )
            _append_audit_block_tx(
                conn,
                action_type="BUDGET_CREATED",
                payload={
                    "budget_id": budget_id,
                    "total_allocation_usd": total_allocation,
                    "allocation_micro_usd": allocation_micro,
                },
            )
        except sqlite3.IntegrityError as exc:
            raise LedgerError(f"Budget {budget_id} already exists") from exc


def get_budget(
    budget_id: int = DEFAULT_BUDGET_ID,
    db_path: Optional[DbPath] = None,
) -> Optional[Dict[str, Any]]:
    """Return {id, total_allocation, cumulative_spent, remaining} in USD, or None."""
    with _connect(db_path) as conn:
        row = conn.execute(
            "SELECT id, total_allocation, cumulative_spent FROM Budget WHERE id = ?;",
            (budget_id,),
        ).fetchone()
    if row is None:
        return None
    return {
        "id": row["id"],
        "total_allocation": _micro_to_usd(row["total_allocation"]),
        "cumulative_spent": _micro_to_usd(row["cumulative_spent"]),
        "remaining": _micro_to_usd(row["total_allocation"] - row["cumulative_spent"]),
    }


def get_remaining_budget(
    budget_id: int = DEFAULT_BUDGET_ID,
    db_path: Optional[DbPath] = None,
) -> float:
    """Remaining USD. Fail-closed: a missing budget reports 0.0 (nothing spendable)."""
    budget = get_budget(budget_id, db_path)
    return float(budget["remaining"]) if budget else 0.0


# -----------------------------------------------------------------------------
# Authorization (the single spend gate)
# -----------------------------------------------------------------------------
def authorize_spend(
    job_id: str,
    target_model: str,
    variant_count: int,
    estimated_cost: float,
    budget_id: int = DEFAULT_BUDGET_ID,
    db_path: Optional[DbPath] = None,
    operator_uid: Optional[str] = None,
    operator_email: Optional[str] = None,
) -> bool:
    """
    Approve or hard-reject a compute job against the Budget table.

    estimated_cost is in USD; it is rounded UP to the nearest micro-USD.

    Within one IMMEDIATE transaction this:
        1. Reads the Budget row (missing row -> reject).
        2. Compares the cost to remaining funds (exact integer math).
        3a. Sufficient: atomically increments cumulative_spent (guarded by a
            conditional UPDATE) and inserts the job as APPROVED -> returns True.
        3b. Insufficient: inserts the job as REJECTED (audit trail) -> returns False.

    Fail-closed: invalid input, a missing/zero/non-finite cost estimate, a missing
    budget, a duplicate job_id, any database error, or any unexpected exception
    returns False. No funds move unless True is returned. Never raises
    (ordinary exceptions are caught and logged).
    """
    try:
        # --- Input validation (no DB access on malformed payloads) ---
        if not _is_non_empty_str(job_id):
            logger.warning("authorize_spend REJECTED: invalid job_id %r", job_id)
            return False
        if not _is_non_empty_str(target_model):
            logger.warning("authorize_spend REJECTED job %r: invalid target_model", job_id)
            return False
        if not _is_positive_int(variant_count):
            logger.warning("authorize_spend REJECTED job %r: invalid variant_count %r",
                           job_id, variant_count)
            return False
        if not _is_positive_finite_number(estimated_cost):
            logger.warning("authorize_spend REJECTED job %r: missing or invalid cost estimate %r",
                           job_id, estimated_cost)
            return False
        if not _is_positive_int(budget_id):
            logger.warning("authorize_spend REJECTED job %r: invalid budget_id %r",
                           job_id, budget_id)
            return False

        cost_micro = _usd_to_micro(estimated_cost, round_up=True)  # >= 1 for any cost > 0

        with _connect(db_path) as conn, _transaction(conn):
            # --- Explicit Budget check ---
            budget = conn.execute(
                "SELECT total_allocation, cumulative_spent FROM Budget WHERE id = ?;",
                (budget_id,),
            ).fetchone()
            if budget is None:
                logger.warning("authorize_spend REJECTED job %r: budget %r does not exist",
                               job_id, budget_id)
                return False

            remaining_micro = budget["total_allocation"] - budget["cumulative_spent"]
            if cost_micro > remaining_micro:
                conn.execute(
                    "INSERT INTO Jobs (job_id, target_model, variant_count, estimated_cost, status) "
                    "VALUES (?, ?, ?, ?, ?);",
                    (job_id, target_model, variant_count, cost_micro, STATUS_REJECTED),
                )
                logger.warning(
                    "authorize_spend REJECTED job %r: cost $%.6f exceeds remaining $%.6f",
                    job_id, _micro_to_usd(cost_micro), _micro_to_usd(remaining_micro),
                )
                return False  # REJECTED row is committed for the audit trail

            # --- Deduct (conditional UPDATE re-asserts the budget invariant) ---
            cur = conn.execute(
                "UPDATE Budget SET cumulative_spent = cumulative_spent + ? "
                "WHERE id = ? AND cumulative_spent + ? <= total_allocation;",
                (cost_micro, budget_id, cost_micro),
            )
            if cur.rowcount != 1:
                raise _SpendRejected("conditional budget deduction did not apply")

            conn.execute(
                "INSERT INTO Jobs (job_id, target_model, variant_count, estimated_cost, status) "
                "VALUES (?, ?, ?, ?, ?);",
                (job_id, target_model, variant_count, cost_micro, STATUS_APPROVED),
            )
            spend_payload = {
                "job_id": job_id,
                "target_model": target_model,
                "variant_count": variant_count,
                "cost_micro_usd": cost_micro,
                "cost_usd": _micro_to_usd(cost_micro),
                "budget_id": budget_id,
            }
            if operator_uid:
                spend_payload["operator_uid"] = operator_uid
            if operator_email:
                spend_payload["operator_email"] = operator_email

            _append_audit_block_tx(
                conn,
                action_type="SPEND_AUTHORIZED",
                payload=spend_payload,
            )

        logger.info("authorize_spend APPROVED job %r: $%.6f for %d variants on %s",
                    job_id, _micro_to_usd(cost_micro), variant_count, target_model)
        return True

    except Exception:  # fail closed on anything unexpected; rollback already done
        logger.exception("authorize_spend REJECTED job %r: error during authorization", job_id)
        return False


# -----------------------------------------------------------------------------
# Jobs
# -----------------------------------------------------------------------------
def _transition_approved_job(job_id: str, new_status: str, db_path: Optional[DbPath]) -> None:
    if not _is_non_empty_str(job_id):
        raise ValueError("job_id must be a non-empty string")
    with _connect(db_path) as conn, _transaction(conn):
        cur = conn.execute(
            "UPDATE Jobs SET status = ? WHERE job_id = ? AND status = ?;",
            (new_status, job_id, STATUS_APPROVED),
        )
        if cur.rowcount != 1:
            raise LedgerError(f"Job {job_id!r} is not APPROVED (or does not exist)")
        _append_audit_block_tx(
            conn,
            action_type="JOB_COMPLETED" if new_status == STATUS_COMPLETED else "JOB_FAILED",
            payload={"job_id": job_id, "new_status": new_status},
        )


def mark_job_completed(job_id: str, db_path: Optional[DbPath] = None) -> None:
    """Transition an APPROVED job to COMPLETED. Raises LedgerError otherwise."""
    _transition_approved_job(job_id, STATUS_COMPLETED, db_path)


def mark_job_failed(job_id: str, db_path: Optional[DbPath] = None) -> None:
    """
    Transition an APPROVED job to FAILED (e.g. Modal execution error).
    Funds are NOT refunded: the compute may already have been billed.
    Raises LedgerError if the job is not APPROVED.
    """
    _transition_approved_job(job_id, STATUS_FAILED, db_path)


def get_job(job_id: str, db_path: Optional[DbPath] = None) -> Optional[Dict[str, Any]]:
    """Return a single job (estimated_cost in USD), or None."""
    with _connect(db_path) as conn:
        row = conn.execute("SELECT * FROM Jobs WHERE job_id = ?;", (job_id,)).fetchone()
    return _job_row_to_dict(row) if row is not None else None


def list_jobs(
    status: Optional[str] = None,
    db_path: Optional[DbPath] = None,
) -> List[Dict[str, Any]]:
    """List jobs (estimated_cost in USD), optionally filtered by status, in insertion order."""
    if status is not None and status not in JOB_STATUSES:
        raise ValueError(f"Unknown job status {status!r}")
    with _connect(db_path) as conn:
        if status is None:
            rows = conn.execute("SELECT * FROM Jobs ORDER BY rowid;").fetchall()
        else:
            rows = conn.execute(
                "SELECT * FROM Jobs WHERE status = ? ORDER BY rowid;", (status,)
            ).fetchall()
    return [_job_row_to_dict(r) for r in rows]


# -----------------------------------------------------------------------------
# Evidence
# -----------------------------------------------------------------------------
def infer_is_mock(compute_engine: str) -> bool:
    """
    Decide whether a compute_engine label denotes synthetic (non-GPU) output.

    Fail-closed: anything not recognised as a real inference engine is treated as
    mock. An undeclared or unfamiliar origin must never read as verified science.
    """
    label = (compute_engine or "").strip().lower()
    return not (label.startswith("esm2_") or label == "esm_cascade")


def record_evidence(
    job_id: str,
    results: Iterable[Tuple[str, float]],
    db_path: Optional[DbPath] = None,
    compute_engine: str = "unknown",
    is_mock: Optional[bool] = None,
) -> int:
    """
    Atomically append (variant_sequence, computed_score) results for an APPROVED job.

    All-or-nothing: if any row is invalid or the job is not APPROVED, nothing is
    written and an exception is raised. Returns the number of rows inserted.

    `compute_engine` records which engine actually produced these scores (e.g.
    'esm2_650m', 'mock_esm', 'mock_esm_fallback_unconfigured'). `is_mock` is
    derived from it via `infer_is_mock` unless passed explicitly. Both are stored
    on every row and sealed into the EVIDENCE_BATCH_RECORDED audit block, so the
    Part 11 chain distinguishes real inference from synthetic output.
    """
    if not _is_non_empty_str(job_id):
        raise ValueError("job_id must be a non-empty string")
    if not _is_non_empty_str(compute_engine):
        raise ValueError("compute_engine must be a non-empty string")

    engine = compute_engine.strip()
    mock_flag = infer_is_mock(engine) if is_mock is None else bool(is_mock)

    rows: List[Tuple[str, str, float, str, int]] = []
    for sequence, score in results:
        if not _is_non_empty_str(sequence):
            raise ValueError("variant_sequence must be a non-empty string")
        if not _is_finite_number(score):
            raise ValueError(f"computed_score must be finite, got {score!r}")
        rows.append((job_id, sequence, float(score), engine, int(mock_flag)))
    if not rows:
        raise ValueError("results must contain at least one (sequence, score) pair")

    with _connect(db_path) as conn, _transaction(conn):
        try:
            conn.executemany(
                "INSERT INTO Evidence (job_id, variant_sequence, computed_score, "
                "compute_engine, is_mock) VALUES (?, ?, ?, ?, ?);",
                rows,
            )
            # Append cryptographic audit block for the recorded evidence batch
            _append_audit_block_tx(
                conn,
                action_type="EVIDENCE_BATCH_RECORDED",
                payload={
                    "job_id": job_id,
                    "variant_count": len(rows),
                    "sample_variants": [r[1] for r in rows[:5]],
                    "compute_engine": engine,
                    "is_mock": bool(mock_flag),
                },
            )
        except sqlite3.DatabaseError as exc:
            raise LedgerError(f"Could not record evidence for job {job_id!r}: {exc}") from exc
    return len(rows)


def get_evidence(
    job_id: Optional[str] = None,
    db_path: Optional[DbPath] = None,
) -> List[Dict[str, Any]]:
    """Return evidence rows (all, or for one job) in insertion order."""
    with _connect(db_path) as conn:
        if job_id is None:
            rows = conn.execute("SELECT * FROM Evidence ORDER BY evidence_id;").fetchall()
        else:
            rows = conn.execute(
                "SELECT * FROM Evidence WHERE job_id = ? ORDER BY evidence_id;", (job_id,)
            ).fetchall()
    return [dict(r) for r in rows]


# -----------------------------------------------------------------------------
# Asynchronous Kriging-Believer Staging & Physical Reconciliation (SPEC-001 Seam 4)
# -----------------------------------------------------------------------------
def stage_speculative_evidence(
    job_id: str,
    fantasy_results: Iterable[Tuple[str, float]],
    db_path: Optional[DbPath] = None,
) -> int:
    """
    Stages provisional Kriging-believer fantasy observations for an APPROVED job and
    transitions its status to STAGED_FOR_WETLAB while preserving the FDA Part 11 Merkle chain.
    """
    if not _is_non_empty_str(job_id):
        raise ValueError("job_id must be a non-empty string")

    rows: List[Tuple[str, str, float, int]] = []
    for sequence, fantasy_score in fantasy_results:
        if not _is_non_empty_str(sequence):
            raise ValueError("variant_sequence must be a non-empty string")
        if not _is_finite_number(fantasy_score):
            raise ValueError(f"fantasy_score must be finite, got {fantasy_score!r}")
        rows.append((job_id, sequence, float(fantasy_score), 1))
    if not rows:
        raise ValueError("fantasy_results must contain at least one (sequence, fantasy_score) pair")

    with _connect(db_path) as conn, _transaction(conn):
        cur = conn.execute(
            "UPDATE Jobs SET status = ? WHERE job_id = ? AND status = ?;",
            (STATUS_STAGED_FOR_WETLAB, job_id, STATUS_APPROVED),
        )
        if cur.rowcount != 1:
            raise LedgerError(f"Job {job_id!r} must be APPROVED before staging speculative evidence")

        conn.executemany(
            "INSERT INTO SpeculativeEvidence (job_id, variant_sequence, fantasy_score, is_provisional) "
            "VALUES (?, ?, ?, ?);",
            rows,
        )
        _append_audit_block_tx(
            conn,
            action_type="SPECULATIVE_BATCH_STAGED",
            payload={
                "job_id": job_id,
                "new_status": STATUS_STAGED_FOR_WETLAB,
                "provisional_count": len(rows),
                "sample_variants": [r[1] for r in rows[:5]],
            },
        )
    return len(rows)


def mark_job_in_flight(
    job_id: str,
    provider: str = "strateos",
    external_run_id: Optional[str] = None,
    db_path: Optional[DbPath] = None,
    simulated: bool = True,
) -> None:
    """
    Transitions an APPROVED or STAGED_FOR_WETLAB job to IN_FLIGHT_SYNTHESIS upon
    cloud-lab protocol dispatch.

    `simulated` records whether the dispatch actually reached a remote provider.
    It defaults to True and is sealed into the audit block, so a locally minted
    run identifier can never be mistaken for a real cloud-lab run.
    """
    if not _is_non_empty_str(job_id):
        raise ValueError("job_id must be a non-empty string")
    with _connect(db_path) as conn, _transaction(conn):
        cur = conn.execute(
            "UPDATE Jobs SET status = ? WHERE job_id = ? AND status IN (?, ?);",
            (STATUS_IN_FLIGHT_SYNTHESIS, job_id, STATUS_APPROVED, STATUS_STAGED_FOR_WETLAB),
        )
        if cur.rowcount != 1:
            raise LedgerError(f"Job {job_id!r} is not in APPROVED or STAGED_FOR_WETLAB state")
        _append_audit_block_tx(
            conn,
            action_type="IN_FLIGHT_DISPATCHED",
            payload={
                "job_id": job_id,
                "new_status": STATUS_IN_FLIGHT_SYNTHESIS,
                "provider": provider,
                "external_run_id": external_run_id or "",
                "simulated": bool(simulated),
            },
        )


def get_speculative_evidence(
    job_id: Optional[str] = None,
    active_only: bool = True,
    db_path: Optional[DbPath] = None,
) -> List[Dict[str, Any]]:
    """Return speculative Kriging-believer rows, optionally filtered by job_id and active provisional state."""
    with _connect(db_path) as conn:
        clauses: List[str] = []
        params: List[Any] = []
        if job_id is not None:
            clauses.append("job_id = ?")
            params.append(job_id)
        if active_only:
            clauses.append("is_provisional = 1")
        where_sql = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        rows = conn.execute(
            f"SELECT * FROM SpeculativeEvidence {where_sql} ORDER BY speculative_id;",
            tuple(params),
        ).fetchall()
    return [dict(r) for r in rows]


def commit_physical_evidence_and_reconcile(
    job_id: str,
    observed_results: Iterable[Tuple[str, float]],
    db_path: Optional[DbPath] = None,
) -> int:
    """
    Atomically retires provisional Kriging-believer rows (`is_provisional = 0`),
    records verified physical Evidence rows, transitions the job to COMPLETED,
    and appends immutable FDA 21 CFR Part 11 Merkle blocks.
    """
    if not _is_non_empty_str(job_id):
        raise ValueError("job_id must be a non-empty string")

    rows: List[Tuple[str, str, float]] = []
    for sequence, score in observed_results:
        if not _is_non_empty_str(sequence):
            raise ValueError("variant_sequence must be a non-empty string")
        if not _is_finite_number(score):
            raise ValueError(f"computed_score must be finite, got {score!r}")
        rows.append((job_id, sequence, float(score)))
    if not rows:
        raise ValueError("observed_results must contain at least one (sequence, score) pair")

    with _connect(db_path) as conn, _transaction(conn):
        conn.execute(
            "UPDATE SpeculativeEvidence SET is_provisional = 0 WHERE job_id = ?;",
            (job_id,),
        )
        conn.executemany(
            "INSERT INTO Evidence (job_id, variant_sequence, computed_score) VALUES (?, ?, ?);",
            rows,
        )
        cur = conn.execute(
            "UPDATE Jobs SET status = ? WHERE job_id = ? AND status IN (?, ?, ?);",
            (
                STATUS_COMPLETED,
                job_id,
                STATUS_APPROVED,
                STATUS_STAGED_FOR_WETLAB,
                STATUS_IN_FLIGHT_SYNTHESIS,
            ),
        )
        if cur.rowcount != 1:
            raise LedgerError(f"Job {job_id!r} cannot be reconciled (invalid status or missing job)")

        _append_audit_block_tx(
            conn,
            action_type="PHYSICAL_EVIDENCE_RECONCILED",
            payload={
                "job_id": job_id,
                "reconciled_count": len(rows),
                "sample_variants": [r[1] for r in rows[:5]],
            },
        )
        _append_audit_block_tx(
            conn,
            action_type="JOB_COMPLETED",
            payload={"job_id": job_id, "new_status": STATUS_COMPLETED},
        )
    return len(rows)


# -----------------------------------------------------------------------------
# FDA 21 CFR Part 11 Audit Ledger Verification & Operations
# -----------------------------------------------------------------------------
def record_lims_sync_audit(
    entity_type: str,
    entity_id: str,
    benchling_id: str,
    details: Optional[Dict[str, Any]] = None,
    db_path: Optional[DbPath] = None,
) -> Dict[str, Any]:
    """
    Record an FDA 21 CFR Part 11 compliant audit block for a Benchling LIMS sync event.
    """
    if not _is_non_empty_str(entity_type):
        raise ValueError("entity_type must be a non-empty string")
    if not _is_non_empty_str(entity_id):
        raise ValueError("entity_id must be a non-empty string")
    if not _is_non_empty_str(benchling_id):
        raise ValueError("benchling_id must be a non-empty string")

    payload = {
        "entity_type": entity_type,
        "entity_id": entity_id,
        "benchling_id": benchling_id,
        "details": details or {},
    }
    with _connect(db_path) as conn, _transaction(conn):
        return _append_audit_block_tx(conn, action_type="LIMS_SYNC", payload=payload)


def record_audit_event(
    action_type: str,
    payload: Dict[str, Any],
    db_path: Optional[DbPath] = None,
) -> Dict[str, Any]:
    """
    Seal an arbitrary event into the hash-linked audit chain.

    A generic appender so callers recording "what actually ran" - which planner
    produced a strategy, which engine scored a batch - do not need a bespoke
    function each. `action_type` is the block label; `payload` is sealed verbatim.
    """
    if not _is_non_empty_str(action_type):
        raise ValueError("action_type must be a non-empty string")
    if not isinstance(payload, dict):
        raise ValueError("payload must be a dict")
    with _connect(db_path) as conn, _transaction(conn):
        return _append_audit_block_tx(conn, action_type=action_type, payload=payload)


def list_audit_blocks(db_path: Optional[DbPath] = None) -> List[Dict[str, Any]]:
    """
    Return all cryptographic audit blocks in chronological order.
    """
    with _connect(db_path) as conn:
        rows = conn.execute("SELECT * FROM AuditLedgerBlocks ORDER BY block_index ASC;").fetchall()
    return [dict(r) for r in rows]


get_audit_blocks = list_audit_blocks


def verify_audit_integrity(db_path: Optional[DbPath] = None) -> Tuple[bool, int, str]:
    """
    Cryptographically verify the entire FDA 21 CFR Part 11 audit block chain.

    Checks:
        1. Monotonic block_index sequence starting at 0.
        2. Block 0 is GENESIS with prev_block_hash = '0'*64.
        3. Strict SHA-256 hash chaining: block[k].prev_block_hash == block[k-1].block_hash.
        4. Block hash authenticity: recomputed SHA-256 matches stored block_hash.
        5. HMAC-SHA256 digital signature authenticity: signature matches stored signature.

    Returns:
        (is_valid: bool, verified_count: int, message: str)
    """
    blocks = list_audit_blocks(db_path)
    if not blocks:
        return False, 0, "Audit ledger is empty (missing Genesis block)"

    for k, block in enumerate(blocks):
        idx = block["block_index"]
        if idx != k:
            return False, k, f"Block index discontinuity: expected {k}, found {idx}"

        prev_hash = block["prev_block_hash"]
        if k == 0:
            if block["action_type"] != "GENESIS":
                return False, 0, f"Block 0 action_type must be GENESIS, found {block['action_type']}"
            if prev_hash != GENESIS_PREV_HASH:
                return False, 0, f"Block 0 prev_block_hash must be {GENESIS_PREV_HASH}, found {prev_hash}"
        else:
            expected_prev = blocks[k - 1]["block_hash"]
            if prev_hash != expected_prev:
                return (
                    False,
                    k,
                    f"Hash chain broken at block {k}: prev_block_hash {prev_hash[:12]}... != parent hash {expected_prev[:12]}...",
                )

        # Verify SHA-256 block hash
        expected_hash = _compute_block_hash(
            idx,
            block["timestamp"],
            block["action_type"],
            block["payload_json"],
            prev_hash,
        )
        if block["block_hash"] != expected_hash:
            return (
                False,
                k,
                f"Data tampering detected at block {k} ({block['action_type']}): block_hash mismatch",
            )

        # Verify HMAC digital signature
        expected_sig = _compute_block_signature(block["block_hash"])
        if block["signature"] != expected_sig:
            return (
                False,
                k,
                f"Digital signature forgery detected at block {k} ({block['action_type']})",
            )

    return True, len(blocks), f"Audit chain verified intact: {len(blocks)} blocks cryptographically validated"


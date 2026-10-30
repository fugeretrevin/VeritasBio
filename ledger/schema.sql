-- =============================================================================
-- ledger/schema.sql — DOMAIN: Immutable Memory
--
-- DDL for the VeritasBio SQLite ledger (ARCHITECTURE.md §6):
--   Budget   — total allocation and cumulative spend.
--   Jobs     — every authorization decision (APPROVED / REJECTED) and its cost.
--   Evidence — ESM fitness scores returned for approved jobs. Each row records
--              the compute engine that produced the score (compute_engine) and
--              whether it is synthetic (is_mock), so a mock score can never be
--              mistaken for a real GPU inference result in the audit trail.
--
-- Columns match §6 exactly. CHECK constraints and triggers are defense-in-depth
-- so the database itself refuses states the Gatekeeper should never produce
-- (overspend, negative cost, rewritten history). The Python API in
-- ledger/database.py remains the only sanctioned write path.
--
-- MONEY UNITS: all money columns (total_allocation, cumulative_spent,
-- estimated_cost) are INTEGER micro-USD (1 USD = 1,000,000) for exact
-- arithmetic. Column names are kept identical to §6; ledger/database.py
-- converts to/from float USD at the API boundary. (Decision: project_context.md)
--
-- STATUS: §6 lists PENDING/APPROVED/REJECTED/COMPLETED. FAILED is added for
-- jobs that were approved and charged but whose Modal execution failed; the
-- funds stay spent. (Decision: project_context.md)
--
-- Idempotent: safe to run on every startup.
-- =============================================================================

PRAGMA foreign_keys = ON;

-- -----------------------------------------------------------------------------
-- Budget
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS Budget (
    id                INTEGER PRIMARY KEY,
    total_allocation  INTEGER NOT NULL   -- micro-USD
                      CHECK (typeof(total_allocation) = 'integer' AND total_allocation >= 0),
    cumulative_spent  INTEGER NOT NULL DEFAULT 0   -- micro-USD
                      CHECK (typeof(cumulative_spent) = 'integer' AND cumulative_spent >= 0),
    -- The ledger can never record an overspend, regardless of caller.
    CHECK (cumulative_spent <= total_allocation)
);

-- -----------------------------------------------------------------------------
-- Jobs
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS Jobs (
    job_id          TEXT    PRIMARY KEY NOT NULL
                    CHECK (typeof(job_id) = 'text' AND length(trim(job_id)) > 0),
    target_model    TEXT    NOT NULL
                    CHECK (typeof(target_model) = 'text' AND length(trim(target_model)) > 0),
    variant_count   INTEGER NOT NULL
                    CHECK (typeof(variant_count) = 'integer' AND variant_count > 0),
    estimated_cost  INTEGER NOT NULL   -- micro-USD
                    CHECK (typeof(estimated_cost) = 'integer' AND estimated_cost >= 0),
    status          TEXT    NOT NULL DEFAULT 'PENDING'
                    CHECK (status IN ('PENDING', 'APPROVED', 'REJECTED', 'COMPLETED', 'FAILED', 'STAGED_FOR_WETLAB', 'IN_FLIGHT_SYNTHESIS'))
);

-- -----------------------------------------------------------------------------
-- Evidence
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS Evidence (
    evidence_id       INTEGER PRIMARY KEY AUTOINCREMENT,
    job_id            TEXT    NOT NULL
                      REFERENCES Jobs (job_id) ON UPDATE RESTRICT ON DELETE RESTRICT,
    variant_sequence  TEXT    NOT NULL
                      CHECK (typeof(variant_sequence) = 'text' AND length(variant_sequence) > 0),
    computed_score    REAL    NOT NULL
                      CHECK (typeof(computed_score) = 'real'),
    -- Provenance of the score: which engine actually produced it. Defaults are
    -- deliberately conservative so a row whose origin was never declared is
    -- reported as unverified mock data rather than silently passing as real.
    compute_engine    TEXT    NOT NULL DEFAULT 'unknown'
                      CHECK (typeof(compute_engine) = 'text' AND length(compute_engine) > 0),
    is_mock           INTEGER NOT NULL DEFAULT 1
                      CHECK (is_mock IN (0, 1)),
    timestamp         TEXT    NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ', 'now'))
);

CREATE INDEX IF NOT EXISTS idx_evidence_job_id ON Evidence (job_id);
CREATE INDEX IF NOT EXISTS idx_jobs_status     ON Jobs (status);

-- -----------------------------------------------------------------------------
-- SpeculativeEvidence (SPEC-001 Seam 4: Kriging-Believer Asynchronous Staging)
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS SpeculativeEvidence (
    speculative_id    INTEGER PRIMARY KEY AUTOINCREMENT,
    job_id            TEXT    NOT NULL
                      REFERENCES Jobs (job_id) ON UPDATE RESTRICT ON DELETE RESTRICT,
    variant_sequence  TEXT    NOT NULL
                      CHECK (typeof(variant_sequence) = 'text' AND length(variant_sequence) > 0),
    fantasy_score     REAL    NOT NULL
                      CHECK (typeof(fantasy_score) = 'real'),
    is_provisional    INTEGER NOT NULL DEFAULT 1
                      CHECK (is_provisional IN (0, 1)),
    timestamp         TEXT    NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ', 'now'))
);

CREATE INDEX IF NOT EXISTS idx_speculative_job_id ON SpeculativeEvidence (job_id);

-- =============================================================================
-- Immutability triggers
-- =============================================================================

-- Budget: rows are permanent and spend only ever increases.
CREATE TRIGGER IF NOT EXISTS trg_budget_no_delete
BEFORE DELETE ON Budget
BEGIN
    SELECT RAISE(ABORT, 'ledger: Budget rows cannot be deleted');
END;

CREATE TRIGGER IF NOT EXISTS trg_budget_spent_monotonic
BEFORE UPDATE OF cumulative_spent ON Budget
WHEN NEW.cumulative_spent < OLD.cumulative_spent
BEGIN
    SELECT RAISE(ABORT, 'ledger: cumulative_spent cannot decrease');
END;

-- Jobs: rows are permanent; only `status` may change, along legal transitions.
CREATE TRIGGER IF NOT EXISTS trg_jobs_no_delete
BEFORE DELETE ON Jobs
BEGIN
    SELECT RAISE(ABORT, 'ledger: Jobs rows cannot be deleted');
END;

CREATE TRIGGER IF NOT EXISTS trg_jobs_no_insert_terminal
BEFORE INSERT ON Jobs
WHEN NEW.status IN ('COMPLETED', 'FAILED')
BEGIN
    SELECT RAISE(ABORT, 'ledger: a job cannot be created as COMPLETED or FAILED');
END;

CREATE TRIGGER IF NOT EXISTS trg_jobs_fields_immutable
BEFORE UPDATE OF job_id, target_model, variant_count, estimated_cost ON Jobs
BEGIN
    SELECT RAISE(ABORT, 'ledger: job identity and cost are immutable');
END;

-- Legal transitions:  PENDING -> APPROVED | REJECTED,  APPROVED -> COMPLETED | FAILED
CREATE TRIGGER IF NOT EXISTS trg_jobs_status_transition
BEFORE UPDATE OF status ON Jobs
WHEN NOT (
       (OLD.status = 'PENDING'             AND NEW.status IN ('APPROVED', 'REJECTED'))
    OR (OLD.status = 'APPROVED'            AND NEW.status IN ('COMPLETED', 'FAILED', 'STAGED_FOR_WETLAB', 'IN_FLIGHT_SYNTHESIS'))
    OR (OLD.status = 'STAGED_FOR_WETLAB'   AND NEW.status IN ('IN_FLIGHT_SYNTHESIS', 'COMPLETED', 'FAILED'))
    OR (OLD.status = 'IN_FLIGHT_SYNTHESIS' AND NEW.status IN ('COMPLETED', 'FAILED'))
)
BEGIN
    SELECT RAISE(ABORT, 'ledger: illegal job status transition');
END;

-- Evidence: append-only, and only for jobs that were actually paid for.
CREATE TRIGGER IF NOT EXISTS trg_evidence_requires_approved_job
BEFORE INSERT ON Evidence
WHEN (SELECT status FROM Jobs WHERE job_id = NEW.job_id) NOT IN ('APPROVED', 'STAGED_FOR_WETLAB', 'IN_FLIGHT_SYNTHESIS')
BEGIN
    SELECT RAISE(ABORT, 'ledger: evidence may only be recorded for APPROVED jobs');
END;

CREATE TRIGGER IF NOT EXISTS trg_evidence_no_update
BEFORE UPDATE ON Evidence
BEGIN
    SELECT RAISE(ABORT, 'ledger: Evidence rows are immutable');
END;

CREATE TRIGGER IF NOT EXISTS trg_evidence_no_delete
BEFORE DELETE ON Evidence
BEGIN
    SELECT RAISE(ABORT, 'ledger: Evidence rows cannot be deleted');
END;

-- -----------------------------------------------------------------------------
-- AuditLedgerBlocks (FDA 21 CFR Part 11 Cryptographic Merkle Hash Chain)
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS AuditLedgerBlocks (
    block_index       INTEGER PRIMARY KEY,
    timestamp         TEXT    NOT NULL,
    action_type       TEXT    NOT NULL
                      CHECK (action_type IN (
                          'GENESIS',
                          'BUDGET_CREATED',
                          'SPEND_AUTHORIZED',
                          'JOB_COMPLETED',
                          'JOB_FAILED',
                          'EVIDENCE_BATCH_RECORDED',
                          'LIMS_SYNC',
                          'SPECULATIVE_BATCH_STAGED',
                          'IN_FLIGHT_DISPATCHED',
                          'PHYSICAL_EVIDENCE_RECONCILED',
                          -- Records which planner produced each round's strategy, so a
                          -- campaign whose live LLM calls fell back to the mock planner
                          -- is distinguishable from a genuine live run.
                          'STRATEGY_PLANNED'
                      )),
    payload_json      TEXT    NOT NULL,
    prev_block_hash   TEXT    NOT NULL,
    block_hash        TEXT    NOT NULL UNIQUE,
    signature         TEXT    NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_audit_prev_hash ON AuditLedgerBlocks (prev_block_hash);

-- Audit blocks are strictly immutable: no updates and no deletes
CREATE TRIGGER IF NOT EXISTS trg_audit_no_update
BEFORE UPDATE ON AuditLedgerBlocks
BEGIN
    SELECT RAISE(ABORT, 'ledger: AuditLedgerBlocks are cryptographically immutable');
END;

CREATE TRIGGER IF NOT EXISTS trg_audit_no_delete
BEFORE DELETE ON AuditLedgerBlocks
BEGIN
    SELECT RAISE(ABORT, 'ledger: AuditLedgerBlocks cannot be deleted');
END;


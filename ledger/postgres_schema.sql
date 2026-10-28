-- =============================================================================
-- VeritasBio Enterprise Multi-Tenant PostgreSQL Schema with Row-Level Security
-- Standard: FDA 21 CFR Part 11 & SOC-2 Type II Regulatory Audit Compliance
-- =============================================================================

CREATE EXTENSION IF NOT EXISTS "pgcrypto";

-- -----------------------------------------------------------------------------
-- 1. Organizations & Tenancy
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS organizations (
    org_id          UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    slug            TEXT UNIQUE NOT NULL,
    name            TEXT NOT NULL,
    tier            TEXT NOT NULL DEFAULT 'ENTERPRISE' CHECK (tier IN ('FREE', 'PRO', 'ENTERPRISE')),
    created_at      TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

CREATE TABLE IF NOT EXISTS users (
    user_id         UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    org_id          UUID NOT NULL REFERENCES organizations(org_id) ON DELETE CASCADE,
    email           TEXT UNIQUE NOT NULL,
    full_name       TEXT NOT NULL,
    role            TEXT NOT NULL CHECK (role IN ('ADMIN', 'SCIENTIST', 'AUDITOR', 'VIEWER')),
    created_at      TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

-- -----------------------------------------------------------------------------
-- 2. Multi-Tenant Project Budgets (Micro-USD integer precision)
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS budgets (
    budget_id                 UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    org_id                    UUID NOT NULL REFERENCES organizations(org_id) ON DELETE CASCADE,
    project_slug              TEXT NOT NULL,
    total_allocation_micro    BIGINT NOT NULL CHECK (total_allocation_micro > 0),
    cumulative_spent_micro    BIGINT NOT NULL DEFAULT 0 CHECK (cumulative_spent_micro >= 0),
    created_at                TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    CONSTRAINT uq_org_project UNIQUE (org_id, project_slug),
    CONSTRAINT chk_spent_lte_alloc CHECK (cumulative_spent_micro <= total_allocation_micro)
);

-- -----------------------------------------------------------------------------
-- 3. Compute Jobs Ledger
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS jobs (
    job_id                    TEXT PRIMARY KEY,
    org_id                    UUID NOT NULL REFERENCES organizations(org_id) ON DELETE CASCADE,
    budget_id                 UUID NOT NULL REFERENCES budgets(budget_id) ON DELETE RESTRICT,
    target_model              TEXT NOT NULL,
    variant_count             INTEGER NOT NULL CHECK (variant_count > 0),
    estimated_cost_micro      BIGINT NOT NULL CHECK (estimated_cost_micro > 0),
    status                    TEXT NOT NULL CHECK (status IN ('PENDING', 'APPROVED', 'REJECTED', 'COMPLETED', 'FAILED')),
    created_at                TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_jobs_org_status ON jobs (org_id, status);

-- -----------------------------------------------------------------------------
-- 4. Evidence (Append-only assayed / evaluated experimental variants)
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS evidence (
    evidence_id               BIGSERIAL PRIMARY KEY,
    org_id                    UUID NOT NULL REFERENCES organizations(org_id) ON DELETE CASCADE,
    job_id                    TEXT NOT NULL REFERENCES jobs(job_id) ON DELETE RESTRICT,
    variant_sequence          TEXT NOT NULL,
    computed_score            DOUBLE PRECISION NOT NULL,
    created_at                TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_evidence_org_job ON evidence (org_id, job_id);

-- -----------------------------------------------------------------------------
-- 5. FDA 21 CFR Part 11 Audit Ledger Blocks (Cryptographic Hash-Chained)
-- -----------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS audit_ledger_blocks (
    block_index               BIGINT NOT NULL,
    org_id                    UUID NOT NULL REFERENCES organizations(org_id) ON DELETE CASCADE,
    timestamp                 TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    action_type               TEXT NOT NULL CHECK (action_type IN (
                                  'GENESIS',
                                  'BUDGET_CREATED',
                                  'SPEND_AUTHORIZED',
                                  'JOB_COMPLETED',
                                  'JOB_FAILED',
                                  'EVIDENCE_BATCH_RECORDED',
                                  'LIMS_SYNC'
                              )),
    payload_json              JSONB NOT NULL,
    prev_block_hash           TEXT NOT NULL,
    block_hash                TEXT NOT NULL,
    signature                 TEXT NOT NULL,
    PRIMARY KEY (org_id, block_index),
    CONSTRAINT uq_audit_block_hash UNIQUE (org_id, block_hash)
);

CREATE INDEX IF NOT EXISTS idx_audit_prev_hash ON audit_ledger_blocks (org_id, prev_block_hash);

-- -----------------------------------------------------------------------------
-- 6. Row-Level Security (RLS) Policies
-- -----------------------------------------------------------------------------
ALTER TABLE budgets ENABLE ROW LEVEL SECURITY;
ALTER TABLE jobs ENABLE ROW LEVEL SECURITY;
ALTER TABLE evidence ENABLE ROW LEVEL SECURITY;
ALTER TABLE audit_ledger_blocks ENABLE ROW LEVEL SECURITY;

CREATE POLICY tenant_isolation_budgets ON budgets
    FOR ALL
    USING (org_id = NULLIF(current_setting('app.current_org_id', true), '')::UUID);

CREATE POLICY tenant_isolation_jobs ON jobs
    FOR ALL
    USING (org_id = NULLIF(current_setting('app.current_org_id', true), '')::UUID);

CREATE POLICY tenant_isolation_evidence ON evidence
    FOR ALL
    USING (org_id = NULLIF(current_setting('app.current_org_id', true), '')::UUID);

CREATE POLICY tenant_isolation_audit ON audit_ledger_blocks
    FOR ALL
    USING (org_id = NULLIF(current_setting('app.current_org_id', true), '')::UUID);

-- -----------------------------------------------------------------------------
-- 7. Immutability Triggers (No UPDATES or DELETES on Evidence or Audit Blocks)
-- -----------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION trg_abort_immutable_operation()
RETURNS TRIGGER AS $$
BEGIN
    RAISE EXCEPTION 'ledger: this table is cryptographically immutable and cannot be modified or deleted';
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER trg_no_update_evidence
    BEFORE UPDATE OR DELETE ON evidence
    FOR EACH ROW EXECUTE FUNCTION trg_abort_immutable_operation();

CREATE TRIGGER trg_no_update_audit
    BEFORE UPDATE OR DELETE ON audit_ledger_blocks
    FOR EACH ROW EXECUTE FUNCTION trg_abort_immutable_operation();

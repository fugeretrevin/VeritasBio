"""
benchmark_end_to_end.py - Full-stack production campaign against the same baselines.

What this adds over benchmark_panoptes_isolated.py:
    That file deliberately exercises only steps 1-3 (and 6b), because step 5 was
    hard-wired to ESM inference and there was no way to substitute a measured assay as
    the outcome source. It therefore cannot say anything about the Gatekeeper, the
    ledger, Modal GPU execution or the audit chain.

    This file drives the real `Orchestrator.run_stream()` loop end to end:

      step 1  fetch_protein_metadata                      (real UniProt/PDB metadata)
      step 2  Orchestrator._generate_strategy              (LIVE Muse / Llama planner)
      step 3  select_next_micro_batch                      (DKL GP + GFlowNet + Pareto + k-DPP)
      step 3b spatial filter / relaxation                  (as configured in production)
      step 4  Gatekeeper.evaluate_and_authorize            (REAL atomic spend against SQLite)
      step 5  Modal A10G ESM-2 650M  +  measured assay     (REAL GPU call, measured outcome)
      step 6  ledger evidence + HMAC-SHA256 Merkle chain   (verified after every campaign)
      step 6b Orchestrator._generate_reflection            (LIVE residual reflection)

    Nothing is mocked. `--allow-mock-fallback` is off by default, so a Modal failure or
    an LLM failure aborts the arm loudly instead of silently degrading to a mock and
    reporting the result as if it were live.

Why step 5 does both ESM and the assay:
    Running the loop on ESM pseudo-likelihoods alone would optimise ESM rather than
    fitness and would consume no experimental budget - the exact reason the isolated
    benchmark refuses to do it. Running it on the assay alone would never touch Modal.
    A real DBTL round does both: ESM scores the batch on GPU to produce the in-silico
    record, then the same batch is synthesised and assayed, and the GP learns from the
    MEASURED numbers. So both run, both are charged, both are recorded, and the ESM
    score is kept alongside the measured value so their correlation is reportable
    rather than assumed.

No hardcoded library:
    The saturated positions and the parent residues are derived from the dataset's own
    variant notations (`derive_library_scope`). Nothing in this file names a position
    or a residue of any particular landscape, so pointing it at another combinatorial
    dataset requires no code change. The library itself is resolved once from
    `--library` and threaded explicitly through every call path; the `ACTIVE_LIBRARY`,
    `ACTIVE_SITE_POSITIONS` and `PARENT_RESIDUES` globals this file used to read are
    gone, and the derived scope is cross-checked against the library's DECLARED
    positions so a dataset and a library specification cannot silently disagree.

Manifest and provenance (migrated to match benchmark_panoptes_isolated):
    The manifest is serialised through that harness's `_json_safe`, NOT through
    `json.dumps(default=str)`. A manifest that needs `default=str` is not
    machine-readable: tuples, sets, Paths and enums all stringify into shapes no
    consumer can parse back, and a value that cannot be serialised becomes an opaque
    string instead of a loud failure. It is then validated against
    `MANIFEST_REQUIRED_PATHS` before it is written.

    Label provenance uses the harness's `LabelProvenance` enum, which is deliberately
    not a `str` subclass, so the old `provenance_status == "UNVERIFIED_INCLUDES_IMPUTED"`
    comparison cannot be written by accident. The legacy string is KEPT beside the new
    stamp - in the CSV stamp, the report and the manifest - so existing run folders and
    the documented greps keep matching.

Comparison:
    The baselines (random, ridge, random_forest, alde) and the scoring harness are
    imported from benchmark_panoptes_isolated rather than reimplemented, so the
    competitor arms are byte-for-byte the same methods under the same protocol, the
    same shared starting sample per seed, and the same evaluator.

Usage:
    python benchmark_end_to_end.py --help
    python benchmark_end_to_end.py --dry-run          # no GPU, no LLM, no spend
    python benchmark_end_to_end.py --provenance-mode acknowledge-unflagged --seeds 3
"""

from __future__ import annotations

import sys

sys.dont_write_bytecode = True

import argparse
import concurrent.futures
import json
import logging
import os
import random
import subprocess
import tempfile
import time
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

# The isolated benchmark is imported read-only, for its dataset loading, evaluator,
# protocol harness and competitor selectors. Importing it does not run it: its CLI is
# guarded under __main__.
import benchmark_panoptes_isolated as iso
from veritas import paths

# Resolved rather than assumed: `paths` answers with the REAL repository root even when
# this file is executed from a git worktree, where the gitignored output roots do not
# exist. `Path(__file__).parent` would silently point a run's artifacts into a worktree
# that nobody then looks in.
REPO_ROOT = paths.repository_root()
DEFAULT_OUTPUT_ROOT = paths.results_dir("benchmark_e2e_results", must_exist=False)
BENCHMARK_VERSION = "1.1.0"

logger = logging.getLogger("panoptes_e2e_benchmark")


class E2EError(RuntimeError):
    """Actionable configuration or execution failure in the end-to-end benchmark."""


# ==============================================================================================
# Label provenance: one vocabulary, taken from the harness rather than restated
# ==============================================================================================
# `LabelProvenance` is read from benchmark_panoptes_isolated, which this file already
# imports for its dataset loading, evaluator and protocol. One definition, one mapping,
# one stamp format. It is deliberately NOT a `str` subclass, so the comparison this file
# used to make - `provenance_status == "UNVERIFIED_INCLUDES_IMPUTED"` - is now a type
# error waiting to happen rather than a quietly false branch.
LabelProvenance = iso.LabelProvenance


def resolved_provenance(universe: "iso.Universe") -> "iso.LabelProvenance":
    """The run's label provenance as an enum member, never as a bare string."""
    return iso.map_legacy_provenance_status(universe.provenance_status)


def provenance_stamp(universe: "iso.Universe") -> str:
    """
    The one-line stamp written into every CSV and the report.

    Identical in format to the isolated harness's: the explicit vocabulary leads and the
    legacy string follows, so documented greps for `UNVERIFIED_INCLUDES_IMPUTED` keep
    matching run folders produced by either file.
    """
    return iso.provenance_stamp(universe)


def label_provenance_block(
    universe: "iso.Universe",
    args: argparse.Namespace,
) -> Dict[str, Any]:
    """
    The manifest's provenance claim, built by the harness so both files emit one shape.

    The only thing added here is what the CLI asked for, which the harness leaves for
    its caller to fill, and the arm's own name - so a reader of this manifest can see
    that an end-to-end run and an isolated run made the same claim about the same data.
    """
    block = dict(iso._label_provenance_block(universe))
    block["provenance_mode_requested"] = args.provenance_mode
    block["benchmark"] = "end_to_end"
    return block


# ==============================================================================================
# Library scope, derived from the dataset rather than declared
# ==============================================================================================
@dataclass(frozen=True)
class LibraryScope:
    """
    What the dataset's library actually varies, recovered from its variant notations.

    `positions` are the saturated sites, ascending. `parent_residues` maps each of
    those positions to the wild-type residue the notations agree on. Both are derived,
    so no landscape-specific literal appears in this file.
    """

    positions: Tuple[int, ...]
    parent_residues: Dict[int, str]
    source: str

    def as_sorted(self) -> List[int]:
        return list(self.positions)


def derive_library_scope(variant_ids: Sequence[str], wt_sequence: str) -> LibraryScope:
    """
    Derives the saturated positions and their parent residues by parsing notations.

    Every notation is of the form `<wt><pos><mut>` joined by ':'. The union of the
    positions mentioned across the dataset is the library's scope; the wild-type letter
    is cross-checked against the supplied WT sequence so a numbering mismatch is caught
    here rather than surfacing later as a silent zero-overlap filter.
    """
    positions: Set[int] = set()
    parents: Dict[int, str] = {}
    inspected = 0

    for variant_id in variant_ids:
        if not variant_id or variant_id == "WT":
            continue
        inspected += 1
        for token in variant_id.split(":"):
            token = token.strip()
            if len(token) < 3:
                raise E2EError(f"unparseable mutation token {token!r} in variant {variant_id!r}")
            wt_aa, mut_aa, pos_text = token[0], token[-1], token[1:-1]
            if not pos_text.isdigit():
                raise E2EError(f"unparseable position in token {token!r} of variant {variant_id!r}")
            position = int(pos_text)
            positions.add(position)
            previous = parents.setdefault(position, wt_aa)
            if previous != wt_aa:
                raise E2EError(
                    f"dataset disagrees on the parent residue at position {position}: "
                    f"saw both {previous!r} and {wt_aa!r}"
                )

    if not positions:
        raise E2EError(
            "no mutated positions found in the dataset's variant ids; cannot derive the "
            "library scope, so there is nothing to benchmark"
        )

    for position, wt_aa in sorted(parents.items()):
        if position > len(wt_sequence):
            raise E2EError(
                f"position {position} exceeds the {len(wt_sequence)}-residue parent sequence"
            )
        actual = wt_sequence[position - 1]
        if actual != wt_aa:
            raise E2EError(
                f"numbering mismatch at position {position}: the dataset's notations say the "
                f"parent residue is {wt_aa!r} but the parent sequence has {actual!r}. The "
                "dataset and the metadata sequence are not on the same numbering."
            )

    return LibraryScope(
        positions=tuple(sorted(positions)),
        parent_residues=dict(sorted(parents.items())),
        source=f"derived from {inspected:,} dataset variant notations, cross-checked against WT",
    )


def check_scope_against_library(scope: LibraryScope, library: Any) -> Dict[str, Any]:
    """
    Cross-checks the DERIVED scope against the library's DECLARED design.

    Deriving the scope from the data is what keeps this file free of landscape literals,
    but a derived scope that mutates a position the library does not declare means the
    dataset and the library specification are not describing the same experiment, and
    every downstream filter would then be computed over the wrong sites. That is an
    error, not a note.

    A derived scope that is a strict SUBSET is reported rather than refused: a sparse
    pool-restricted library can legitimately leave a declared position wild type in
    every record it happens to contain.
    """
    declared = tuple(int(p) for p in library.positions)
    derived = tuple(scope.positions)
    unexpected = sorted(set(derived) - set(declared))
    if unexpected:
        raise E2EError(
            f"the dataset mutates position(s) {unexpected} that library "
            f"{library.key!r} does not declare (declared: {list(declared)}). The dataset "
            "and the library specification are not describing the same design; nothing "
            "downstream can be trusted to filter on the right sites."
        )
    absent = sorted(set(declared) - set(derived))
    if absent:
        logger.warning(
            "library %s declares position(s) %s that never vary in the eligible universe",
            library.key, absent,
        )
    return {
        "library": library.key,
        "declared_positions": list(declared),
        "derived_positions": list(derived),
        "declared_but_never_varying": absent,
        "derived_source": scope.source,
    }


def notation_from_sequence(wt_sequence: str, sequence: str) -> str:
    """Builds canonical ascending notation by diffing a full sequence against wild type."""
    if len(sequence) != len(wt_sequence):
        raise E2EError(
            f"sequence length {len(sequence)} does not match the {len(wt_sequence)}-residue parent"
        )
    parts = [
        f"{wt_aa}{index}{mut_aa}"
        for index, (wt_aa, mut_aa) in enumerate(zip(wt_sequence, sequence), start=1)
        if wt_aa != mut_aa
    ]
    return ":".join(parts)


# ==============================================================================================
# Step 5: real Modal ESM-2 scoring plus the measured assay outcome
# ==============================================================================================
@dataclass
class _RoundRecord:
    """What one authorised round actually did."""

    round_number: int
    variant_ids: List[str]
    measured: List[float]
    esm_scores: List[float]
    esm_engine: str
    esm_seconds: float


class ModalEsmPlusAssaySource:
    """
    Step 5 outcome source: real GPU inference AND the real measured outcome.

    `score()` does, in order, exactly what a DBTL round does:
      1. Sends the authorised batch to Modal A10G for ESM-2 masked-marginal scoring.
         `use_mock_fallback=False`, so an unavailable Modal aborts the arm rather than
         quietly substituting synthetic numbers.
      2. Looks up each variant's MEASURED fitness and returns those as the learning
         signal, so the GP is fitted on experimental data rather than on ESM.

    The ESM scores are retained per round so the report can state their rank
    correlation with the measured values instead of assuming ESM is informative.

    The truth table is supplied through the evaluator's documented
    `export_truth_table()` hook - the same mechanism ALDE uses - because this arm runs
    its own internal loop. It reads an entry only after the Gatekeeper has authorised
    that batch, preserving commit-then-reveal order.
    """

    engine_name = "modal_esm2_650m+measured_assay"
    target_model = "esm_triage_plus_wetlab_assay"

    def __init__(
        self,
        truth: Dict[str, float],
        wt_sequence: str,
        esm_model: str,
        allow_mock_fallback: bool,
        run_esm: bool = True,
    ) -> None:
        self._truth = truth
        self._wt = wt_sequence
        self._esm_model = esm_model
        self._allow_mock_fallback = allow_mock_fallback
        self._run_esm = run_esm
        self.rounds: List[_RoundRecord] = []
        self.unknown_variants: List[str] = []

    def score(
        self,
        wt_sequence: str,
        variant_sequences: List[str],
        round_number: int,
    ) -> List[float]:
        variant_ids = [notation_from_sequence(self._wt, seq) for seq in variant_sequences]

        # Fail loudly on anything the assay cannot measure. With the library scope
        # enforced upstream this should be unreachable; if it fires, the scope is wrong
        # and fabricating a value would corrupt the comparison.
        missing = [vid for vid in variant_ids if vid not in self._truth]
        if missing:
            self.unknown_variants.extend(missing)
            raise E2EError(
                f"round {round_number}: {len(missing)} authorised variant(s) are absent from the "
                f"measured dataset, e.g. {missing[:3]}. The enforced library scope and the "
                "dataset's eligible universe disagree."
            )

        esm_scores: List[float] = []
        esm_engine = "not_run"
        esm_seconds = 0.0
        if self._run_esm:
            from workers.esm_worker import score_variants_with_provenance

            started = time.perf_counter()
            outcome = score_variants_with_provenance(
                wt_sequence=wt_sequence,
                variant_sequences=variant_sequences,
                target_model=self._esm_model,
                use_mock_fallback=self._allow_mock_fallback,
            )
            esm_seconds = time.perf_counter() - started
            esm_scores = [float(s) for s in outcome.scores]
            esm_engine = str(outcome.engine)
            if getattr(outcome, "is_mock", False) and not self._allow_mock_fallback:
                raise E2EError(
                    f"round {round_number}: Modal returned mock scores while "
                    "--allow-mock-fallback was not set"
                )

        measured = [float(self._truth[vid]) for vid in variant_ids]
        self.rounds.append(_RoundRecord(
            round_number=round_number,
            variant_ids=variant_ids,
            measured=measured,
            esm_scores=esm_scores,
            esm_engine=esm_engine,
            esm_seconds=round(esm_seconds, 3),
        ))
        return measured


# ==============================================================================================
# The production arm, plugged into the isolated benchmark's protocol harness
# ==============================================================================================
class ProductionCampaignSelector(iso.Selector):
    """
    Drives the real `run_stream()` loop for one seed and reports the batches it chose.

    Implements the harness's `run_campaign` / `attach_truth` contract, which exists for
    methods that own their internal loop. The harness then re-reveals and re-scores
    every returned pick with its own evaluator, so the metrics are computed identically
    to every other arm.
    """

    def __init__(
        self,
        seed: int,
        scope: LibraryScope,
        metadata: Any,
        wt_sequence: str,
        pdb_path: Optional[Path],
        args: argparse.Namespace,
    ) -> None:
        super().__init__(seed)
        self.name = "panoptes_end_to_end"
        self.scope = scope
        self.metadata = metadata
        self.wt_sequence = wt_sequence
        self.pdb_path = pdb_path
        self.args = args
        self._truth: Dict[str, float] = {}
        self.source: Optional[ModalEsmPlusAssaySource] = None
        self.ledger_summary: Dict[str, Any] = {}
        self.token_usage: Dict[str, int] = {}
        self.config = {
            "planner": "live Muse" if args.enable_live_muse else "annotation-derived (no live LLM)",
            "outcome": ModalEsmPlusAssaySource.engine_name if not args.dry_run else "measured assay only",
            "priced_as": ModalEsmPlusAssaySource.target_model,
            "assayable_positions": scope.as_sorted(),
            "assayable_positions_source": scope.source,
            "steps_exercised": ["1_init", "2_strategy", "3_selection", "4_authorization",
                                "5_execution", "6_evidence", "6b_reflection"],
            "steps_skipped": [],
        }

    def attach_truth(self, truth: Dict[str, float]) -> None:
        self._truth = dict(truth)

    def run_campaign(
        self,
        view: iso.CandidateView,
        initial_ids: Sequence[str],
        n_rounds: int,
        batch_size: int,
    ) -> List[List[str]]:
        from controller.gatekeeper import Gatekeeper
        from controller.orchestrator import Orchestrator, OrchestratorConfig
        from ledger import database as ledger_db

        if not self._truth:
            raise E2EError("the evaluator did not supply a truth table to the production arm")

        with tempfile.TemporaryDirectory(prefix="panoptes_e2e_") as tmp_dir:
            db_path = str(Path(tmp_dir) / f"ledger_seed{self.seed}.db")
            ledger_db.init_db(db_path)

            # Size the allocation from the real cost model rather than a guess, with a
            # margin for the cold-start pass, so the campaign is not halted by an
            # under-provisioned budget that says nothing about selection quality.
            gatekeeper = Gatekeeper(db_path=Path(db_path))
            # Round 1 may run a larger cold-start batch than `batch_size`, and the
            # structural filter can widen a round, so price every round at the largest
            # batch the orchestrator might submit. Under-provisioning here would halt
            # the campaign on budget and say nothing about selection quality.
            widest_batch = max(batch_size, 2 * batch_size)
            per_round, _warm = gatekeeper.calculate_cost(
                variant_count=widest_batch,
                target_model=ModalEsmPlusAssaySource.target_model,
                sequence_length_aa=len(self.wt_sequence),
                force_warm=False,
            )
            # The shared starting plate is a real measured batch, so it is priced and
            # authorised like any other. Budget covers it plus every planned round.
            initial_cost, _ = gatekeeper.calculate_cost(
                variant_count=len(initial_ids),
                target_model="wetlab_assay",
                sequence_length_aa=len(self.wt_sequence),
                force_warm=False,
            )
            # The Tier-3 sweep is a real authorised batch, so the allocation has to cover
            # it or the Gatekeeper fail-closes and the stage silently never runs.
            refine_cost = 0.0
            if self.args.enable_lead_refinement:
                refine_variants = max(40, batch_size * 2 * self.args.lead_refinement_top_n)
                refine_cost, _ = gatekeeper.calculate_cost(
                    variant_count=refine_variants,
                    target_model=ModalEsmPlusAssaySource.target_model,
                    sequence_length_aa=len(self.wt_sequence),
                    force_warm=False,
                )
            budget = round((per_round * n_rounds + initial_cost + refine_cost) * 1.15, 2)
            ledger_db.create_budget(total_allocation=budget, db_path=db_path)

            # Seed the shared starting sample as real evidence, so this arm plans its
            # first round from the same measurements the baselines fit on. Without it
            # the production arm would choose round 1 blind while ridge and random
            # forest train on those 25 points.
            from controller.gatekeeper import Variant as GkVariant

            initial_pairs = [
                (GkVariant(notation=vid).apply_to(self.wt_sequence), float(self._truth[vid]))
                for vid in initial_ids
            ]
            seed_job = f"e2e_seed{self.seed}_initial_{uuid.uuid4().hex[:8]}"
            approved = ledger_db.authorize_spend(
                job_id=seed_job,
                target_model="wetlab_assay",
                variant_count=len(initial_pairs),
                estimated_cost=float(initial_cost),
                db_path=db_path,
            )
            if not approved:
                raise E2EError(
                    f"the ledger refused to authorise the shared starting plate "
                    f"(${initial_cost:,.2f} against a ${budget:,.2f} allocation)"
                )
            ledger_db.record_evidence(
                seed_job, initial_pairs, db_path=db_path,
                compute_engine="shared_initial_sample",
            )
            ledger_db.mark_job_completed(seed_job, db_path=db_path)

            self.source = ModalEsmPlusAssaySource(
                truth=self._truth,
                wt_sequence=self.wt_sequence,
                esm_model=self.args.esm_model,
                allow_mock_fallback=self.args.allow_mock_fallback,
                run_esm=not self.args.dry_run,
            )

            config = OrchestratorConfig(
                target_accession=self.args.accession,
                total_budget_usd=budget,
                max_iterations=n_rounds,
                batch_size=batch_size,
                # Pin round 1 to the same size as every other round. Left unset this
                # expands to 2x batch_size to amortise container spin-up, which is
                # sensible in production and ruinous for a benchmark: the orchestrator
                # would fit its GP on 50 round-1 measurements while every baseline
                # learns from 25, an information advantage the scored measurement
                # count does not reveal.
                cold_start_batch_size=batch_size,
                db_path=db_path,
                # Live planner unless explicitly disabled; never silently mocked.
                use_mock_llm=not self.args.enable_live_muse,
                enable_live_llm=bool(self.args.enable_live_muse),
                require_live_llm=bool(self.args.enable_live_muse) and not self.args.allow_mock_fallback,
                enable_reflection=bool(self.args.enable_live_muse),
                allow_mock_scoring_fallback=bool(self.args.allow_mock_fallback),
                # Run the configured number of rounds so every arm gets the same
                # measurement budget; early stopping is reported separately when on.
                enable_regret_stopping=bool(self.args.enable_stopping),
                enable_stagnation_stopping=bool(self.args.enable_stopping),
                # Tier-3 Hamming-1 polish around the top leads: the engine's only
                # exploitation stage, and the obvious lever for the flat peak-fitness
                # result. Deliberately OPT-IN rather than on by default, because it
                # measures up to `max(40, batch_size*2*top_n)` ADDITIONAL variants that
                # the benchmark harness does not count but the GP does learn from. On by
                # default it would silently hand this arm a larger measurement budget
                # than the baselines - the same unfairness as the cold-start batch.
                # Enable it only alongside a reduced --rounds so totals stay matched.
                enable_lead_refinement=bool(self.args.enable_lead_refinement),
                lead_refinement_top_n=self.args.lead_refinement_top_n,
                enable_structural_filter=bool(self.args.enable_structural_filter),
                structure_pdb_id=self.args.structure_pdb_id,
                enable_sidechain_relaxation=bool(self.args.enable_relaxation),
                random_seed=self.seed,
            )

            orchestrator = Orchestrator(
                config,
                outcome_source=self.source,
                # Derived from the dataset, not declared here.
                assayable_positions=self.scope.positions,
                # Never MockLLMWorker: its TrpB branch hardcodes this dataset's four
                # positions and names three of the four global-optimum residues, so
                # using it as the non-live planner would encode the answer. The
                # annotation-driven planner derives its search space from the target's
                # own annotations instead.
                planner=iso._AnnotationDrivenPlanner(),
            )

            events: List[Any] = []
            try:
                for event in orchestrator.run_stream():
                    events.append(event)
                    if event.step in ("halted", "stopped"):
                        self.notes.append(f"run_stream {event.step}: {event.message[:160]}")
            except Exception as exc:
                raise E2EError(f"run_stream failed: {type(exc).__name__}: {exc}") from exc

            self.model_calls = len(getattr(orchestrator, "strategy_engines_used", []) or [])
            engines = list(getattr(orchestrator, "strategy_engines_used", []) or [])
            if engines:
                self.notes.append("planners: " + ",".join(engines))
                # Only a FALLBACK is an error. When --enable-live-muse is off the
                # injected annotation-driven planner is the intended planner for every
                # round, not a degradation, so its appearance here is expected.
                if self.args.enable_live_muse and not self.args.allow_mock_fallback:
                    fell_back = [e for e in engines if "mock" in e.lower()]
                    if fell_back:
                        raise E2EError(
                            f"the live planner fell back in {len(fell_back)} round(s) "
                            f"({fell_back}) while --allow-mock-fallback was not set"
                        )

            self.token_usage = dict(getattr(orchestrator, "token_usage", {}) or {})
            self.ledger_summary = self._summarise_ledger(db_path, budget, ledger_db)

        return [list(record.variant_ids) for record in (self.source.rounds if self.source else [])]

    def _summarise_ledger(self, db_path: str, budget: float, ledger_db: Any) -> Dict[str, Any]:
        """Records what the campaign actually spent and whether the audit chain holds."""
        remaining = float(ledger_db.get_remaining_budget(1, db_path))
        evidence = ledger_db.get_evidence(db_path=db_path)
        engines = sorted({str(row.get("compute_engine")) for row in evidence})

        # verify_audit_integrity returns (is_valid, blocks_checked, detail).
        try:
            chain_ok, blocks_checked, chain_detail = ledger_db.verify_audit_integrity(db_path=db_path)
            chain_ok = bool(chain_ok)
            chain_detail = f"{blocks_checked} block(s): {chain_detail}"
        except Exception as exc:
            chain_ok = False
            chain_detail = f"verify_audit_integrity raised {type(exc).__name__}: {exc}"

        return {
            "budget_allocated_usd": budget,
            "budget_remaining_usd": round(remaining, 6),
            "campaign_spend_usd": round(budget - remaining, 6),
            "evidence_rows": len(evidence),
            "compute_engines_recorded": engines,
            "audit_chain_valid": chain_ok,
            "audit_chain_detail": chain_detail,
        }


# ==============================================================================================
# Reporting
# ==============================================================================================
# ==============================================================================================
# Per-seed execution, serial or parallel
# ==============================================================================================
# Set in worker processes so a child can never fan out again. Parallelism is one level
# deep by construction; without this a mis-set --parallel-seeds in a child would spawn
# grandchildren and the process count would grow without bound.
CHILD_SENTINEL_ENV = "VERITAS_E2E_WORKER"

# Hard ceiling regardless of what is asked for. Protects against a typo like
# --parallel-seeds 500 launching 500 interpreters, and against tripping Meta API rate
# limits or Modal concurrency caps.
MAX_WORKERS = 4

# Each worker rebuilds the whole eligible universe: the variant id list, the combo and
# fitness dicts for ~160k entries, and the one-hot matrix.
#
# Calibrated from an observed failure rather than a guess: 8 workers on a machine with
# 3.2GB free exhausted RAM and OpenBLAS killed the processes mid-run (BrokenProcessPool,
# and a silent 0/0 result table - see the guard below). That is ~0.4GB per worker, and
# it was measured BEFORE BLAS was pinned to one thread per process, which removed the
# per-core workspace multiplier. 0.5GB leaves margin over that figure.
#
# If workers are OOM-killed again, raise this rather than lowering --parallel-seeds:
# this constant is what makes the clamp degrade to serial instead of crashing.
WORKER_MEMORY_BUDGET_GB = 0.5

# BLAS threads are counter-productive here: the seeds are already parallel, so each
# worker spawning one thread per core multiplies both memory and CPU contention for no
# gain. Set in the parent so spawned children inherit it before numpy initialises.
_THREAD_LIMIT_ENV = (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
)


def limit_blas_threads() -> None:
    """Pins BLAS to one thread per process. Must run before children are spawned."""
    for name in _THREAD_LIMIT_ENV:
        os.environ[name] = "1"


def _free_memory_gb() -> Optional[float]:
    """Best-effort free physical memory, or None when it cannot be determined."""
    try:
        if sys.platform == "win32":
            import ctypes

            class _MemStatus(ctypes.Structure):
                _fields_ = [
                    ("dwLength", ctypes.c_ulong),
                    ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
                ]

            status = _MemStatus()
            status.dwLength = ctypes.sizeof(_MemStatus)
            if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
                return None
            return status.ullAvailPhys / (1024 ** 3)
        return (os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")) / (1024 ** 3)
    except Exception:
        return None


def resolve_worker_count(args: argparse.Namespace, n_seeds: int) -> int:
    """
    Bounded worker count: never more than seeds, never more than MAX_WORKERS, never
    more than free memory can hold, and never inside a child process.
    """
    requested = int(getattr(args, "parallel_seeds", 1) or 1)
    if requested <= 1:
        return 1
    if os.getenv(CHILD_SENTINEL_ENV):
        logger.warning(
            "ignoring --parallel-seeds %d: already running inside a worker process",
            requested,
        )
        return 1

    workers = max(1, min(requested, n_seeds, MAX_WORKERS))

    free_gb = _free_memory_gb()
    if free_gb is not None:
        # Leave a gigabyte for the parent and the OS rather than filling RAM exactly.
        affordable = max(1, int((free_gb - 1.0) / WORKER_MEMORY_BUDGET_GB))
        if affordable < workers:
            logger.warning(
                "reducing workers %d -> %d: only %.1fGB free and each worker needs "
                "about %.1fGB for the eligible universe",
                workers, affordable, free_gb, WORKER_MEMORY_BUDGET_GB,
            )
            workers = affordable

    if workers < requested:
        logger.warning("clamping --parallel-seeds %d to %d (seeds=%d, max=%d)",
                       requested, workers, n_seeds, MAX_WORKERS)
    return workers


@dataclass
class SeedContext:
    """Everything one seed needs. Rebuilt inside each worker rather than pickled."""

    universe: Any
    evaluator: Any
    view: Any
    scope: LibraryScope
    metadata: Any
    wt_sequence: str
    pdb_path: Optional[Path]


def build_seed_context(args: argparse.Namespace) -> SeedContext:
    """Loads the dataset and derives the library scope. Deterministic for fixed args."""
    library = iso.resolve_library(args.library)
    universe = iso.load_universe(
        library=library,
        dataset_path=paths.explicit_path(args.dataset, what="--dataset", must_be_dir=False),
        provenance_mode=args.provenance_mode,
        imputed_ids_path=(
            paths.explicit_path(args.imputed_ids, what="--imputed-ids", must_be_dir=False)
            if args.imputed_ids else None
        ),
        max_variants=args.max_variants,
        subsample_seed=args.subsample_seed,
    )
    metadata, wt_sequence, pdb_path, _design = iso._load_panoptes_context(library, pdb=True)
    scope = derive_library_scope(universe.variant_ids, wt_sequence)
    check_scope_against_library(scope, library)
    return SeedContext(
        universe=universe,
        evaluator=iso.Evaluator(universe),
        view=iso.CandidateView(
            variant_ids=list(universe.variant_ids),
            one_hot=universe.one_hot,
            index_of=dict(universe.index_of),
            combos=dict(universe.combos),
            library=library,
            batch_size=args.batch_size,
            n_rounds=args.rounds,
        ),
        scope=scope,
        metadata=metadata,
        wt_sequence=wt_sequence,
        pdb_path=pdb_path,
    )


def run_seed(
    ctx: SeedContext,
    seed: int,
    methods: Sequence[str],
    status: Dict[str, Dict[str, str]],
    args: argparse.Namespace,
) -> Tuple[List[Any], List[Dict[str, Any]]]:
    """
    Runs every available method for one seed against the shared starting sample.

    Returns (outcomes, production_facts). Production facts are returned as plain dicts
    rather than selector objects so the result crosses a process boundary cleanly.
    """
    # Identical draw to benchmark_panoptes_isolated.execute().
    sampler = random.Random(seed)
    initial_ids = sampler.sample(ctx.universe.variant_ids, args.initial_size)

    outcomes: List[Any] = []
    production: List[Dict[str, Any]] = []

    for method in methods:
        if status.get(method, {}).get("status") != "available":
            continue
        logger.info("seed %d | %s", seed, method)
        if method == "random":
            selector: Any = iso.RandomSelector(seed)
        elif method == "ridge":
            selector = iso.RidgeSelector(seed)
        elif method == "random_forest":
            selector = iso.RandomForestSelector(seed)
        else:
            selector = ProductionCampaignSelector(
                seed=seed, scope=ctx.scope, metadata=ctx.metadata,
                wt_sequence=ctx.wt_sequence, pdb_path=ctx.pdb_path, args=args,
            )
        outcome = iso.run_single(selector, ctx.evaluator, ctx.view, initial_ids, seed)
        outcomes.append(outcome)
        if isinstance(selector, ProductionCampaignSelector):
            production.append({
                "seed": seed,
                "ledger": dict(selector.ledger_summary),
                "token_usage": dict(selector.token_usage or {}),
                "rounds": [
                    {
                        "round_number": r.round_number,
                        "measured": r.measured,
                        "esm_scores": r.esm_scores,
                        "esm_engine": r.esm_engine,
                        "esm_seconds": r.esm_seconds,
                    }
                    for r in (selector.source.rounds if selector.source else [])
                ],
            })
        logger.info(
            "  -> best %.6f | top1 %d | %d measurements | status %s",
            outcome.rounds[-1].best_fitness if outcome.rounds else float("nan"),
            outcome.rounds[-1].cum_top1_including_initial if outcome.rounds else -1,
            outcome.rounds[-1].n_measurements if outcome.rounds else -1,
            outcome.status,
        )
    return outcomes, production


def _run_seed_worker(
    seed: int,
    methods: Sequence[str],
    status: Dict[str, Dict[str, str]],
    args: argparse.Namespace,
) -> Tuple[List[Any], List[Dict[str, Any]]]:
    """
    Worker entry point: one process per seed.

    Processes rather than threads deliberately. `select_next_micro_batch` draws its
    Thompson sample from numpy's GLOBAL RNG, which the harness seeds per call; threads
    would race on that shared state and the run would stop being reproducible. A
    separate process has its own RNG, its own ledger file and its own Modal client.

    The context is rebuilt here instead of inherited, which keeps a ~51MB one-hot
    matrix off the pickle path.
    """
    os.environ[CHILD_SENTINEL_ENV] = "1"
    logging.basicConfig(
        level=logging.INFO,
        format=f"%(asctime)s [%(levelname)s] [seed {seed}] %(message)s",
        datefmt="%H:%M:%S",
    )
    logging.getLogger("workers").setLevel(logging.WARNING)
    logging.getLogger("controller").setLevel(logging.WARNING)
    iso._load_repo_dotenv()
    ctx = build_seed_context(args)
    return run_seed(ctx, seed, methods, status, args)


def spearman(a: Sequence[float], b: Sequence[float]) -> Optional[float]:
    """Rank correlation without a SciPy dependency; None when undefined."""
    if len(a) != len(b) or len(a) < 3:
        return None
    ra = np.argsort(np.argsort(np.asarray(a, dtype=np.float64)))
    rb = np.argsort(np.argsort(np.asarray(b, dtype=np.float64)))
    if np.std(ra) == 0 or np.std(rb) == 0:
        return None
    return float(np.corrcoef(ra, rb)[0, 1])


def summarise_production(facts: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """Aggregates what the production arm did beyond the shared selection metrics."""
    spends, chains, engines, esm_engines = [], [], set(), set()
    esm_rho: List[float] = []
    esm_seconds = 0.0
    rounds_total = 0

    tokens = {"calls": 0, "prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0}
    for entry in facts:
        for key, value in (entry.get("token_usage") or {}).items():
            if key in tokens and isinstance(value, (int, float)):
                tokens[key] += int(value)
        ledger = entry.get("ledger") or {}
        if ledger:
            spends.append(float(ledger.get("campaign_spend_usd", 0.0)))
            chains.append(ledger.get("audit_chain_valid"))
            engines.update(ledger.get("compute_engines_recorded") or [])
        for record in entry.get("rounds") or []:
            rounds_total += 1
            esm_seconds += float(record.get("esm_seconds") or 0.0)
            esm_engines.add(record.get("esm_engine"))
            if record.get("esm_scores"):
                rho = spearman(record["esm_scores"], record["measured"])
                if rho is not None:
                    esm_rho.append(rho)

    return {
        "campaigns": len(facts),
        "rounds_executed": rounds_total,
        "mean_campaign_spend_usd": round(float(np.mean(spends)), 4) if spends else None,
        "total_spend_usd": round(float(np.sum(spends)), 4) if spends else None,
        "audit_chain_valid_all": (all(c is True for c in chains) if chains else None),
        "compute_engines_recorded": sorted(engines),
        "esm_engines": sorted(esm_engines),
        "esm_gpu_seconds": round(esm_seconds, 2),
        # Provider-reported, not estimated. Zero calls means no live model ran.
        "llm_calls_with_usage": tokens["calls"],
        "llm_prompt_tokens": tokens["prompt_tokens"],
        "llm_completion_tokens": tokens["completion_tokens"],
        "llm_total_tokens": tokens["total_tokens"],
        "esm_vs_measured_spearman_mean": round(float(np.mean(esm_rho)), 4) if esm_rho else None,
        "esm_vs_measured_spearman_per_round": [round(r, 4) for r in esm_rho],
    }


def write_report(
    path: Path,
    args: argparse.Namespace,
    universe: iso.Universe,
    evaluator: iso.Evaluator,
    scope: LibraryScope,
    rows: Sequence[Dict[str, Any]],
    production: Dict[str, Any],
    status: Dict[str, Dict[str, str]],
    started: datetime,
) -> str:
    lines: List[str] = []
    add = lines.append
    rule = "-" * 78

    resolved = resolved_provenance(universe)
    add("PanoptesBio END-TO-END benchmark")
    add("=" * 78)
    add(f"Started (UTC): {started.isoformat()}")
    # Both vocabularies on the line: the explicit claim, then the legacy stamp that
    # existing greps and older run folders carry.
    add(f"Provenance:    {resolved.value} (legacy {universe.provenance_status})")
    add(f"Dataset:       {universe.fingerprint.path}")
    add(f"               sha256 {universe.fingerprint.sha256[:16]}..., "
        f"{len(universe):,} eligible variants")
    add("")
    add("Library scope (DERIVED, not declared)")
    add(rule)
    add(f"  positions:       {scope.as_sorted()}")
    add(f"  parent residues: "
        + ", ".join(f"{aa}{pos}" for pos, aa in scope.parent_residues.items()))
    add(f"  source:          {scope.source}")
    add("")
    add("Protocol")
    add(rule)
    add(f"  {args.seeds} paired seeds; {args.initial_size} shared starting variants, "
        f"then {args.rounds} rounds x {args.batch_size}.")
    add(f"  Top-1% threshold: {evaluator.top1_threshold:.6f} "
        f"({evaluator.n_top1:,} of {len(universe):,} variants qualify).")
    add(f"  Global optimum in the eligible set: {evaluator.global_best:.6g}")
    add("")
    add("Results")
    add(rule)
    add(f"  {'method':<22}{'status':<12}{'done':>6}{'measured':>10}{'best fit':>11}{'sd':>9}"
        f"{'top1':>7}{'sd':>8}{'secs':>9}")
    for row in rows:
        done = f"{row['seeds_completed']}/{row['seeds_attempted']}"
        add(f"  {row['method']:<22}{row['status']:<12}{done:>6}"
            f"{row['mean_measurements_used']:>10}{row['mean_best_fitness']:>11}"
            f"{row['std_best_fitness']:>9}{row['mean_top1_hits_including_initial']:>7}"
            f"{row['std_top1_hits_including_initial']:>8}{row['mean_runtime_seconds']:>9}")
    add("")
    add("  'best fit' is the best MEASURED fitness found within the budget, averaged over seeds.")
    add("")
    add("Production arm: what actually executed")
    add(rule)
    for key, value in production.items():
        add(f"  {key}: {value}")
    add("")
    add("Method status")
    add(rule)
    for method, info in status.items():
        add(f"  {method}: {info.get('status')}")
        detail = info.get("detail")
        if detail:
            add(f"      - {detail}")
    add("")
    add("Limitations")
    add(rule)
    for item in _limitations(args, universe, production):
        add(f"  - {item}")
    add("")

    text = "\n".join(lines)
    path.write_text(text, encoding="utf-8")
    return text


def _limitations(
    args: argparse.Namespace,
    universe: iso.Universe,
    production: Dict[str, Any],
) -> List[str]:
    items: List[str] = []
    # Compared against the enum, never against a bare string. `LabelProvenance` is not a
    # `str` subclass, so `universe.provenance_status == "UNVERIFIED_INCLUDES_IMPUTED"`
    # tested one spelling of one of six distinct claims and silently said nothing for
    # the other five. The legacy string is still quoted inside the message so a reader
    # of an old report and a reader of this one are looking at the same stamp.
    resolved = resolved_provenance(universe)
    legacy = universe.provenance_status
    if resolved is LabelProvenance.HISTORICAL_UNKNOWN_PROVENANCE:
        # The imputed count belongs to the library that was loaded, which the universe
        # carries. It used to come from a module global that only ever held TrpB's.
        imputed = universe.library.published_imputed_count
        items.append(
            f"Label provenance is {resolved.value} (legacy stamp {legacy}): up to "
            f"{imputed} published KNN-imputed values may be included in the eligible "
            "set, and they cannot be located, so no individual label in this run can be "
            "claimed as measured."
            if imputed else
            f"Label provenance is {resolved.value} (legacy stamp {legacy}) for library "
            f"{universe.library.key!r}: {universe.library.provenance_note}"
        )
    elif resolved is LabelProvenance.SIMULATED_NOISY_OUTCOME:
        items.append(
            f"Label provenance is {resolved.value} (legacy stamp {legacy}): the labels "
            "are generated, not measured. This run validates mechanism only and says "
            "nothing about biological performance."
        )
    elif resolved is LabelProvenance.HISTORICAL_MEASURED:
        items.append(
            f"Label provenance is {resolved.value} (legacy stamp {legacy}): the labels "
            "are retrospective measurements from a published landscape. They support no "
            "prospective claim; only PROSPECTIVE_PHYSICAL_MEASUREMENT does, and no run "
            "in this benchmark produces one."
        )
    else:
        items.append(
            f"Label provenance is {resolved.value} (legacy stamp {legacy}). See "
            "manifest.label_provenance for what that claim does and does not support."
        )
    if args.dry_run:
        items.append(
            "DRY RUN: no Modal GPU call and no live LLM call was made. Step 5 returned the "
            "measured assay only, so these numbers say nothing about GPU execution."
        )
    if not args.enable_live_muse:
        items.append(
            "The live planner was NOT enabled, so step 2 used the annotation-derived planner. "
            "Pass --enable-live-muse to exercise the real Muse planner."
        )
    items.append(
        "The measured outcome is a table lookup standing in for synthesis and assay. The "
        "physical steps are priced and authorised as if performed, but no DNA was ordered, "
        "so turnaround time and replicate noise are absent."
    )
    items.append(
        "The full fitness table is handed to the production arm through the evaluator's "
        "documented export hook, because run_stream() owns its loop. Every returned pick is "
        "re-scored by the evaluator, and the top-1% threshold is never exported."
    )
    if production.get("esm_vs_measured_spearman_mean") is not None:
        items.append(
            "ESM-2 scores are recorded alongside measured fitness but do NOT drive selection; "
            "the GP learns from the measured values. Their rank correlation is reported above."
        )
    if getattr(args, "enable_lead_refinement", False):
        items.append(
            "LEAD REFINEMENT WAS ENABLED: the Tier-3 Hamming-1 sweep measured additional "
            "variants beyond the per-round batches, and those measurements inform the model "
            "but are NOT counted in this run's measurement column. This arm therefore had a "
            "larger effective budget than the baselines; the comparison is not budget-matched "
            "unless --rounds was reduced to compensate."
        )
    calls = int(production.get("llm_calls_with_usage") or 0)
    if calls:
        items.append(
            f"Token usage is MEASURED, not estimated: {calls} live model call(s) reported "
            f"{production.get('llm_total_tokens')} total tokens "
            f"({production.get('llm_prompt_tokens')} prompt / "
            f"{production.get('llm_completion_tokens')} completion). Dollar cost still "
            "depends on the provider's rate card, which the client does not return."
        )
    elif args.enable_live_muse:
        items.append(
            "Live planning was enabled but no call returned a usage block, so token counts "
            "are unavailable for this run."
        )
    return items


# ==============================================================================================
# Manifest
# ==============================================================================================
# Every field a reader needs to identify this run and know what its numbers may be used
# for. Checked before the manifest is written, so an incomplete record fails the run
# instead of being discovered by whoever tries to reconstruct the arm months later.
MANIFEST_REQUIRED_PATHS: Tuple[str, ...] = (
    "run_id",
    "benchmark",
    "benchmark_version",
    "started_utc",
    "finished_utc",
    "library.key",
    "library.declared_positions",
    "library_scope.positions",
    "library_scope.parent_residues",
    "library_scope.source",
    "library_scope.declared_vs_derived",
    "config",
    "dataset",
    "label_provenance.status",
    "label_provenance.legacy_status",
    "label_provenance.definition",
    "label_provenance.limitations",
    "provenance_status",
    "methods",
    "summary",
    "production",
    "harness_exports",
    "packages",
    "repository_state",
    "environment.repository_root",
    "environment.in_worktree",
    "limitations",
)


def build_manifest(
    *,
    run_id: str,
    args: argparse.Namespace,
    library: Any,
    universe: "iso.Universe",
    scope: LibraryScope,
    scope_check: Dict[str, Any],
    evaluator: "iso.Evaluator",
    status: Dict[str, Dict[str, str]],
    rows: Sequence[Dict[str, Any]],
    production: Dict[str, Any],
    production_facts: Sequence[Dict[str, Any]],
    started: datetime,
) -> Dict[str, Any]:
    """
    The run's machine-readable record, built in one place.

    `provenance_status` is retained at the top level as the LEGACY key, unchanged, so
    the documented greps and any existing reader keep working; `label_provenance` is the
    explicit claim beside it. Keeping both is deliberate - the point of the enum is that
    one opaque stamp stopped standing in for six different claims, not that the old
    stamp disappears from the record.
    """
    return {
        "run_id": run_id,
        "benchmark": "end_to_end",
        "benchmark_version": BENCHMARK_VERSION,
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        # The library that was resolved from --library and threaded through the run,
        # recorded so the manifest names its own landscape rather than implying one.
        "library": {
            "key": library.key,
            "label": library.label,
            "accession": library.accession,
            "declared_positions": list(library.positions),
            "parent_residues": library.parent_residues,
            "selection_mode": library.selection_mode,
            "source_publication": library.source_publication,
            "source_doi": library.source_doi,
        },
        "library_scope": {
            "positions": scope.as_sorted(),
            "parent_residues": {str(k): v for k, v in scope.parent_residues.items()},
            "source": scope.source,
            "declared_vs_derived": scope_check,
        },
        "config": vars(args),
        "dataset": universe.fingerprint.to_dict(),
        # LEGACY key, retained verbatim beside the explicit vocabulary below.
        "provenance_status": universe.provenance_status,
        "label_provenance": label_provenance_block(universe, args),
        "methods": status,
        "summary": list(rows),
        "production": production,
        "harness_exports": evaluator.harness_exports,
        "packages": iso._package_versions(),
        "repository_state": iso.repository_state(),
        # Which roots this process actually resolved, so a later reader can tell whether
        # the run was executed from a worktree instead of inferring it.
        "environment": paths.describe_environment(),
        "per_campaign_ledger": [
            {"seed": e.get("seed"), **(e.get("ledger") or {})} for e in production_facts
        ],
        "limitations": _limitations(args, universe, production),
    }


def validate_manifest(manifest: Dict[str, Any]) -> None:
    """Refuses a manifest that cannot identify its own run or state its provenance."""
    missing: List[str] = []
    for path in MANIFEST_REQUIRED_PATHS:
        node: Any = manifest
        for key in path.split("."):
            if not isinstance(node, dict) or key not in node:
                missing.append(path)
                break
            node = node[key]
    if missing:
        raise iso.ManifestIncompleteError(
            "the end-to-end run manifest is missing required field(s): "
            + ", ".join(missing)
        )


# ==============================================================================================
# CLI
# ==============================================================================================
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="End-to-end production campaign benchmarked against the same baselines.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--library", choices=tuple(sorted(iso.LIBRARIES)),
                        default=iso.DEFAULT_LIBRARY_KEY,
                        help="Which registered library to run (default: %(default)s). The "
                             "library is threaded through every call path; nothing is read "
                             "from module state.")
    parser.add_argument("--dataset", default=str(iso.DEFAULT_DATASET))
    parser.add_argument("--accession",
                        default=iso.LIBRARIES[iso.DEFAULT_LIBRARY_KEY].accession)
    parser.add_argument("--methods", default="random,ridge,random_forest,panoptes_end_to_end",
                        help="Comma-separated. Default: %(default)s")
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--base-seed", type=int, default=0)
    parser.add_argument("--initial-size", type=int, default=25)
    parser.add_argument("--rounds", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=25)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--provenance-mode", choices=("require", "acknowledge-unflagged"),
                        default="require")
    parser.add_argument("--imputed-ids", default=None)
    parser.add_argument("--max-variants", type=int, default=None)
    parser.add_argument("--subsample-seed", type=int, default=0)

    parser.add_argument("--enable-live-muse", action="store_true",
                        help="Use the real Muse/Llama planner for step 2. Makes paid API calls.")
    parser.add_argument("--esm-model", default="esm2_650m",
                        help="Modal ESM model for step 5 (default: %(default)s).")
    parser.add_argument("--allow-mock-fallback", action="store_true",
                        help="Permit mock substitution when Modal or the LLM fails. OFF by "
                             "default so a degraded run fails instead of being reported as live.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Exercise steps 1-4 and 6 with no Modal GPU call and no live LLM. "
                             "Costs nothing; proves the wiring.")
    parser.add_argument("--enable-stopping", action="store_true",
                        help="Leave regret/stagnation early stopping on (production default). "
                             "Off here so every arm spends the same measurement budget.")
    parser.add_argument("--enable-lead-refinement", action="store_true",
                        help="Enable the Tier-3 Hamming-1 polish around top leads - the "
                             "engine's only exploitation stage, and the lever for peak "
                             "fitness. WARNING: it measures up to max(40, batch_size*2*top_n) "
                             "extra variants that the harness does not count, so the arm gets "
                             "a larger effective budget than the baselines. Pair it with a "
                             "reduced --rounds to keep totals matched.")
    parser.add_argument("--lead-refinement-top-n", type=int, default=3,
                        help="Number of top leads to sweep Hamming-1 neighbours around "
                             "(default: %(default)s).")
    parser.add_argument("--enable-structural-filter", action="store_true")
    parser.add_argument("--structure-pdb-id", default=None)
    parser.add_argument("--enable-relaxation", action="store_true")
    parser.add_argument("--parallel-seeds", type=int, default=1,
                        help="Run this many seeds concurrently, one process each "
                             f"(clamped to the seed count and to {MAX_WORKERS}). Seeds are "
                             "independent campaigns. Processes, not threads: "
                             "select_next_micro_batch draws from numpy's global RNG, which "
                             "threads would race on. Note this inflates the reported "
                             "per-method wall times, so keep it at 1 when timing matters.")
    parser.add_argument("--parallel-timeout", type=int, default=3600,
                        help="Hard ceiling in seconds on the whole parallel phase. On expiry "
                             "unfinished seeds are cancelled and reported missing rather than "
                             "hanging the run (default: %(default)s).")
    parser.add_argument("--verbose", action="store_true")
    return parser


def execute(args: argparse.Namespace) -> int:
    started = datetime.now(timezone.utc)
    run_id = f"e2e_{started.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    out_dir = Path(args.output_dir) / run_id
    out_dir.mkdir(parents=True, exist_ok=True)

    handler = logging.FileHandler(out_dir / "benchmark.log", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
    logging.getLogger().addHandler(handler)
    logger.info("Run %s -> %s", run_id, out_dir)

    library = iso.resolve_library(args.library)
    # An operator-supplied path is a constraint, not a hint: it resolves or the run
    # stops here, naming the flag and the value. Nothing is substituted, because a run
    # that quietly scored a different file than the one named would still report the
    # named one.
    dataset_path = paths.explicit_path(args.dataset, what="--dataset", must_be_dir=False)
    imputed_ids_path = (
        paths.explicit_path(args.imputed_ids, what="--imputed-ids", must_be_dir=False)
        if args.imputed_ids else None
    )
    universe = iso.load_universe(
        library=library,
        dataset_path=dataset_path,
        provenance_mode=args.provenance_mode,
        imputed_ids_path=imputed_ids_path,
        max_variants=args.max_variants,
        subsample_seed=args.subsample_seed,
    )
    evaluator = iso.Evaluator(universe)

    metadata, wt_sequence, pdb_path, _design = iso._load_panoptes_context(library, pdb=True)
    scope = derive_library_scope(universe.variant_ids, wt_sequence)
    # Derived from the data, then checked against what the library declares. A silent
    # disagreement here would mean every downstream site filter was computed over the
    # wrong positions.
    scope_check = check_scope_against_library(scope, library)
    logger.info("Derived library scope: %s (library %s declares %s)",
                scope.as_sorted(), library.key, list(library.positions))

    view = iso.CandidateView(
        variant_ids=list(universe.variant_ids),
        one_hot=universe.one_hot,
        index_of=dict(universe.index_of),
        combos=dict(universe.combos),
        library=library,
        batch_size=args.batch_size,
        n_rounds=args.rounds,
    )

    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    seeds = [args.base_seed + i for i in range(args.seeds)]

    production_facts: List[Dict[str, Any]] = []
    status: Dict[str, Dict[str, str]] = {}
    outcomes: List[iso.RunOutcome] = []

    for method in methods:
        if method == "panoptes_end_to_end":
            status[method] = {
                "status": "available",
                "detail": (
                    "real run_stream(): steps 1-6 including Gatekeeper authorisation, "
                    f"{'Modal GPU ESM + ' if not args.dry_run else ''}measured assay, ledger "
                    "evidence and Merkle chain verification"
                ),
            }
        elif method in ("random", "ridge", "random_forest"):
            status[method] = {"status": "available", "detail": f"baseline: {method}"}
        else:
            status[method] = {"status": "UNSUPPORTED", "detail": f"unknown method '{method}'"}

    workers = resolve_worker_count(args, len(seeds))
    if workers > 1:
        logger.info("Running %d seeds across %d worker processes", len(seeds), workers)
        limit_blas_threads()
        failed_seeds: List[int] = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_run_seed_worker, seed, methods, status, args): seed
                for seed in seeds
            }
            try:
                for future in concurrent.futures.as_completed(
                    futures, timeout=args.parallel_timeout
                ):
                    seed = futures[future]
                    try:
                        seed_outcomes, seed_prod = future.result()
                    except Exception as exc:
                        logger.error("seed %d failed in worker: %s: %s",
                                     seed, type(exc).__name__, exc)
                        failed_seeds.append(seed)
                        continue
                    outcomes.extend(seed_outcomes)
                    production_facts.extend(seed_prod)
                    logger.info("seed %d complete (%d outcomes)", seed, len(seed_outcomes))
            except concurrent.futures.TimeoutError:
                # Bounded: whatever has not finished by now is abandoned and reported
                # as missing rather than allowed to hang the run indefinitely.
                unfinished = [s for f, s in futures.items() if not f.done()]
                logger.error(
                    "parallel timeout after %ds; abandoning unfinished seed(s) %s",
                    args.parallel_timeout, unfinished,
                )
                for future in futures:
                    future.cancel()
                failed_seeds.extend(unfinished)
        if failed_seeds and not outcomes:
            # Every seed died. Previously this printed a 0/0 table and exited 0, which
            # reads exactly like a completed run with no advantage - the worst possible
            # failure mode for a benchmark. Refuse instead.
            raise E2EError(
                f"all {len(seeds)} seed(s) failed in worker processes, so no results were "
                f"produced. Seeds: {sorted(failed_seeds)}. A BrokenProcessPool here usually "
                f"means the workers exhausted RAM: each needs about "
                f"{WORKER_MEMORY_BUDGET_GB}GB for the eligible universe. Retry with "
                f"--parallel-seeds 1, or free memory."
            )
        if failed_seeds:
            logger.warning("%d of %d seeds failed: %s",
                           len(failed_seeds), len(seeds), sorted(failed_seeds))
    else:
        for seed in seeds:
            seed_outcomes, seed_prod = run_seed(
                ctx=SeedContext(universe, evaluator, view, scope, metadata,
                                wt_sequence, pdb_path),
                seed=seed, methods=methods, status=status, args=args,
            )
            outcomes.extend(seed_outcomes)
            production_facts.extend(seed_prod)

    outcomes.sort(key=lambda o: (o.method, o.seed))
    rows = iso.summarise(outcomes, status)
    production = summarise_production(production_facts)

    # Explicit vocabulary first, legacy string after it, exactly as the isolated harness
    # writes it, so a grep that finds one file's run folders finds the other's too.
    stamp = provenance_stamp(universe)
    iso.write_trace_csv(out_dir / "selection_trace.csv", outcomes, stamp)
    iso.write_round_metrics_csv(out_dir / "round_metrics.csv", outcomes, stamp)
    iso.write_summary_csv(out_dir / "summary.csv", rows, stamp)

    report = write_report(
        out_dir / "report.txt", args, universe, evaluator, scope, rows, production, status, started,
    )
    print("\n" + report)

    manifest = build_manifest(
        run_id=run_id,
        args=args,
        library=library,
        universe=universe,
        scope=scope,
        scope_check=scope_check,
        evaluator=evaluator,
        status=status,
        rows=rows,
        production=production,
        production_facts=production_facts,
        started=started,
    )
    validate_manifest(manifest)
    # No `default=str` fallback: a manifest that needs one is not machine-readable, and a
    # value that silently stringifies is a value no consumer can parse back. `_json_safe`
    # is the harness's own converter, so both benchmarks write the same shapes.
    (out_dir / "manifest.json").write_text(
        json.dumps(iso._json_safe(manifest), indent=2), encoding="utf-8",
    )

    print(f"Artifacts: {out_dir}")
    failed = [o for o in outcomes if o.status == "failed"]
    return 1 if failed else 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    iso._load_repo_dotenv()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    logging.getLogger("workers").setLevel(logging.WARNING)
    logging.getLogger("controller").setLevel(logging.INFO)

    for name, value in (("--seeds", args.seeds), ("--rounds", args.rounds),
                        ("--batch-size", args.batch_size), ("--initial-size", args.initial_size),
                        ("--parallel-seeds", args.parallel_seeds),
                        ("--parallel-timeout", args.parallel_timeout)):
        if value < 1:
            parser.error(f"{name} must be >= 1")

    # --dataset and --accession still carry the default library's values, so a user who
    # named a library without naming a file gets that library's own dataset and parent
    # rather than one library's variants scored against another's table.
    default_spec = iso.LIBRARIES[iso.DEFAULT_LIBRARY_KEY]
    if args.library != iso.DEFAULT_LIBRARY_KEY:
        chosen = iso.LIBRARIES[args.library]
        if args.dataset == str(default_spec.dataset):
            args.dataset = str(chosen.dataset)
        if args.accession == default_spec.accession:
            args.accession = chosen.accession

    try:
        return execute(args)
    except iso.ProvenanceError as exc:
        print(f"\nPROVENANCE ERROR\n{'=' * 78}\n{exc}\n", file=sys.stderr)
        return 3
    except paths.PathResolutionError as exc:
        # Printed rather than raised: the message already names every path that was
        # checked, and a traceback over the top of it only buries the useful part.
        print(f"\nPATH ERROR\n{'=' * 78}\n{exc}\n", file=sys.stderr)
        return 2
    except (E2EError, iso.BenchmarkError) as exc:
        print(f"\nERROR: {exc}\n", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

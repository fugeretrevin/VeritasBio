"""
benchmark_panoptes_isolated.py - Standalone comparison benchmark on the TrpB landscape.

Purpose:
    Compares PanoptesBio's production selection path against credible alternatives on the
    existing four-site TrpB active-site dataset (Johnston et al. 2024), under a single shared
    hidden-results protocol. Nothing in the existing project is imported for its side effects,
    and nothing in the existing project is modified: this file is self-contained, writes only
    into a fresh run folder under `benchmark_comparison_results/`, and imports production
    modules read-only.

Methods compared:
    A. random          - uniform sampling of untested eligible variants, without replacement.
    B. ridge           - one-hot (4 positions x 20 AA = 80 features) ridge regression, refit each round.
    C. random_forest   - same encoding, RandomForestRegressor, refit each round.
    D. alde            - the published ALDE implementation (jsunn-y/ALDE). Availability-probed;
                         marked UNAVAILABLE with a precise reason when the checkout or its
                         dependencies are absent. Never silently replaced by a homemade method.
    E. panoptes_default- the real production planning + selection path with its supported default
                         (non-Muse) strategy.
    F. panoptes_muse   - the same production path with live Muse planning and reflection feedback.
                         Requires explicit --enable-live-muse opt-in.

Protocol (shared by every method):
    - One evaluator owns the complete measured fitness table. Selectors never see it.
    - Per seed: the same 25 starting variants for every method, then 4 rounds of 25 newly
      selected variants. 125 unique measurements per method per seed, including the start.
    - A method commits its batch before the evaluator reveals those outcomes.
    - Duplicate, already-observed, unknown and invalid selections are rejected and logged.
      A short batch stays short; it is never backfilled with random picks.

Provenance (read this before quoting any number from this benchmark):
    The published TrpB four-site landscape is NOT fully measured. Johnston et al. (2024) state
    that 159,129 variants (99.45%) had sufficient sequencing coverage for quantification and
    that "the missing 871 fitness values were imputed" with the scikit-learn KNN imputer.
    The prepared CSV in this repository carries no flag distinguishing the two, so by default
    this benchmark refuses to run on it and tells you how to supply the imputed-variant list.
    See --provenance-mode.

Usage:
    python benchmark_panoptes_isolated.py --help
    python benchmark_panoptes_isolated.py --smoke-test
"""

from __future__ import annotations

# --- Isolation preamble: must precede every project import. --------------------------------
import sys

sys.dont_write_bytecode = True  # no __pycache__ anywhere, including under the project tree

import argparse
import csv
import hashlib
import json
import logging
import os
import platform
import random
import statistics
import subprocess
import time
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_DATASET = REPO_ROOT / "scripts" / "data" / "TRPB_THEMA_Johnston_2024.csv"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "benchmark_comparison_results"
PARENT_ACCESSION = "TRPB_TM9D8"
PARENT_METADATA_JSON = REPO_ROOT / "workers" / "data" / f"{PARENT_ACCESSION}.json"
STRUCTURE_PDB = REPO_ROOT / "workers" / "data" / "8EH0.pdb"

# The four saturated active-site positions of the published library (1-indexed, Tm9D8* numbering).
ACTIVE_SITE_POSITIONS: Tuple[int, ...] = (183, 184, 227, 228)
PARENT_RESIDUES = "VFVS"
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(AA_ALPHABET)}

# Published provenance of the upstream landscape (bioRxiv 10.1101/2024.06.23.600144v1).
PUBLISHED_MEASURED_COUNT = 159_129
PUBLISHED_IMPUTED_COUNT = 871
PUBLISHED_GLOBAL_OPTIMUM = "V183A:F184I:V227K:S228G"

# Higher fitness is better: the authors normalise each landscape's maximum fitness to 1.
FITNESS_DIRECTION = "higher_is_better"

BENCHMARK_VERSION = "1.0.0"

logger = logging.getLogger("panoptes_isolated_benchmark")


# ==============================================================================================
# Errors
# ==============================================================================================
class BenchmarkError(RuntimeError):
    """Base class for actionable benchmark configuration failures."""


class ProvenanceError(BenchmarkError):
    """Raised when measurement provenance cannot be established for the eligible dataset."""


# ==============================================================================================
# Dataset loading, provenance and the eligible universe
# ==============================================================================================
@dataclass(frozen=True)
class DatasetFingerprint:
    """Identity of the dataset file the run was scored against."""

    path: str
    sha256: str
    size_bytes: int
    rows_parsed: int
    unique_variants: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "path": self.path,
            "sha256": self.sha256,
            "size_bytes": self.size_bytes,
            "rows_parsed": self.rows_parsed,
            "unique_variants": self.unique_variants,
        }


@dataclass
class Universe:
    """
    The eligible candidate universe and its evaluator-only truth table.

    `fitness` is deliberately the ONLY place complete outcomes live. Selectors are handed
    `variant_ids` and `one_hot` (both sequence-derived) plus their own revealed measurements.
    """

    variant_ids: List[str]                      # canonical, position-ordered mutation notation
    combos: Dict[str, str]                      # variant id -> 4 residues at the saturated sites
    fitness: Dict[str, float]                   # EVALUATOR ONLY
    one_hot: np.ndarray                         # (n_variants, 80), row order == variant_ids
    index_of: Dict[str, int]
    fingerprint: DatasetFingerprint
    provenance_status: str
    provenance_detail: str
    excluded_imputed: int
    duplicates_dropped: int
    synthetic: bool = False

    def __len__(self) -> int:
        return len(self.variant_ids)


def _sha256_file(path: Path) -> Tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def _canonical_variant_id(combo: str) -> str:
    """Builds position-ordered ProteinGym notation for a 4-residue combo (parent residues dropped)."""
    parts = [
        f"{PARENT_RESIDUES[i]}{pos}{combo[i]}"
        for i, pos in enumerate(ACTIVE_SITE_POSITIONS)
        if combo[i] != PARENT_RESIDUES[i]
    ]
    return ":".join(parts)


def _combo_from_notation(notation: str) -> Optional[str]:
    """
    Maps 'V183A:F184I' to the 4-residue combo 'AIVS'. Returns None when the notation touches a
    position outside the saturated library or disagrees with the parent residue.
    """
    residues = list(PARENT_RESIDUES)
    if notation in ("", "WT"):
        return "".join(residues)
    for token in notation.split(":"):
        token = token.strip()
        if len(token) < 3:
            return None
        wt_aa, mut_aa, pos_text = token[0], token[-1], token[1:-1]
        if not pos_text.isdigit():
            return None
        position = int(pos_text)
        if position not in ACTIVE_SITE_POSITIONS:
            return None
        slot = ACTIVE_SITE_POSITIONS.index(position)
        if wt_aa != PARENT_RESIDUES[slot] or mut_aa not in AA_INDEX:
            return None
        residues[slot] = mut_aa
    return "".join(residues)


def _one_hot_matrix(combos: Sequence[str]) -> np.ndarray:
    """One-hot encodes combos as 4 positions x 20 amino acids = 80 binary features."""
    matrix = np.zeros((len(combos), len(ACTIVE_SITE_POSITIONS) * len(AA_ALPHABET)), dtype=np.float32)
    for row, combo in enumerate(combos):
        for slot, residue in enumerate(combo):
            matrix[row, slot * len(AA_ALPHABET) + AA_INDEX[residue]] = 1.0
    return matrix


def _load_imputed_ids(
    dataset_path: Path,
    explicit_path: Optional[Path],
) -> Tuple[Optional[Set[str]], str]:
    """
    Resolves the list of imputed (non-measured) variants, in order of preference:
      1. --imputed-ids FILE               (one variant id or 4-residue combo per line, or a CSV column)
      2. <dataset>.provenance.json        ({"imputed_variants": [...]})
    Returns (ids, source-description), or (None, reason) when no list could be resolved.
    """
    def _parse_lines(path: Path) -> Set[str]:
        ids: Set[str] = set()
        with open(path, "r", encoding="utf-8") as handle:
            for raw in handle:
                token = raw.strip().split(",")[0].strip().strip('"')
                if not token or token.lower() in ("mutant", "variant", "variant_id", "combo"):
                    continue
                ids.add(token)
        return ids

    if explicit_path is not None:
        if not explicit_path.is_file():
            raise ProvenanceError(f"--imputed-ids file not found: {explicit_path}")
        if explicit_path.suffix.lower() == ".json":
            payload = json.loads(explicit_path.read_text(encoding="utf-8"))
            ids = set(payload.get("imputed_variants") or payload.get("imputed") or [])
        else:
            ids = _parse_lines(explicit_path)
        return ids, f"--imputed-ids {explicit_path}"

    sidecar = dataset_path.with_suffix(dataset_path.suffix + ".provenance.json")
    if sidecar.is_file():
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
        ids = set(payload.get("imputed_variants") or payload.get("imputed") or [])
        return ids, f"sidecar {sidecar.name}"

    return None, "no imputed-variant list supplied and the dataset carries no provenance column"


PROVENANCE_HELP = f"""
Measurement provenance could not be established for this dataset.

Why this matters:
  Johnston et al. (2024) report that {PUBLISHED_MEASURED_COUNT:,} of the 160,000 four-site TrpB
  variants had sufficient sequencing coverage to be quantified, and that the remaining
  {PUBLISHED_IMPUTED_COUNT} fitness values were IMPUTED with the scikit-learn KNN imputer.
  The prepared CSV in this repository re-encodes the published `label` column verbatim and
  carries no flag separating measured from imputed records, so a run over it would mix
  {PUBLISHED_IMPUTED_COUNT} model-filled values into a benchmark presented as experimental.

How to fix it (pick one):
  1. Supply the imputed-variant list:
       --imputed-ids <file>
     The file may be one variant per line (either 'V183A:F184I:V227K:S228G' notation or a
     4-residue combo such as 'AIKG'), a single-column CSV, or a JSON object shaped
     {{"imputed_variants": [...]}}. The published list is derivable from the authors'
     supplementary material (SI Appendix, Fig. S16) or the upstream per-replicate read counts.
  2. Drop a sidecar next to the dataset named '<dataset>.csv.provenance.json' with the same
     JSON shape; it is picked up automatically.
  3. Add a provenance column to a copy of the dataset: one of 'provenance', 'is_imputed',
     'imputed' or 'measured', and point --dataset at the copy.
  4. Acknowledge the limitation explicitly and proceed anyway:
       --provenance-mode acknowledge-unflagged
     Every output file, the plot and the report are then stamped
     'provenance=UNVERIFIED_INCLUDES_IMPUTED'. Results from such a run must not be described
     as a clean experimentally-grounded comparison.
"""


def load_universe(
    dataset_path: Path,
    provenance_mode: str,
    imputed_ids_path: Optional[Path],
    max_variants: Optional[int] = None,
    subsample_seed: int = 0,
) -> Universe:
    """
    Loads the TrpB CSV into an eligible universe, enforcing the provenance policy.

    Filtering rules, all reported in the manifest:
      - a record must map onto the four saturated positions and agree with the parent residues;
      - duplicate variant ids keep the first occurrence;
      - non-finite scores are dropped;
      - imputed records are dropped when a provenance list is available.
    """
    if not dataset_path.is_file():
        raise BenchmarkError(
            f"Dataset not found: {dataset_path}\n"
            "Prepare it first with:  python scripts/download_trpb_benchmark.py"
        )

    sha, size = _sha256_file(dataset_path)

    rows: List[Tuple[str, float, Optional[str]]] = []
    provenance_column: Optional[str] = None
    with open(dataset_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = [f.strip() for f in (reader.fieldnames or [])]
        if "mutant" not in fieldnames or "DMS_score" not in fieldnames:
            raise BenchmarkError(
                f"{dataset_path} does not look like the prepared TrpB benchmark CSV "
                f"(expected 'mutant' and 'DMS_score' columns, found {fieldnames})."
            )
        for candidate in ("provenance", "is_imputed", "imputed", "measured"):
            if candidate in fieldnames:
                provenance_column = candidate
                break
        for row in reader:
            notation = (row.get("mutant") or "").strip()
            raw_score = (row.get("DMS_score") or "").strip()
            if not notation or not raw_score:
                continue
            try:
                score = float(raw_score)
            except ValueError:
                continue
            if not np.isfinite(score):
                continue
            flag = (row.get(provenance_column) or "").strip() if provenance_column else None
            rows.append((notation, score, flag))

    if not rows:
        raise BenchmarkError(f"No usable rows parsed from {dataset_path}")

    # ---- Provenance resolution -------------------------------------------------------------
    imputed_lookup: Optional[Set[str]] = None
    provenance_detail: str
    if provenance_column:
        provenance_status = "VERIFIED_COLUMN"
        provenance_detail = f"dataset column '{provenance_column}'"
    else:
        imputed_lookup, source = _load_imputed_ids(dataset_path, imputed_ids_path)
        if imputed_lookup is not None:
            provenance_status = "VERIFIED_EXTERNAL_LIST"
            provenance_detail = f"{len(imputed_lookup)} imputed ids from {source}"
        elif provenance_mode == "acknowledge-unflagged":
            provenance_status = "UNVERIFIED_INCLUDES_IMPUTED"
            provenance_detail = (
                f"{source}; up to {PUBLISHED_IMPUTED_COUNT} published KNN-imputed values "
                "are included because --provenance-mode acknowledge-unflagged was passed"
            )
        else:
            raise ProvenanceError(PROVENANCE_HELP.strip())

    # ---- Build the universe ----------------------------------------------------------------
    seen: Set[str] = set()
    duplicates = 0
    excluded_imputed = 0
    unmapped = 0
    ids: List[str] = []
    combos: List[str] = []
    scores: List[float] = []

    for notation, score, flag in rows:
        combo = _combo_from_notation(notation)
        if combo is None:
            unmapped += 1
            continue
        variant_id = _canonical_variant_id(combo)
        if not variant_id:
            continue  # the parent itself is not a library member
        if variant_id in seen:
            duplicates += 1
            continue
        if provenance_column:
            is_imputed = flag is not None and (
                flag.lower() in ("1", "true", "yes", "imputed")
                if provenance_column != "measured"
                else flag.lower() in ("0", "false", "no")
            )
            if is_imputed:
                excluded_imputed += 1
                continue
        elif imputed_lookup is not None and (variant_id in imputed_lookup or combo in imputed_lookup):
            excluded_imputed += 1
            continue
        seen.add(variant_id)
        ids.append(variant_id)
        combos.append(combo)
        scores.append(score)

    if unmapped:
        logger.warning("%d record(s) did not map onto positions %s and were dropped.", unmapped, list(ACTIVE_SITE_POSITIONS))
    if len(ids) < 1000:
        raise BenchmarkError(
            f"Only {len(ids)} eligible variants survived filtering; refusing to benchmark on so few records."
        )

    if max_variants is not None and max_variants < len(ids):
        rng = random.Random(subsample_seed)
        keep = sorted(rng.sample(range(len(ids)), max_variants))
        ids = [ids[i] for i in keep]
        combos = [combos[i] for i in keep]
        scores = [scores[i] for i in keep]
        provenance_detail += f"; candidate universe subsampled to {max_variants} variants (seed {subsample_seed})"

    order = sorted(range(len(ids)), key=lambda i: ids[i])
    ids = [ids[i] for i in order]
    combos = [combos[i] for i in order]
    scores = [scores[i] for i in order]

    universe = Universe(
        variant_ids=ids,
        combos=dict(zip(ids, combos)),
        fitness=dict(zip(ids, scores)),
        one_hot=_one_hot_matrix(combos),
        index_of={v: i for i, v in enumerate(ids)},
        fingerprint=DatasetFingerprint(
            path=str(dataset_path),
            sha256=sha,
            size_bytes=size,
            rows_parsed=len(rows),
            unique_variants=len(ids),
        ),
        provenance_status=provenance_status,
        provenance_detail=provenance_detail,
        excluded_imputed=excluded_imputed,
        duplicates_dropped=duplicates,
    )

    best_id = max(universe.fitness, key=lambda v: universe.fitness[v])
    if best_id != PUBLISHED_GLOBAL_OPTIMUM:
        logger.warning(
            "Global optimum in the eligible set is %s, not the published %s. Check the dataset and "
            "the fitness direction before interpreting results.",
            best_id, PUBLISHED_GLOBAL_OPTIMUM,
        )
    return universe


def build_synthetic_universe(n_variants: int = 0, seed: int = 7) -> Universe:
    """
    Builds a clearly-labelled SYNTHETIC four-site landscape for mechanism validation.

    Values here are generated additive-plus-pairwise numbers, NOT biology. Every output of a
    synthetic run is stamped so these numbers can never be mistaken for assay performance.

    `n_variants <= 0` (the default) covers the complete 20^4 - 1 combo space, which is what the
    production adapter needs: it proposes variants from the full saturated library, so a sparse
    synthetic subsample would reject nearly every one of its selections and leave that code path
    untested. A positive value subsamples instead, for a faster but less representative check.
    """
    rng = np.random.default_rng(seed)
    singles = rng.normal(0.0, 0.25, size=(len(ACTIVE_SITE_POSITIONS), len(AA_ALPHABET)))
    pair = rng.normal(0.0, 0.12, size=(len(ACTIVE_SITE_POSITIONS), len(AA_ALPHABET),
                                       len(ACTIVE_SITE_POSITIONS), len(AA_ALPHABET)))

    complete = [
        a + b + c + d
        for a in AA_ALPHABET for b in AA_ALPHABET for c in AA_ALPHABET for d in AA_ALPHABET
        if a + b + c + d != PARENT_RESIDUES
    ]
    if n_variants and n_variants < len(complete):
        combos: Sequence[str] = random.Random(seed).sample(complete, n_variants)
    else:
        combos = complete

    ids: List[str] = []
    combo_list: List[str] = []
    scores: List[float] = []
    for combo in sorted(combos):
        value = 0.0
        for i, residue_i in enumerate(combo):
            value += singles[i, AA_INDEX[residue_i]]
            for j, residue_j in enumerate(combo):
                if j > i:
                    value += pair[i, AA_INDEX[residue_i], j, AA_INDEX[residue_j]]
        ids.append(_canonical_variant_id(combo))
        combo_list.append(combo)
        scores.append(float(value))

    order = sorted(range(len(ids)), key=lambda i: ids[i])
    ids = [ids[i] for i in order]
    combo_list = [combo_list[i] for i in order]
    scores = [scores[i] for i in order]

    return Universe(
        variant_ids=ids,
        combos=dict(zip(ids, combo_list)),
        fitness=dict(zip(ids, scores)),
        one_hot=_one_hot_matrix(combo_list),
        index_of={v: i for i, v in enumerate(ids)},
        fingerprint=DatasetFingerprint(
            path=f"<synthetic:seed={seed}>",
            sha256=hashlib.sha256(f"synthetic-{seed}-{n_variants}".encode()).hexdigest(),
            size_bytes=0,
            rows_parsed=len(ids),
            unique_variants=len(ids),
        ),
        provenance_status="SYNTHETIC_NOT_BIOLOGICAL",
        provenance_detail="generated additive + pairwise values for mechanism validation only",
        excluded_imputed=0,
        duplicates_dropped=0,
        synthetic=True,
    )


# ==============================================================================================
# Evaluator and per-method observation state
# ==============================================================================================
@dataclass
class RevealedState:
    """A single method's private view: only what the evaluator has already revealed to it."""

    measured: Dict[str, float] = field(default_factory=dict)

    @property
    def ids(self) -> List[str]:
        return list(self.measured.keys())

    def best(self) -> float:
        return max(self.measured.values()) if self.measured else float("nan")


class Evaluator:
    """
    Owns the complete measured fitness table and the top-1% definition.

    Selection methods receive this object's `reveal()` output only after they have committed a
    batch. The threshold and top-1% membership stay evaluator-side and are never passed to a
    selector, logged into selector-visible structures, or exposed through candidate metadata.
    """

    def __init__(self, universe: Universe) -> None:
        self._universe = universe
        self.harness_exports: List[str] = []
        values = np.array([universe.fitness[v] for v in universe.variant_ids], dtype=np.float64)
        # Defined once over the eligible measured dataset. Ties: every variant at or above the
        # threshold counts as a hit, so the hit set can exceed 1% of the pool when values repeat.
        self.top1_threshold = float(np.percentile(values, 99.0))
        self._top1: Set[str] = {v for v in universe.variant_ids if universe.fitness[v] >= self.top1_threshold}
        self.global_best = float(values.max())
        self.n_top1 = len(self._top1)

    @property
    def eligible_ids(self) -> List[str]:
        return self._universe.variant_ids

    def reveal(self, variant_ids: Sequence[str]) -> Dict[str, float]:
        """Reveals measured outcomes for an already-committed batch."""
        return {v: self._universe.fitness[v] for v in variant_ids}

    def is_top1(self, variant_id: str) -> bool:
        return variant_id in self._top1

    def export_truth_table(self, method: str) -> Dict[str, float]:
        """
        Hands the complete fitness table to a method that runs its own simulation loop.

        This is the one documented exception to the hidden-results separation, and it is
        recorded here so it shows up in the manifest rather than happening quietly. It exists
        because ALDE is a simulation harness: its loop owns the lookup table and reads an entry
        only after its acquisition function has committed to that point, which preserves the
        commit-then-reveal order even though the array is in its memory. The top-1% threshold
        and membership are NOT exported, and this benchmark re-scores every returned pick with
        its own evaluator.
        """
        self.harness_exports.append(method)
        logger.info("Evaluator exported the full fitness table to harness method '%s'.", method)
        return dict(self._universe.fitness)


@dataclass
class CandidateView:
    """Everything a selector is allowed to see. Deliberately free of outcome information."""

    variant_ids: List[str]
    one_hot: np.ndarray
    index_of: Dict[str, int]
    combos: Dict[str, str]
    batch_size: int
    n_rounds: int

    def untested(self, revealed: RevealedState) -> List[str]:
        observed = set(revealed.measured)
        return [v for v in self.variant_ids if v not in observed]


# ==============================================================================================
# Selection methods
# ==============================================================================================
class Selector:
    """Interface every method implements. `propose` must not consult the evaluator."""

    name = "selector"
    version = BENCHMARK_VERSION
    config: Dict[str, Any] = {}

    def __init__(self, seed: int) -> None:
        self.seed = seed
        self.model_calls = 0
        self.token_usage: Optional[int] = None  # None == unknown, never guessed
        self.notes: List[str] = []

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        raise NotImplementedError


def _rank_select(
    predictions: np.ndarray,
    candidate_ids: Sequence[str],
    n: int,
) -> List[str]:
    """
    Picks the n highest-predicted candidates.

    Tie-breaking is reproducible: ties are broken by the candidate's position in the universe's
    fixed lexicographic order, which `candidate_ids` preserves. np.lexsort is stable on that key.
    """
    order = np.lexsort((np.arange(len(candidate_ids)), -predictions))
    return [candidate_ids[int(i)] for i in order[:n]]


class RandomSelector(Selector):
    """A. Uniform random sampling of untested eligible variants, without replacement."""

    name = "random"
    config = {"strategy": "uniform without replacement"}

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        pool = view.untested(revealed)
        rng = random.Random((self.seed + 1) * 100_003 + round_index)
        return rng.sample(pool, min(n, len(pool)))


class RidgeSelector(Selector):
    """
    B. One-hot ridge regression over the 80 binary features, refit on revealed measurements
    each round, selecting the highest predicted untested variants.

    alpha=1.0 is scikit-learn's documented default and was NOT tuned against held-out results.
    """

    name = "ridge"
    config = {"model": "sklearn.linear_model.Ridge", "alpha": 1.0, "encoding": "onehot_4x20"}

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        from sklearn.linear_model import Ridge

        pool = view.untested(revealed)
        if not revealed.measured:
            return RandomSelector(self.seed).propose(round_index, n, view, revealed)
        train_rows = [view.index_of[v] for v in revealed.measured]
        model = Ridge(alpha=1.0)
        model.fit(view.one_hot[train_rows], np.array([revealed.measured[v] for v in revealed.measured]))
        self.model_calls += 1
        predictions = model.predict(view.one_hot[[view.index_of[v] for v in pool]])
        return _rank_select(predictions, pool, min(n, len(pool)))


class RandomForestSelector(Selector):
    """
    C. One-hot random forest over the same encoding and the same revealed observations,
    refit each round. Nonlinear, so it can represent interactions between the four sites.

    Settings are scikit-learn defaults apart from n_estimators=300 (variance reduction) and a
    per-seed random_state for reproducibility. Nothing here was tuned on hidden results.
    """

    name = "random_forest"
    config = {
        "model": "sklearn.ensemble.RandomForestRegressor",
        "n_estimators": 300,
        "max_depth": None,
        "min_samples_leaf": 1,
        "encoding": "onehot_4x20",
    }

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        from sklearn.ensemble import RandomForestRegressor

        pool = view.untested(revealed)
        if not revealed.measured:
            return RandomSelector(self.seed).propose(round_index, n, view, revealed)
        train_rows = [view.index_of[v] for v in revealed.measured]
        model = RandomForestRegressor(
            n_estimators=300,
            random_state=self.seed,
            n_jobs=-1,
        )
        model.fit(view.one_hot[train_rows], np.array([revealed.measured[v] for v in revealed.measured]))
        self.model_calls += 1
        predictions = model.predict(view.one_hot[[view.index_of[v] for v in pool]])
        return _rank_select(predictions, pool, min(n, len(pool)))


# ----------------------------------------------------------------------------------------------
# D. Published ALDE
# ----------------------------------------------------------------------------------------------
ALDE_REQUIRED_MODULES = ("torch", "gpytorch", "botorch", "xgboost", "sklearn", "pandas")
ALDE_REQUIRED_FILES = ("execute_simulation.py", "src/optimize.py", "src/objectives.py",
                       "src/acquisition.py", "src/encoding_utils.py")

# Worker executed by the isolated ALDE interpreter. It is written to the isolated workspace at
# run time so that the ALDE checkout itself stays byte-for-byte unmodified, and so that all
# integration logic lives in this one benchmark file.
#
# Faithfulness to the authors' simulation (execute_simulation.py, commit d0b3593):
#   - encoding   : src.encoding_utils.generate_onehot, their own function, reshaped to [n, 80]
#   - objective  : src.objectives.Combo, unmodified (it max-normalises y exactly as in their runs)
#   - model      : DNN_ENSEMBLE with architecture [d, 30, 30, 1], the branch their script selects
#                  for onehot encodings; lr 1e-3, activation lrelu, train_iter 300, dropout 0
#   - acquisition: TS (Thompson sampling)
#   - loop       : src.optimize.BayesianOptimization.run, unmodified
# The only substitution is the initial sample: their utils.samp_discrete(n, obj, seed) draw is
# replaced by this benchmark's shared starting variants, so every method starts from the same 25.
ALDE_WORKER_SOURCE = '''
"""Generated by benchmark_panoptes_isolated.py. Runs the authors' ALDE loop unmodified."""
import json, os, sys, time, warnings
warnings.filterwarnings("ignore")

payload_path, result_path = sys.argv[1], sys.argv[2]
with open(payload_path, "r", encoding="utf-8") as handle:
    cfg = json.load(handle)

sys.path.insert(0, cfg["alde_checkout"])
os.chdir(cfg["work_dir"])

import numpy as np
import pandas as pd
import torch
import random

from src.encoding_utils import generate_onehot
from src.optimize import BayesianOptimization, BO_ARGS
import src.objectives as objectives

combos = np.load(cfg["data_npz"], allow_pickle=True)
combo_list = [str(c) for c in combos["combos"]]
fitness = combos["fitness"].astype(np.float64)

# Lay the data out the way Combo expects, inside the isolated work dir.
data_dir = os.path.join("data", cfg["protein"])
os.makedirs(data_dir, exist_ok=True)
pd.DataFrame({"Combo": combo_list, "fitness": fitness}).to_csv(
    os.path.join(data_dir, "fitness.csv"), index=False)
if not os.path.exists(os.path.join(data_dir, "onehot_x.pt")):
    encoded = generate_onehot(combo_list)
    encoded = torch.reshape(encoded, (encoded.shape[0], -1))
    torch.save(encoded, os.path.join(data_dir, "onehot_x.pt"))

obj = objectives.Combo(cfg["protein"], "onehot")
domain = obj.get_domain()
disc_X = obj.get_points()[0]
disc_y = obj.get_points()[1]

seed = int(cfg["seed"])
torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

# The shared starting sample replaces utils.samp_discrete, so all methods begin identically.
start_indices = torch.tensor(cfg["start_indices"], dtype=torch.long)
start_x = disc_X[start_indices]
start_y = disc_y[start_indices]

arc = [domain[0].size(-1), 30, 30, 1]   # the onehot branch of execute_simulation.py
savedir = os.path.join(cfg["savedir_subdir"], "run_")
os.makedirs(os.path.dirname(savedir), exist_ok=True)

args = BO_ARGS(
    mtype="DNN_ENSEMBLE",
    kernel="RBF",
    acq_fn="TS",
    xi=4,
    architecture=arc,
    activation="lrelu",
    min_noise=1e-6,
    trainlr=1e-3,
    train_iter=int(cfg["train_iter"]),
    dropout=0,
    mcdropout=0,
    verbose=int(cfg["verbose"]),
    bb_fn=obj.objective,
    domain=domain,
    disc_X=disc_X,
    disc_y=disc_y,
    noise_std=0,
    n_rand_init=0,
    budget=int(cfg["budget"]),
    query_cost=1,
    queries_x=start_x,
    queries_y=start_y,
    indices=start_indices,
    savedir=savedir,
    batch_size=int(cfg["batch_size"]),
)

started = time.time()
BayesianOptimization.run(args, seed)
runtime = time.time() - started

selected = torch.load(savedir + "indices.pt").reshape(-1).tolist()
n_initial = len(cfg["start_indices"])
queried = [int(i) for i in selected[n_initial:]]

with open(result_path, "w", encoding="utf-8") as handle:
    json.dump({
        "queried_indices": queried,
        "n_initial": n_initial,
        "runtime_seconds": runtime,
        "versions": {
            "torch": torch.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }, handle)
print("ALDE_WORKER_OK", len(queried))
'''


def probe_alde(checkout: Optional[Path], python_exe: Optional[Path]) -> Dict[str, Any]:
    """
    Reports precisely what is and is not available for running the authors' ALDE implementation
    (https://github.com/jsunn-y/ALDE, Nat Commun 2025, doi:10.1038/s41467-025-55987-8).

    Dependencies are probed in ALDE's OWN interpreter, not this one: the two environments are
    deliberately separate, so importing torch here would prove nothing. No part of ALDE is
    reimplemented: when something is missing, the method is reported UNAVAILABLE with the
    specific cause (checkout, dependency, data or adapter) and nothing is substituted.
    """
    report: Dict[str, Any] = {
        "checkout_path": str(checkout) if checkout else None,
        "python_exe": str(python_exe) if python_exe else None,
        "checkout_present": False,
        "interpreter_present": False,
        "missing_files": [],
        "dependencies": {},
        "missing_dependencies": [],
        "imports_alde": None,
        "import_error": None,
        "commit": None,
        "available": False,
        "reason": "",
        "pinned_versions": "alde.yml pins python=3.11, pytorch=2.1.1, gpytorch==1.11, botorch==0.9.4",
    }

    if checkout is None:
        report["reason"] = (
            "No ALDE checkout supplied. Pass --alde-path /path/to/ALDE (clone of "
            "https://github.com/jsunn-y/ALDE) or set ALDE_PATH."
        )
        return report
    if not checkout.is_dir():
        report["reason"] = f"--alde-path {checkout} is not a directory."
        return report

    report["checkout_present"] = True
    report["missing_files"] = [f for f in ALDE_REQUIRED_FILES if not (checkout / f).exists()]
    try:
        report["commit"] = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(checkout),
            capture_output=True, text=True, timeout=15,
        ).stdout.strip() or None
    except Exception:
        report["commit"] = None

    if python_exe is None or not Path(python_exe).is_file():
        report["reason"] = (
            f"ALDE interpreter not found at {python_exe}. Create an isolated environment for the "
            "checkout and point --alde-python at its python executable."
        )
        return report
    report["interpreter_present"] = True

    probe_source = (
        "import json,sys\n"
        "out={'versions':{},'missing':[],'import_error':None}\n"
        "for m in ('torch','gpytorch','botorch','xgboost','sklearn','pandas','numpy'):\n"
        "    try:\n"
        "        out['versions'][m]=getattr(__import__(m),'__version__','unknown')\n"
        "    except Exception as e:\n"
        "        out['missing'].append(m); out['versions'][m]=None\n"
        "try:\n"
        "    import src.optimize, src.objectives, src.acquisition, src.encoding_utils\n"
        "    out['imports_alde']=True\n"
        "except Exception as e:\n"
        "    out['imports_alde']=False; out['import_error']='%s: %s' % (type(e).__name__, e)\n"
        "print(json.dumps(out))\n"
    )
    try:
        completed = subprocess.run(
            [str(python_exe), "-c", probe_source],
            cwd=str(checkout), capture_output=True, text=True, timeout=300,
        )
        parsed = json.loads(completed.stdout.strip().splitlines()[-1])
        report["dependencies"] = parsed["versions"]
        report["missing_dependencies"] = parsed["missing"]
        report["imports_alde"] = parsed["imports_alde"]
        report["import_error"] = parsed["import_error"]
    except Exception as exc:
        report["reason"] = f"could not probe the ALDE interpreter: {type(exc).__name__}: {exc}"
        return report

    problems = []
    if report["missing_files"]:
        problems.append(f"missing repository file(s): {report['missing_files']}")
    if report["missing_dependencies"]:
        problems.append("missing Python dependenc(ies) in the ALDE environment: "
                        + ", ".join(report["missing_dependencies"]))
    if report["imports_alde"] is False:
        problems.append(f"ALDE's own modules do not import: {report['import_error']}")
    report["reason"] = "; ".join(problems) or "checkout, interpreter and dependencies all present"
    report["available"] = not problems
    return report


class AldeSelector(Selector):
    """
    D. The authors' ALDE implementation, run unmodified in its own interpreter.

    ALDE owns its active-learning loop: BayesianOptimization.optimize() retrains and then picks
    `batch_size` points per round until `budget` queries are spent. So this adapter runs one
    campaign per seed through `run_campaign` rather than proposing round by round, and this
    benchmark's evaluator independently re-scores whatever ALDE picked.

    Hidden-results discipline inside ALDE (verified in src/optimize.py at commit d0b3593):
      - the surrogate trains only on `queries_x` / `norm_y`, the points already queried;
      - `disc_y[idx]` is read in `update_trajectory` only AFTER a point has been selected, which
        is the same commit-then-reveal order this benchmark enforces for every other method;
      - the full `disc_y` is otherwise used only for regret logging (`obj_max`).

    Budget convention: ALDE's `budget` EXCLUDES the initial sample, so 25 shared starting
    measurements plus budget=100 gives the same 125 total measurements as every other method.
    """

    name = "alde"

    def __init__(
        self,
        seed: int,
        checkout: Path,
        python_exe: Path,
        workspace: Path,
        probe: Dict[str, Any],
        timeout_seconds: int,
        train_iter: int = 300,
        verbose: int = 1,
    ) -> None:
        super().__init__(seed)
        self.checkout = checkout
        self.python_exe = python_exe
        self.workspace = workspace
        self.probe = probe
        self.timeout_seconds = timeout_seconds
        self.train_iter = train_iter
        self.verbose = verbose
        self._truth: Optional[Dict[str, float]] = None
        self.config = {
            "source": "https://github.com/jsunn-y/ALDE",
            "commit": probe.get("commit"),
            "paper": "doi:10.1038/s41467-025-55987-8",
            "encoding": "onehot (src.encoding_utils.generate_onehot)",
            "model": "DNN_ENSEMBLE, architecture [d, 30, 30, 1]",
            "acquisition": "TS",
            "train_iter": train_iter,
            "lr": 1e-3,
            "activation": "lrelu",
            "budget_convention": "budget excludes the initial sample",
            "execution": "subprocess, separate interpreter",
            "interpreter": str(python_exe),
            "dependency_versions": probe.get("dependencies"),
            "deviation_from_alde_yml": probe.get("deviation", ""),
        }

    def run_campaign(
        self,
        view: CandidateView,
        initial_ids: Sequence[str],
        n_rounds: int,
        batch_size: int,
    ) -> List[List[str]]:
        """Runs one full ALDE campaign and returns its picks split into per-round batches."""
        # The encoded design space is identical across seeds, so it is built once in a shared
        # work dir; only each seed's outputs and logs are per-run. Without this, every seed
        # would rewrite the same ~100 MB one-hot tensor.
        work_dir = self.workspace / "alde_work"
        run_root = work_dir / "runs" / f"seed_{self.seed}_{uuid.uuid4().hex[:6]}"
        run_root.mkdir(parents=True, exist_ok=True)

        data_npz = work_dir / "universe.npz"
        universe_key = hashlib.sha256(
            "".join(view.variant_ids).encode()
        ).hexdigest()
        key_path = work_dir / "universe.key"
        if not data_npz.is_file() or not key_path.is_file() or key_path.read_text().strip() != universe_key:
            np.savez_compressed(
                data_npz,
                combos=np.array([view.combos[v] for v in view.variant_ids], dtype=object),
                fitness=np.array([self._truth_placeholder(v) for v in view.variant_ids], dtype=np.float64),
            )
            key_path.write_text(universe_key, encoding="utf-8")
            for stale in (work_dir / "data").glob("**/onehot_x.pt"):
                stale.unlink()  # the design space changed, so the cached encoding is invalid
        payload = {
            "alde_checkout": str(self.checkout),
            "work_dir": str(work_dir),
            "savedir_subdir": str(run_root),
            "data_npz": str(data_npz),
            "protein": "TrpB_panoptes",
            "seed": self.seed,
            "budget": n_rounds * batch_size,   # excludes the initial sample, per ALDE's convention
            "batch_size": batch_size,
            "train_iter": self.train_iter,
            "verbose": self.verbose,
            "start_indices": [view.index_of[v] for v in initial_ids],
        }
        payload_path = run_root / "payload.json"
        payload_path.write_text(json.dumps(payload), encoding="utf-8")
        worker_path = run_root / "alde_worker_generated.py"
        worker_path.write_text(ALDE_WORKER_SOURCE, encoding="utf-8")
        result_path = run_root / "result.json"

        environment = dict(os.environ)
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        environment["OMP_NUM_THREADS"] = environment.get("OMP_NUM_THREADS", "4")
        started = time.perf_counter()
        completed = subprocess.run(
            [str(self.python_exe), str(worker_path), str(payload_path), str(result_path)],
            cwd=str(run_root), capture_output=True, text=True,
            timeout=self.timeout_seconds, env=environment,
        )
        elapsed = time.perf_counter() - started
        (run_root / "stdout.log").write_text(completed.stdout or "", encoding="utf-8")
        (run_root / "stderr.log").write_text(completed.stderr or "", encoding="utf-8")

        if completed.returncode != 0 or not result_path.is_file():
            tail = (completed.stderr or completed.stdout or "").strip().splitlines()[-12:]
            raise BenchmarkError(
                f"ALDE subprocess failed (exit {completed.returncode}) after {elapsed:.1f}s. "
                f"Logs: {run_root}. Last output:\n  " + "\n  ".join(tail)
            )

        result = json.loads(result_path.read_text(encoding="utf-8"))
        queried = result["queried_indices"]
        self.model_calls += n_rounds  # one surrogate refit per round
        self.notes.append(
            f"ALDE commit {self.probe.get('commit', 'unknown')[:12]} returned {len(queried)} "
            f"queries in {result['runtime_seconds']:.1f}s (torch {result['versions']['torch']})"
        )
        batches: List[List[str]] = []
        for start in range(0, n_rounds * batch_size, batch_size):
            chunk = queried[start:start + batch_size]
            batches.append([
                view.variant_ids[i] for i in chunk if 0 <= i < len(view.variant_ids)
            ])
        return batches

    def _truth_placeholder(self, variant_id: str) -> float:
        """
        Supplies the measured value for ALDE's `disc_y`.

        ALDE is a simulation harness: its loop reads `disc_y[idx]` to reveal the outcome of a
        point it has just chosen, exactly as this benchmark's evaluator reveals a committed
        batch. The surrogate never trains on unqueried entries. The injected table is set by
        `attach_truth` before the campaign runs.
        """
        if self._truth is None:
            raise BenchmarkError("ALDE campaign started before the evaluator exported its table.")
        return self._truth[variant_id]

    def attach_truth(self, truth: Dict[str, float]) -> None:
        """Called by the runner only for harness-style methods that run their own loop."""
        self._truth = truth

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        raise BenchmarkError("ALDE runs as a campaign; run_campaign() should have been used.")


# ----------------------------------------------------------------------------------------------
# E/F. The existing PanoptesBio production path
# ----------------------------------------------------------------------------------------------
class _LiveCallRecorder(logging.Handler):
    """
    Observes the orchestrator's own log stream to record whether Muse calls actually went live.

    Production falls back to the mock planner when the Meta API call fails. That fallback is
    silent to the caller, so a run could otherwise be mislabelled 'live'. This handler only
    reads log records; it changes no production behaviour.
    """

    def __init__(self) -> None:
        super().__init__(level=logging.INFO)
        self.live_attempts = 0
        self.fallbacks = 0
        self.messages: List[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        try:
            message = record.getMessage()
        except Exception:
            return
        if "Querying Meta Model API" in message:
            self.live_attempts += 1
        if "Meta Model API" in message and "failed" in message:
            self.fallbacks += 1
            self.messages.append(message[:200])


class _AnnotationDrivenPlanner:
    """
    The non-live planner used by `panoptes_default`.

    Why this exists instead of `workers.mock_workers.MockLLMWorker`:
        MockLLMWorker is a scripted demo planner. For any accession matching TrpB it
        returns the literal positions 183/184/227/228 - exactly the four saturated
        sites of this benchmark's library - and `suggested_chemistries` naming three
        of the four residues of the published global optimum V183A:F184I:V227K:S228G
        (183:[A,V,C], 184:[I,V,L], 228:[G,A,S]). That is appropriate for a product
        demo and disqualifying for a baseline: it encodes the answer this benchmark
        is supposed to measure the discovery of.

    What this does instead:
        Derives the search space from the target's own annotations, which is how a
        real campaign knows what its ordered library can assay. Positions come from
        `metadata.binding_positions` (populated from the entry's annotated active
        sites), catalytic machinery from `metadata.catalytic_positions` is forbidden,
        and chemistries are a fixed neutral alphabet-spanning set chosen for physical
        diversity, never for fitness on this landscape. Nothing here mentions TrpB,
        so the same code path serves any target with annotated sites.
    """

    # A neutral, target-independent probe set spanning small/polar/charged/branched/
    # aromatic classes. Deliberately NOT tuned to any landscape: it must not encode
    # which residues win. GFlowNet sampling in step 3 explores beyond these.
    NEUTRAL_CHEMISTRIES: Tuple[str, ...] = ("A", "S", "D", "K", "L", "F")

    def generate_strategy(
        self,
        metadata: Any,
        round_number: int = 1,
        top_evidence: Optional[Sequence[Dict[str, Any]]] = None,
        current_hop: int = 0,
    ) -> Any:
        from controller.prompts import BiophysicalHypothesis, ResidueCluster, SearchStrategy

        annotated = sorted({int(p) for p in (getattr(metadata, "binding_positions", None) or [])})
        catalytic = sorted({int(p) for p in (getattr(metadata, "catalytic_positions", None) or [])})
        targetable = [p for p in annotated if p not in set(catalytic)]
        if not targetable:
            raise BenchmarkError(
                "annotation-driven planner found no targetable annotated positions for "
                f"{getattr(metadata, 'accession', '<unknown>')}: binding_positions={annotated}, "
                f"catalytic_positions={catalytic}. panoptes_default cannot plan without "
                "annotated sites; supply them in the target's metadata."
            )

        # Split the annotated sites into two coupled clusters so the production
        # selector sees interacting groups rather than one undifferentiated block.
        midpoint = (len(targetable) + 1) // 2
        groups = [targetable[:midpoint], targetable[midpoint:]]
        clusters = [
            ResidueCluster(
                cluster_id=f"annotated_site_group_{i + 1}",
                positions=group,
                rationale="Annotated active/binding-site positions reported for this target.",
            )
            for i, group in enumerate(groups)
            if group
        ]

        chemistries = {str(p): list(self.NEUTRAL_CHEMISTRIES) for p in targetable}
        hypothesis = (
            f"Round {round_number}: combinatorial exploration of the {len(targetable)} "
            "annotated active-site positions, preserving annotated catalytic machinery."
        )
        return SearchStrategy(
            target_protein=metadata.protein_name,
            target_accession=metadata.accession,
            hypothesis=hypothesis,
            target_clusters=clusters,
            forbidden_positions=catalytic,
            # Depth tracks the annotated site count so a combinatorially complete
            # library is reachable, rather than being fixed to a literal.
            max_mutation_depth=min(len(targetable), 2 if round_number == 1 else 4),
            mechanistic_hypotheses=[
                BiophysicalHypothesis(
                    pocket_name="annotated_active_site",
                    target_residues=targetable,
                    mechanism="Combinatorial substitution across annotated active-site positions.",
                    structural_rationale=(
                        f"Positions {targetable} are the annotated substrate/cofactor contact "
                        f"sites for this entry; catalytic positions {catalytic} are withheld."
                    ),
                    suggested_chemistries=chemistries,
                    epistatic_coupling="Annotated co-located sites are treated as potentially coupled.",
                )
            ],
        )


class PanoptesSelector(Selector):
    """
    E/F. The real production planning + selection path.

    Production functions called, in production order:
        workers.rag_worker.fetch_protein_metadata      (step 1, cached, no network)
        _AnnotationDrivenPlanner.generate_strategy     (step 2, annotation-derived search space) [E]
        controller.orchestrator.Orchestrator._generate_strategy (step 2, live Muse planner)  [F]
        workers.math_worker.select_next_micro_batch    (step 3, DKL + GFlowNet + Pareto + k-DPP)
        controller.orchestrator.Orchestrator._generate_reflection (step 6b, feedback)        [F]

    Stages NOT exercised, and why:
        Step 4 (Gatekeeper spend authorisation) and step 5 (Modal ESM-2 scoring) are skipped.
        `Orchestrator.run_stream()` hard-wires step 5 to `workers.esm_worker.score_variants`,
        and the loop learns from those ESM pseudo-likelihoods via the ledger's evidence table.
        There is no supported interface for substituting an experimental assay as the outcome
        source, so running the full loop would optimise ESM scores rather than measured fitness
        and would consume none of the experimental budget. Pricing an ESM call this benchmark
        never makes would also fabricate spend. The planning, selection and feedback stages that
        decide WHICH variants are measured are exercised in full.

    Evidence handed to production contains only measurements already revealed to this method.
    """

    def __init__(
        self,
        seed: int,
        wt_sequence: str,
        metadata: Any,
        use_live_muse: bool,
        pdb_path: Optional[Path],
        call_budget: Optional[int] = None,
    ) -> None:
        super().__init__(seed)
        self.name = "panoptes_muse" if use_live_muse else "panoptes_default"
        self.wt_sequence = wt_sequence
        self.metadata = metadata
        self.use_live_muse = use_live_muse
        self.pdb_path = pdb_path
        self.call_budget = call_budget
        self.config = {
            "strategy_source": (
                "live Muse (Meta API)" if use_live_muse
                else "annotation-driven planner (metadata.binding_positions; no scripted per-target literals)"
            ),
            "selection": "workers.math_worker.select_next_micro_batch",
            "acquisition": "hybrid_thompson (production default)",
            "reflection_feedback": bool(use_live_muse),
            "stages_exercised": ["1_init", "2_strategy", "3_selection"] + (["6b_reflection"] if use_live_muse else []),
            "stages_skipped": ["4_authorization", "5_esm_execution"],
        }
        self._orchestrator = None
        self._mock_planner = None
        self.recorder = _LiveCallRecorder()
        self.live_calls = 0
        self.reflection_calls = 0
        self.last_batch_predictions: Dict[str, float] = {}

    # -- production handles -------------------------------------------------------------------
    def _planner(self):
        if self.use_live_muse:
            if self._orchestrator is None:
                from controller.orchestrator import Orchestrator, OrchestratorConfig

                self._orchestrator = Orchestrator(
                    OrchestratorConfig(
                        target_accession=PARENT_ACCESSION,
                        use_mock_llm=False,
                        enable_live_llm=True,
                        enable_reflection=True,
                        random_seed=self.seed,
                    ),
                    # The four saturated sites ARE this benchmark's ordered library:
                    # declaring them is experimental design, not a hint about which
                    # residues win. Without it the planner proposes distal allosteric
                    # residues that the four-site library cannot measure, and whole
                    # rounds are discarded as not_in_eligible_universe.
                    assayable_positions=ACTIVE_SITE_POSITIONS,
                )
            return self._orchestrator
        if self._mock_planner is None:
            self._mock_planner = _AnnotationDrivenPlanner()
        return self._mock_planner

    def _strategy(self, round_index: int, revealed: RevealedState):
        """Step 2, through whichever production planner this method is configured for."""
        top_evidence = [
            {"variant_sequence": self._sequence_for(v), "computed_score": score}
            for v, score in sorted(revealed.measured.items(), key=lambda kv: kv[1], reverse=True)[:5]
        ]
        planner = self._planner()
        if self.use_live_muse:
            if self.call_budget is not None and self.live_calls >= self.call_budget:
                raise BenchmarkError(
                    f"live model-call cap of {self.call_budget} reached (--max-model-calls)"
                )
            logging.getLogger("controller.orchestrator").addHandler(self.recorder)
            before = self.recorder.live_attempts
            try:
                strategy = planner._generate_strategy(
                    self.metadata, round_number=round_index + 1, top_evidence=top_evidence, current_hop=0,
                )
            finally:
                logging.getLogger("controller.orchestrator").removeHandler(self.recorder)
            if self.recorder.live_attempts > before:
                self.live_calls += 1
                self.model_calls += 1
            return strategy
        return planner.generate_strategy(
            self.metadata, round_number=round_index + 1, top_evidence=top_evidence, current_hop=0,
        )

    def _sequence_for(self, variant_id: str) -> str:
        from controller.gatekeeper import Variant

        return Variant(notation=variant_id).apply_to(self.wt_sequence)

    def _reflect(self, round_index: int, revealed_batch: Dict[str, float]) -> None:
        """Step 6b: production residual reflection, fed back exactly as run_stream() does."""
        from controller.prompts import ResidualOutlier

        outliers: List[ResidualOutlier] = []
        for variant_id, observed in revealed_batch.items():
            predicted = self.last_batch_predictions.get(variant_id, 0.0)
            outliers.append(ResidualOutlier(
                variant_notation=variant_id,
                predicted_score=round(float(predicted), 4),
                observed_score=round(float(observed), 4),
                residual=round(float(observed) - float(predicted), 4),
                mutation_details=variant_id.split(":"),
            ))
        outliers.sort(key=lambda o: o.residual, reverse=True)
        positive = [o for o in outliers if o.residual > 0.05][:3]
        negative = [o for o in reversed(outliers) if o.residual < -0.05][:3]
        if not positive and not negative:
            return
        if self.call_budget is not None and self.live_calls >= self.call_budget:
            self.notes.append("reflection skipped: live model-call cap reached")
            return
        orchestrator = self._planner()
        logging.getLogger("controller.orchestrator").addHandler(self.recorder)
        before = self.recorder.live_attempts
        try:
            reflection = orchestrator._generate_reflection(
                metadata=self.metadata,
                round_number=round_index + 1,
                positive_outliers=positive,
                negative_outliers=negative,
                pdb_id=getattr(self.metadata, "pdb_id", None),
            )
        finally:
            logging.getLogger("controller.orchestrator").removeHandler(self.recorder)
        # Production stores the diagnosis on the orchestrator; the next strategy merges it in.
        orchestrator.last_reflection = reflection
        self.reflection_calls += 1
        if self.recorder.live_attempts > before:
            self.live_calls += 1
            self.model_calls += 1

    # -- selection ----------------------------------------------------------------------------
    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        from workers.math_worker import select_next_micro_batch

        strategy = self._strategy(round_index, revealed)
        evidence = [
            {"variant_sequence": self._sequence_for(v), "computed_score": float(score)}
            for v, score in revealed.measured.items()
        ]
        # select_next_micro_batch draws its Thompson sample from numpy's global RNG; seeding it
        # here is a caller-side action that makes the production path reproducible.
        np.random.seed((self.seed * 7919 + round_index) % (2**31 - 1))
        batch = select_next_micro_batch(
            strategy=strategy,
            wt_sequence=self.wt_sequence,
            evidence=evidence,
            batch_size=n,
            catalytic_invariants=set(getattr(self.metadata, "catalytic_positions", []) or []),
            pdb_path=self.pdb_path,
            random_seed=(self.seed or 0) * 1000 + round_index,
        )
        self.last_batch_predictions = {}
        proposals: List[str] = []
        for variant, mean in zip(batch.variants, batch.predicted_means or [0.0] * len(batch.variants)):
            combo = _combo_from_notation(variant.notation)
            variant_id = _canonical_variant_id(combo) if combo else variant.notation
            self.last_batch_predictions[variant_id] = float(mean)
            proposals.append(variant_id)
        return proposals

    def observe(self, round_index: int, revealed_batch: Dict[str, float]) -> None:
        """Production's between-round feedback stage; only live-Muse runs implement one."""
        if self.use_live_muse and revealed_batch:
            try:
                self._reflect(round_index, revealed_batch)
            except Exception as exc:
                self.notes.append(f"reflection failed: {exc}")


# ==============================================================================================
# Protocol
# ==============================================================================================
@dataclass
class RoundMetrics:
    method: str
    seed: int
    round_index: int
    n_measurements: int
    best_fitness: float
    new_top1_hits: int
    cum_top1_excluding_initial: int
    cum_top1_including_initial: int
    found_top1: bool
    round_seconds: float
    cumulative_seconds: float
    model_calls: int
    token_usage: str
    requested: int
    accepted: int
    rejected: int
    status: str
    failure_reason: str


@dataclass
class TraceRow:
    method: str
    seed: int
    round_index: int
    rank: int
    variant_id: str
    status: str
    reason: str
    revealed_fitness: Optional[float]


@dataclass
class RunOutcome:
    method: str
    seed: int
    status: str                       # completed | failed | unavailable | unsupported | skipped
    failure_reason: str = ""
    rounds: List[RoundMetrics] = field(default_factory=list)
    trace: List[TraceRow] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    model_calls: int = 0
    token_usage: Optional[int] = None


def run_single(
    selector: Selector,
    evaluator: Evaluator,
    view: CandidateView,
    initial_ids: Sequence[str],
    seed: int,
) -> RunOutcome:
    """
    Runs one method for one seed under the shared protocol.

    The initial sample is identical across methods within a seed; everything after it diverges,
    because each method learns only from the outcomes of its own choices.
    """
    outcome = RunOutcome(method=selector.name, seed=seed, status="completed")
    revealed = RevealedState()
    started = time.perf_counter()

    # Round 0: the shared starting sample.
    initial = evaluator.reveal(initial_ids)
    revealed.measured.update(initial)
    initial_hits = sum(1 for v in initial_ids if evaluator.is_top1(v))
    for rank, variant_id in enumerate(initial_ids):
        outcome.trace.append(TraceRow(
            method=selector.name, seed=seed, round_index=0, rank=rank, variant_id=variant_id,
            status="initial", reason="", revealed_fitness=initial[variant_id],
        ))
    elapsed = time.perf_counter() - started
    outcome.rounds.append(RoundMetrics(
        method=selector.name, seed=seed, round_index=0,
        n_measurements=len(revealed.measured), best_fitness=revealed.best(),
        new_top1_hits=initial_hits, cum_top1_excluding_initial=0,
        cum_top1_including_initial=initial_hits, found_top1=initial_hits > 0,
        round_seconds=round(elapsed, 4), cumulative_seconds=round(elapsed, 4),
        model_calls=0, token_usage="unknown", requested=len(initial_ids),
        accepted=len(initial_ids), rejected=0, status="completed", failure_reason="",
    ))

    # Harness-style methods (ALDE) run their own internal loop and return every round's picks
    # at once. The outcomes are still revealed round by round below, and the metrics, rejection
    # rules and trace are identical to a per-round method's.
    campaign: Optional[List[List[str]]] = None
    campaign_runner = getattr(selector, "run_campaign", None)
    if callable(campaign_runner):
        attach = getattr(selector, "attach_truth", None)
        if callable(attach):
            attach(evaluator.export_truth_table(selector.name))
        try:
            campaign = campaign_runner(view, initial_ids, view.n_rounds, view.batch_size)
        except Exception as exc:
            outcome.status = "failed"
            outcome.failure_reason = f"{type(exc).__name__}: {exc}"
            logger.error("%s seed %d campaign failed: %s", selector.name, seed, exc)
            outcome.model_calls = selector.model_calls
            outcome.notes.extend(selector.notes)
            return outcome

    cumulative_hits_excluding_initial = 0
    for round_index in range(1, view.n_rounds + 1):
        round_started = time.perf_counter()
        try:
            if campaign is not None:
                proposals = campaign[round_index - 1] if round_index - 1 < len(campaign) else []
            else:
                proposals = selector.propose(round_index, view.batch_size, view, revealed)
        except Exception as exc:
            outcome.status = "failed"
            outcome.failure_reason = f"round {round_index}: {type(exc).__name__}: {exc}"
            logger.error("%s seed %d failed at round %d: %s", selector.name, seed, round_index, exc)
            break

        accepted: List[str] = []
        seen_this_batch: Set[str] = set()
        for rank, variant_id in enumerate(proposals):
            if len(accepted) >= view.batch_size:
                reason = "over_batch_size"
            elif variant_id in seen_this_batch:
                reason = "duplicate_within_batch"
            elif variant_id in revealed.measured:
                reason = "already_observed"
            elif variant_id not in view.index_of:
                reason = "not_in_eligible_universe"
            else:
                reason = ""
            if reason:
                outcome.trace.append(TraceRow(
                    method=selector.name, seed=seed, round_index=round_index, rank=rank,
                    variant_id=variant_id, status="rejected", reason=reason, revealed_fitness=None,
                ))
                continue
            seen_this_batch.add(variant_id)
            accepted.append(variant_id)

        # The batch is committed here; only now does the evaluator reveal its outcomes.
        batch_results = evaluator.reveal(accepted)
        revealed.measured.update(batch_results)
        new_hits = sum(1 for v in accepted if evaluator.is_top1(v))
        cumulative_hits_excluding_initial += new_hits
        for rank, variant_id in enumerate(accepted):
            outcome.trace.append(TraceRow(
                method=selector.name, seed=seed, round_index=round_index, rank=rank,
                variant_id=variant_id, status="accepted", reason="",
                revealed_fitness=batch_results[variant_id],
            ))

        observe = getattr(selector, "observe", None)
        if callable(observe):
            try:
                observe(round_index, batch_results)
            except Exception as exc:
                selector.notes.append(f"observe failed at round {round_index}: {exc}")

        round_seconds = time.perf_counter() - round_started
        rejected = len(proposals) - len(accepted)
        status = "completed" if len(accepted) == view.batch_size else "short_batch"
        outcome.rounds.append(RoundMetrics(
            method=selector.name, seed=seed, round_index=round_index,
            n_measurements=len(revealed.measured), best_fitness=revealed.best(),
            new_top1_hits=new_hits,
            cum_top1_excluding_initial=cumulative_hits_excluding_initial,
            cum_top1_including_initial=cumulative_hits_excluding_initial + initial_hits,
            found_top1=(cumulative_hits_excluding_initial + initial_hits) > 0,
            round_seconds=round(round_seconds, 4),
            cumulative_seconds=round(time.perf_counter() - started, 4),
            model_calls=selector.model_calls,
            token_usage="unknown" if selector.token_usage is None else str(selector.token_usage),
            requested=view.batch_size, accepted=len(accepted), rejected=rejected,
            status=status,
            failure_reason="" if status == "completed" else f"{rejected} selection(s) rejected",
        ))
        if status == "short_batch":
            logger.warning(
                "%s seed %d round %d: only %d/%d selections were usable.",
                selector.name, seed, round_index, len(accepted), view.batch_size,
            )

    outcome.model_calls = selector.model_calls
    outcome.token_usage = selector.token_usage
    outcome.notes.extend(selector.notes)
    return outcome


# ==============================================================================================
# Self-checks (run as part of --smoke-test)
# ==============================================================================================
class _ScriptedSelector(Selector):
    """A selector that returns a fixed list, used to prove the rejection rules actually bite."""

    name = "scripted"

    def __init__(self, seed: int, scripted: Sequence[str]) -> None:
        super().__init__(seed)
        self.scripted = list(scripted)

    def propose(self, round_index: int, n: int, view: CandidateView, revealed: RevealedState) -> List[str]:
        return list(self.scripted)


def run_self_checks() -> List[Tuple[str, bool, str]]:
    """
    Validates the benchmark's own mechanics: encoding round trips, evaluator isolation,
    selection rejection rules and budget accounting. Returns (check, passed, detail) rows.
    """
    results: List[Tuple[str, bool, str]] = []

    def check(name: str, condition: bool, detail: str = "") -> None:
        results.append((name, bool(condition), detail))

    # 1. Notation <-> combo round trip, including parent-residue slots.
    combo = "AIKG"
    notation = _canonical_variant_id(combo)
    check("notation round trip", _combo_from_notation(notation) == combo, f"{combo} -> {notation}")
    partial = _canonical_variant_id("AFVS")  # only position 183 differs from the parent
    check("partial variant maps back", _combo_from_notation(partial) == "AFVS", f"'{partial}'")
    check("off-library notation rejected", _combo_from_notation("S65T") is None, "position outside the library")
    check("parent-mismatch rejected", _combo_from_notation("A183V") is None, "wrong parent residue")

    universe = build_synthetic_universe(n_variants=600, seed=3)
    evaluator = Evaluator(universe)
    view = CandidateView(
        variant_ids=universe.variant_ids, one_hot=universe.one_hot, index_of=universe.index_of,
        combos=universe.combos, batch_size=5, n_rounds=1,
    )

    # 2. The selector-visible view carries no outcome information.
    check(
        "candidate view holds no fitness",
        not any("fitness" in field_name for field_name in vars(view)),
        f"fields: {sorted(vars(view))}",
    )

    # 3. Rejection rules: duplicate, already-observed, unknown, and over-budget selections.
    initial = universe.variant_ids[:5]
    fresh = universe.variant_ids[10:13]
    scripted = [fresh[0], fresh[0], initial[0], "V183Z", fresh[1], fresh[2]]
    outcome = run_single(_ScriptedSelector(0, scripted), evaluator, view, initial, seed=0)
    reasons = {row.reason for row in outcome.trace if row.status == "rejected"}
    accepted = [row.variant_id for row in outcome.trace if row.status == "accepted"]
    check("duplicate rejected", "duplicate_within_batch" in reasons, str(sorted(reasons)))
    check("already-observed rejected", "already_observed" in reasons, str(sorted(reasons)))
    check("unknown variant rejected", "not_in_eligible_universe" in reasons, str(sorted(reasons)))
    check("accepted only valid picks", accepted == fresh, f"accepted={accepted}")
    check("short batch is not backfilled", outcome.rounds[-1].status == "short_batch",
          f"accepted {outcome.rounds[-1].accepted}/{outcome.rounds[-1].requested}")

    # 4. Budget accounting: unique measurements == initial + accepted, never more.
    final = outcome.rounds[-1]
    check("budget counted from unique measurements",
          final.n_measurements == len(initial) + len(fresh),
          f"{final.n_measurements} == {len(initial)} + {len(fresh)}")

    # 5. A full-budget run consumes exactly the configured budget.
    view_full = CandidateView(
        variant_ids=universe.variant_ids, one_hot=universe.one_hot, index_of=universe.index_of,
        combos=universe.combos, batch_size=10, n_rounds=3,
    )
    full = run_single(RandomSelector(1), evaluator, view_full, universe.variant_ids[:10], seed=1)
    check("full budget consumed", full.rounds[-1].n_measurements == 10 + 3 * 10,
          f"{full.rounds[-1].n_measurements} measurements")
    check("top-1% threshold is evaluator-side",
          evaluator.top1_threshold > min(universe.fitness.values()),
          f"threshold {evaluator.top1_threshold:.4g}")
    return results


# ==============================================================================================
# Outputs
# ==============================================================================================
def _write_csv(path: Path, header: Sequence[str], rows: Iterable[Sequence[Any]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def write_trace_csv(path: Path, outcomes: Sequence[RunOutcome], stamp: str) -> None:
    _write_csv(
        path,
        ["provenance", "method", "seed", "round", "rank", "variant_id", "status", "reason", "revealed_fitness"],
        [
            [stamp, r.method, r.seed, r.round_index, r.rank, r.variant_id, r.status, r.reason,
             "" if r.revealed_fitness is None else f"{r.revealed_fitness:.9g}"]
            for outcome in outcomes for r in outcome.trace
        ],
    )


def write_round_metrics_csv(path: Path, outcomes: Sequence[RunOutcome], stamp: str) -> None:
    _write_csv(
        path,
        ["provenance", "method", "seed", "round", "n_unique_measurements", "best_measured_fitness",
         "new_top1_hits", "cum_top1_hits_excluding_initial", "cum_top1_hits_including_initial",
         "found_top1", "round_seconds", "cumulative_seconds", "model_calls", "token_usage",
         "requested", "accepted", "rejected", "status", "failure_reason"],
        [
            [stamp, m.method, m.seed, m.round_index, m.n_measurements, f"{m.best_fitness:.9g}",
             m.new_top1_hits, m.cum_top1_excluding_initial, m.cum_top1_including_initial,
             int(m.found_top1), m.round_seconds, m.cumulative_seconds, m.model_calls, m.token_usage,
             m.requested, m.accepted, m.rejected, m.status, m.failure_reason]
            for outcome in outcomes for m in outcome.rounds
        ],
    )


def summarise(outcomes: Sequence[RunOutcome], method_status: Dict[str, Dict[str, str]]) -> List[Dict[str, Any]]:
    """Aggregates across seeds. Unavailable and failed methods stay visible with zero successes."""
    by_method: Dict[str, List[RunOutcome]] = {}
    for outcome in outcomes:
        by_method.setdefault(outcome.method, []).append(outcome)

    rows: List[Dict[str, Any]] = []
    for method, status in method_status.items():
        runs = by_method.get(method, [])
        completed = [r for r in runs if r.status == "completed" and r.rounds]
        finals = [r.rounds[-1] for r in completed]
        best = [m.best_fitness for m in finals]
        hits_incl = [m.cum_top1_including_initial for m in finals]
        hits_excl = [m.cum_top1_excluding_initial for m in finals]
        measurements = [m.n_measurements for m in finals]
        runtimes = [m.cumulative_seconds for m in finals]

        def _mean(values: Sequence[float]) -> str:
            return f"{statistics.fmean(values):.6g}" if values else ""

        def _stdev(values: Sequence[float]) -> str:
            return f"{statistics.stdev(values):.6g}" if len(values) > 1 else ("0" if values else "")

        rows.append({
            "method": method,
            "status": status.get("status", "unknown"),
            "status_detail": status.get("detail", ""),
            "seeds_attempted": len(runs),
            "seeds_completed": len(completed),
            "completion_rate": f"{len(completed) / len(runs):.3f}" if runs else "0.000",
            "mean_measurements_used": _mean(measurements),
            "mean_best_fitness": _mean(best),
            "std_best_fitness": _stdev(best),
            "mean_top1_hits_including_initial": _mean(hits_incl),
            "std_top1_hits_including_initial": _stdev(hits_incl),
            "mean_top1_hits_excluding_initial": _mean(hits_excl),
            "std_top1_hits_excluding_initial": _stdev(hits_excl),
            "success_rate_found_top1": f"{sum(1 for m in finals if m.found_top1) / len(finals):.3f}" if finals else "",
            "mean_runtime_seconds": _mean(runtimes),
            "total_model_calls": sum(r.model_calls for r in runs),
            "token_usage": "unknown",
            "estimated_cost_usd": "unknown",
        })
    return rows


def write_summary_csv(path: Path, rows: Sequence[Dict[str, Any]], stamp: str) -> None:
    if not rows:
        return
    header = ["provenance"] + list(rows[0].keys())
    _write_csv(path, header, [[stamp] + [row[k] for k in rows[0]] for row in rows])


def write_plot(
    path: Path,
    outcomes: Sequence[RunOutcome],
    stamp: str,
    synthetic: bool,
    cache_dir: Path,
) -> Optional[str]:
    """Best measured fitness vs total measurements, mean with a +/-1 s.d. band across seeds."""
    os.environ["MPLCONFIGDIR"] = str(cache_dir)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        return f"plot skipped: matplotlib unavailable ({exc})"

    series: Dict[str, Dict[int, List[float]]] = {}
    for outcome in outcomes:
        if outcome.status != "completed":
            continue
        for metric in outcome.rounds:
            series.setdefault(outcome.method, {}).setdefault(metric.n_measurements, []).append(metric.best_fitness)
    if not series:
        return "plot skipped: no completed runs"

    figure, axes = plt.subplots(figsize=(9, 5.5))
    for method in sorted(series):
        points = sorted(series[method])
        means = np.array([statistics.fmean(series[method][p]) for p in points])
        devs = np.array([statistics.stdev(series[method][p]) if len(series[method][p]) > 1 else 0.0 for p in points])
        axes.plot(points, means, marker="o", label=method)
        axes.fill_between(points, means - devs, means + devs, alpha=0.15)

    axes.set_xlabel("Unique experimental measurements used")
    axes.set_ylabel("Best measured fitness found")
    title = "TrpB four-site landscape - best measured fitness vs budget"
    if synthetic:
        title = "SYNTHETIC landscape (NOT biological) - best value vs budget"
    axes.set_title(title)
    axes.legend()
    axes.grid(alpha=0.25)
    figure.text(0.01, 0.01, stamp, fontsize=7, color="#666666")
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return None


def write_report(
    path: Path,
    manifest: Dict[str, Any],
    summary_rows: Sequence[Dict[str, Any]],
    outcomes: Sequence[RunOutcome],
    stamp: str,
    plot_note: Optional[str],
) -> None:
    lines: List[str] = []
    add = lines.append
    add("PanoptesBio isolated comparison benchmark")
    add("=" * 78)
    add(f"Run id:        {manifest['run_id']}")
    add(f"Started (UTC): {manifest['started_utc']}")
    add(f"Provenance:    {stamp}")
    add(f"Dataset:       {manifest['dataset']['path']}")
    add(f"               sha256 {manifest['dataset']['sha256'][:16]}..., "
        f"{manifest['dataset']['unique_variants']:,} eligible variants")
    add("")
    add("Protocol")
    add("-" * 78)
    config = manifest["config"]
    add(f"  {config['seeds']} paired seeds; {config['initial_size']} shared starting variants, "
        f"then {config['rounds']} rounds x {config['batch_size']}.")
    add(f"  Total budget per method per seed: {config['total_budget']} unique measurements, including the start.")
    add(f"  Fitness direction: {FITNESS_DIRECTION}.")
    add(f"  Top-1% threshold: {manifest['top1']['threshold']:.6g} (99th percentile over the eligible set).")
    add(f"  Ties: every variant at or above the threshold counts as a hit - "
        f"{manifest['top1']['n_hits']:,} of {manifest['dataset']['unique_variants']:,} variants qualify.")
    add(f"  Global optimum in the eligible set: {manifest['top1']['global_best']:.6g}.")
    add("")
    add("Results")
    add("-" * 78)
    def _num(value: Any, places: int = 4) -> str:
        try:
            return f"{float(value):.{places}f}"
        except (TypeError, ValueError):
            return "-"

    add(f"  {'method':<17} {'status':<13} {'done':>6} {'best fit':>10} {'sd':>8} "
        f"{'top1':>7} {'sd':>7} {'secs':>8}")
    for row in summary_rows:
        add(f"  {row['method']:<17} {row['status']:<13} "
            f"{str(row['seeds_completed']) + '/' + str(row['seeds_attempted']):>6} "
            f"{_num(row['mean_best_fitness']):>10} {_num(row['std_best_fitness']):>8} "
            f"{_num(row['mean_top1_hits_including_initial'], 2):>7} "
            f"{_num(row['std_top1_hits_including_initial'], 2):>7} "
            f"{_num(row['mean_runtime_seconds'], 1):>8}")
    add("")
    add("  'best fit' is the best measured fitness found within the budget, averaged over seeds.")
    add("  'top1' counts top-1% variants found, INCLUDING those present in the shared starting sample")
    add("  (the per-round CSV reports the excluding-initial count separately).")
    add("")
    add("Method status")
    add("-" * 78)
    for row in summary_rows:
        add(f"  {row['method']}: {row['status']}")
        if row["status_detail"]:
            for chunk in str(row["status_detail"]).split("; "):
                add(f"      - {chunk}")
    add("")
    failures = [o for o in outcomes if o.status not in ("completed",)]
    if failures:
        add("Failed runs (not excluded from the denominators above)")
        add("-" * 78)
        for outcome in failures:
            add(f"  {outcome.method} seed {outcome.seed}: {outcome.status} - {outcome.failure_reason}")
        add("")
    notes = [(o.method, n) for o in outcomes for n in o.notes]
    if notes:
        add("Run notes")
        add("-" * 78)
        for method, note in notes:
            add(f"  {method}: {note}")
        add("")
    add("Limitations")
    add("-" * 78)
    for limitation in manifest["limitations"]:
        add(f"  - {limitation}")
    if plot_note:
        add(f"  - {plot_note}")
    add("")
    add("Cost")
    add("-" * 78)
    add("  No dollar figures are reported. Local baselines make no paid calls; the production")
    add("  adapters make no cloud calls in this benchmark; and the Meta API client in")
    add("  controller/orchestrator.py discards the provider's usage block, so token counts are")
    add("  genuinely unknown rather than estimated.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ==============================================================================================
# Integrity checks on the existing repository
# ==============================================================================================
def repository_state(ignore_prefixes: Sequence[str] = ()) -> Dict[str, Any]:
    """
    Records git status and a hash of tracked source/config files so drift is detectable.

    `ignore_prefixes` drops this benchmark's own outputs from the comparison: the run folder is
    created before the first snapshot and filled afterwards, which is expected, not drift.
    """
    state: Dict[str, Any] = {"git_available": False}

    def _is_ignored(line: str) -> bool:
        path = line[3:].strip().strip('"')
        return any(path.startswith(prefix) for prefix in ignore_prefixes)

    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT),
                              capture_output=True, text=True, timeout=30)
        status = subprocess.run(["git", "status", "--porcelain=v1", "--untracked-files=all"],
                                cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=60)
        if head.returncode == 0:
            state["git_available"] = True
            state["head"] = head.stdout.strip()
            state["status_lines"] = sorted(
                line for line in status.stdout.splitlines() if line.strip() and not _is_ignored(line)
            )
    except Exception as exc:
        state["error"] = str(exc)

    digest = hashlib.sha256()
    counted = 0
    for pattern in ("*.py", "*.txt", "*.json", "*.sql", "*.md"):
        for file_path in sorted(REPO_ROOT.glob(f"**/{pattern}")):
            parts = set(file_path.parts)
            if parts & {".git", "__pycache__", "node_modules", ".pytest_cache", "benchmark_comparison_results"}:
                continue
            if file_path.name == Path(__file__).name:
                continue
            relative = file_path.relative_to(REPO_ROOT).as_posix()
            if any(relative.startswith(prefix) for prefix in ignore_prefixes):
                continue
            try:
                digest.update(file_path.relative_to(REPO_ROOT).as_posix().encode())
                digest.update(file_path.read_bytes())
                counted += 1
            except OSError:
                continue
    state["source_tree_sha256"] = digest.hexdigest()
    state["files_hashed"] = counted
    return state


# ==============================================================================================
# Orchestration
# ==============================================================================================
ALL_METHODS = ("random", "ridge", "random_forest", "alde", "panoptes_default", "panoptes_muse")
DEFAULT_METHODS = ("random", "ridge", "random_forest", "panoptes_default")


def _load_panoptes_context(pdb: bool = True) -> Tuple[Any, str, Optional[Path]]:
    """Step 1 of the production loop: cached protein metadata, no network and no writes."""
    from workers.rag_worker import fetch_protein_metadata

    if not PARENT_METADATA_JSON.is_file():
        raise BenchmarkError(
            f"Parent metadata cache {PARENT_METADATA_JSON} is missing. Generate it with:\n"
            "  python scripts/download_trpb_benchmark.py"
        )
    metadata = fetch_protein_metadata(PARENT_ACCESSION, use_cache=True)
    sequence = metadata.sequence
    observed = "".join(sequence[p - 1] for p in ACTIVE_SITE_POSITIONS)
    if observed != PARENT_RESIDUES:
        raise BenchmarkError(
            f"Parent residues at {list(ACTIVE_SITE_POSITIONS)} are '{observed}', expected "
            f"'{PARENT_RESIDUES}'. Residue numbering does not match the dataset; refusing to run."
        )
    pdb_path = STRUCTURE_PDB if (pdb and STRUCTURE_PDB.is_file()) else None
    return metadata, sequence, pdb_path


def build_selector_factories(
    args: argparse.Namespace,
    methods: Sequence[str],
) -> Tuple[Dict[str, Callable[[int], Selector]], Dict[str, Dict[str, str]]]:
    """
    Resolves each requested method into a factory, or a visible unavailable/unsupported status.

    Optional integrations are resolved defensively: a failure here marks that one method and
    never prevents the local baselines from running.
    """
    factories: Dict[str, Callable[[int], Selector]] = {}
    status: Dict[str, Dict[str, str]] = {}

    for method in methods:
        if method == "random":
            factories[method] = lambda seed: RandomSelector(seed)
            status[method] = {"status": "available", "detail": "uniform random, no model"}
        elif method == "ridge":
            factories[method] = lambda seed: RidgeSelector(seed)
            status[method] = {"status": "available", "detail": "sklearn Ridge(alpha=1.0) on 80 one-hot features"}
        elif method == "random_forest":
            factories[method] = lambda seed: RandomForestSelector(seed)
            status[method] = {"status": "available", "detail": "sklearn RandomForestRegressor(n_estimators=300)"}
        elif method == "alde":
            checkout = Path(args.alde_path).expanduser().resolve() if args.alde_path else None
            if checkout is None and os.getenv("ALDE_PATH"):
                checkout = Path(os.environ["ALDE_PATH"]).expanduser().resolve()
            python_exe = Path(args.alde_python).expanduser().resolve() if args.alde_python else None
            if python_exe is None and checkout is not None:
                for candidate in (checkout.parent / ".venv-alde" / "Scripts" / "python.exe",
                                  checkout.parent / ".venv-alde" / "bin" / "python"):
                    if candidate.is_file():
                        python_exe = candidate
                        break
            probe = probe_alde(checkout, python_exe)
            if probe["available"]:
                workspace = (Path(args.alde_workspace).expanduser().resolve()
                             if args.alde_workspace else checkout.parent / "benchmark_runs")
                workspace.mkdir(parents=True, exist_ok=True)
                factories[method] = (
                    lambda seed, c=checkout, pe=python_exe, w=workspace, p=probe:
                    AldeSelector(seed, c, pe, w, p, timeout_seconds=args.alde_timeout,
                                 train_iter=args.alde_train_iter)
                )
                versions = probe.get("dependencies", {})
                status[method] = {
                    "status": "available",
                    "detail": (
                        f"checkout {probe['checkout_path']} @ {(probe['commit'] or 'unknown')[:12]}; "
                        f"onehot + DNN_ENSEMBLE + TS, run unmodified in {python_exe}; "
                        f"torch {versions.get('torch')}, botorch {versions.get('botorch')}, "
                        f"gpytorch {versions.get('gpytorch')}, numpy {versions.get('numpy')}; "
                        f"outputs under {workspace}"
                    ),
                }
            else:
                status[method] = {"status": "UNAVAILABLE", "detail": probe["reason"]}
            status[method]["probe"] = json.dumps(probe)
        elif method in ("panoptes_default", "panoptes_muse"):
            live = method == "panoptes_muse"
            if live and not args.enable_live_muse:
                status[method] = {
                    "status": "SKIPPED",
                    "detail": "live Muse comparison requires explicit opt-in: pass --enable-live-muse "
                              "(and set LLAMA_API_KEY or META_API_KEY)",
                }
                continue
            if live and not (os.getenv("LLAMA_API_KEY") or os.getenv("META_API_KEY")):
                status[method] = {
                    "status": "UNAVAILABLE",
                    "detail": "--enable-live-muse was passed but neither LLAMA_API_KEY nor META_API_KEY is "
                              "set, so no live call can be made. Refusing to run the mock planner and label "
                              "the result 'live'.",
                }
                continue
            try:
                metadata, sequence, pdb_path = _load_panoptes_context()
            except Exception as exc:
                status[method] = {"status": "UNSUPPORTED", "detail": f"{type(exc).__name__}: {exc}"}
                continue
            factories[method] = (
                lambda seed, m=metadata, s=sequence, p=pdb_path, l=live:
                PanoptesSelector(seed, s, m, l, p, call_budget=args.max_model_calls)
            )
            detail = (
                "production path: fetch_protein_metadata -> "
                + ("Orchestrator._generate_strategy (live Muse)" if live
                   else "_AnnotationDrivenPlanner.generate_strategy (search space from annotated sites)")
                + " -> select_next_micro_batch"
                + (" -> Orchestrator._generate_reflection" if live else "")
                + "; steps 4 (gatekeeper) and 5 (ESM execution) are not exercised - "
                  "run_stream() offers no interface to substitute the experimental assay for ESM scoring"
            )
            status[method] = {"status": "available", "detail": detail}
        else:
            status[method] = {"status": "UNSUPPORTED", "detail": f"unknown method '{method}'"}
    return factories, status


def execute(args: argparse.Namespace) -> int:
    started = datetime.now(timezone.utc)
    run_id = f"run_{started.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    output_root = Path(args.output_dir).expanduser().resolve()
    run_dir = output_root / run_id
    run_dir.mkdir(parents=True, exist_ok=False)  # never overwrites an existing run
    cache_dir = run_dir / "cache"
    cache_dir.mkdir()

    file_handler = logging.FileHandler(run_dir / "benchmark.log", encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
    logging.getLogger().addHandler(file_handler)
    logger.info("Run %s -> %s", run_id, run_dir)

    try:
        ignore_prefixes = (output_root.relative_to(REPO_ROOT).as_posix(),)
    except ValueError:
        ignore_prefixes = ()
    before_state = repository_state(ignore_prefixes)

    # ---- Dataset ---------------------------------------------------------------------------
    self_checks: List[Tuple[str, bool, str]] = []
    if args.smoke_test:
        self_checks = run_self_checks()
        failed = [name for name, passed, _ in self_checks if not passed]
        print("\nSelf-checks")
        print("-" * 78)
        for name, passed, detail in self_checks:
            print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
        if failed:
            raise BenchmarkError(f"self-checks failed: {failed}")
        universe = build_synthetic_universe(n_variants=args.synthetic_size, seed=7)
        logger.warning("SMOKE TEST: using a SYNTHETIC landscape. These numbers are not biology.")
    else:
        universe = load_universe(
            dataset_path=Path(args.dataset).expanduser().resolve(),
            provenance_mode=args.provenance_mode,
            imputed_ids_path=Path(args.imputed_ids).expanduser().resolve() if args.imputed_ids else None,
            max_variants=args.max_variants,
            subsample_seed=args.subsample_seed,
        )

    stamp = f"provenance={universe.provenance_status}"
    if universe.synthetic:
        stamp = "SYNTHETIC_NOT_BIOLOGICAL - mechanism validation only"
    logger.info("Eligible universe: %d variants (%s)", len(universe), universe.provenance_status)

    evaluator = Evaluator(universe)
    view = CandidateView(
        variant_ids=universe.variant_ids,
        one_hot=universe.one_hot,
        index_of=universe.index_of,
        combos=universe.combos,
        batch_size=args.batch_size,
        n_rounds=args.rounds,
    )

    total_budget = args.initial_size + args.rounds * args.batch_size
    if total_budget > len(universe):
        raise BenchmarkError(
            f"Budget of {total_budget} measurements exceeds the eligible universe of {len(universe)}."
        )

    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    unknown = [m for m in methods if m not in ALL_METHODS]
    if unknown:
        raise BenchmarkError(f"Unknown method(s) {unknown}. Choose from: {', '.join(ALL_METHODS)}")
    factories, method_status = build_selector_factories(args, methods)

    # ---- Run -------------------------------------------------------------------------------
    outcomes: List[RunOutcome] = []
    seeds = list(range(args.base_seed, args.base_seed + args.seeds))
    for seed in seeds:
        sampler = random.Random(seed)
        initial_ids = sampler.sample(universe.variant_ids, args.initial_size)  # shared across methods
        for method in methods:
            factory = factories.get(method)
            if factory is None:
                continue
            logger.info("seed %d | %s", seed, method)
            try:
                selector = factory(seed)
                outcome = run_single(selector, evaluator, view, initial_ids, seed)
            except Exception as exc:
                logger.error("%s seed %d could not start: %s", method, seed, exc)
                outcome = RunOutcome(method=method, seed=seed, status="failed",
                                     failure_reason=f"{type(exc).__name__}: {exc}")
            outcomes.append(outcome)
            if outcome.status == "completed" and outcome.rounds:
                final = outcome.rounds[-1]
                logger.info(
                    "  -> best %.6g | top1 %d (incl. start) | %d measurements | %.1fs",
                    final.best_fitness, final.cum_top1_including_initial,
                    final.n_measurements, final.cumulative_seconds,
                )

    for method in methods:
        if method in factories and any(o.method == method and o.status == "failed" for o in outcomes):
            method_status[method]["detail"] += "; one or more seeds failed (see report)"

    # ---- Outputs ---------------------------------------------------------------------------
    summary_rows = summarise(outcomes, method_status)
    write_trace_csv(run_dir / "selection_trace.csv", outcomes, stamp)
    write_round_metrics_csv(run_dir / "round_metrics.csv", outcomes, stamp)
    write_summary_csv(run_dir / "summary.csv", summary_rows, stamp)
    plot_note = write_plot(run_dir / "best_fitness_vs_measurements.png", outcomes, stamp,
                           universe.synthetic, cache_dir)
    if plot_note:
        logger.warning(plot_note)

    after_state = repository_state(ignore_prefixes)
    unchanged = (
        before_state.get("source_tree_sha256") == after_state.get("source_tree_sha256")
        and before_state.get("status_lines") == after_state.get("status_lines")
    )

    limitations = [
        "Methods diverge after the shared starting sample, so per-round comparisons are paired "
        "only at round 0.",
        "panoptes_default / panoptes_muse exercise production planning, selection and (live only) "
        "reflection. Gatekeeper authorisation and Modal ESM-2 execution are not exercised, because "
        "Orchestrator.run_stream() has no interface for substituting an experimental assay for ESM "
        "scoring in step 5.",
        "panoptes_default plans from the target's annotated active sites "
        "(metadata.binding_positions), not from a per-target script. It deliberately does NOT "
        "use workers.mock_workers.MockLLMWorker, whose TrpB branch hardcodes positions "
        "183/184/227/228 and suggests three of the four global-optimum residues - that would "
        "encode the answer this benchmark measures the discovery of.",
        "panoptes_muse is given this library's saturated positions as its assayable scope, "
        "which is enforced deterministically on the planner's output. Without it the planner "
        "proposed distal allosteric residues the four-site library cannot measure and whole "
        "rounds were discarded as not_in_eligible_universe. Check the accepted/requested "
        "columns in round_metrics.csv: a short batch still means proposals were rejected.",
        "Token usage is unknown: controller/orchestrator.py discards the provider's usage block.",
        "Top-1% membership is defined over the eligible set in this run; a different eligibility "
        "filter changes the threshold.",
    ]
    if universe.synthetic:
        limitations.insert(0, "SYNTHETIC RUN: values are generated, not measured. These results say "
                              "nothing about biological performance.")
    if universe.provenance_status == "UNVERIFIED_INCLUDES_IMPUTED":
        limitations.insert(0, f"Measurement provenance is UNVERIFIED: up to {PUBLISHED_IMPUTED_COUNT} "
                              "published KNN-imputed values may be included in the eligible set.")
    if not unchanged:
        limitations.insert(0, "WARNING: the repository source hash changed during this run - see manifest.")

    manifest = {
        "run_id": run_id,
        "benchmark_version": BENCHMARK_VERSION,
        "benchmark_file": str(Path(__file__).resolve()),
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "config": {
            "smoke_test": bool(args.smoke_test),
            "seeds": args.seeds,
            "base_seed": args.base_seed,
            "seed_values": seeds,
            "initial_size": args.initial_size,
            "rounds": args.rounds,
            "batch_size": args.batch_size,
            "total_budget": total_budget,
            "methods_requested": methods,
            "provenance_mode": args.provenance_mode,
            "enable_live_muse": bool(args.enable_live_muse),
            "max_model_calls": args.max_model_calls,
            "max_variants": args.max_variants,
        },
        "dataset": {
            **universe.fingerprint.to_dict(),
            "provenance_status": universe.provenance_status,
            "provenance_detail": universe.provenance_detail,
            "excluded_imputed": universe.excluded_imputed,
            "duplicates_dropped": universe.duplicates_dropped,
            "synthetic": universe.synthetic,
            "fitness_direction": FITNESS_DIRECTION,
            "active_site_positions": list(ACTIVE_SITE_POSITIONS),
            "parent_residues": PARENT_RESIDUES,
            "sequence_mapping": "variant id == position-ordered notation vs the Tm9D8* parent; "
                                "selections are matched to the library by their 4-residue combo",
            "published_split": {
                "measured": PUBLISHED_MEASURED_COUNT,
                "imputed_knn": PUBLISHED_IMPUTED_COUNT,
                "source": "bioRxiv 10.1101/2024.06.23.600144v1",
            },
        },
        "top1": {
            "threshold": evaluator.top1_threshold,
            "definition": "99th percentile of measured fitness over the eligible set",
            "tie_handling": "fitness >= threshold counts as a hit; ties all count",
            "n_hits": evaluator.n_top1,
            "global_best": evaluator.global_best,
            "visibility": "evaluator-only; never passed to a selector",
        },
        "methods": {
            method: {
                **method_status.get(method, {}),
                "config": (factories[method](0).config if method in factories else {}),
            }
            for method in methods
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "packages": _package_versions(),
        },
        "repository_state_before": before_state,
        "repository_state_after": after_state,
        "repository_unchanged": unchanged,
        "self_checks": [{"check": n, "passed": p, "detail": d} for n, p, d in self_checks],
        "limitations": limitations,
    }
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    write_report(run_dir / "report.txt", manifest, summary_rows, outcomes, stamp, plot_note)

    print()
    print((run_dir / "report.txt").read_text(encoding="utf-8"))
    print(f"Artifacts: {run_dir}")
    if not unchanged:
        print("WARNING: repository files changed during this run. Inspect manifest.json "
              "(repository_state_before / repository_state_after).", file=sys.stderr)
        return 2
    return 0


def _package_versions() -> Dict[str, Optional[str]]:
    versions: Dict[str, Optional[str]] = {}
    for name in ("sklearn", "scipy", "matplotlib", "torch", "gpytorch", "botorch", "pydantic"):
        try:
            versions[name] = getattr(__import__(name), "__version__", "unknown")
        except Exception:
            versions[name] = None
    return versions


# ==============================================================================================
# CLI
# ==============================================================================================
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="benchmark_panoptes_isolated.py",
        description=(
            "Isolated comparison benchmark for PanoptesBio on the TrpB four-site landscape. "
            "Writes only into a fresh folder under benchmark_comparison_results/ and never "
            "modifies the existing project."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # mechanism validation on synthetic data, a few seconds, no project data needed
  python benchmark_panoptes_isolated.py --smoke-test

  # real data, local baselines plus the production path (no paid calls)
  python benchmark_panoptes_isolated.py --imputed-ids data/trpb_imputed.txt

  # optional: the published ALDE implementation from a local checkout
  python benchmark_panoptes_isolated.py --methods alde --alde-path ~/src/ALDE

  # optional: live Muse planning (makes paid Meta API calls; capped)
  python benchmark_panoptes_isolated.py --methods panoptes_muse --enable-live-muse --max-model-calls 20
""",
    )
    parser.add_argument("--smoke-test", action="store_true",
                        help="Run a small bounded CPU test on a clearly-labelled SYNTHETIC landscape.")
    parser.add_argument("--dataset", default=str(DEFAULT_DATASET),
                        help="Path to the prepared TrpB CSV (default: %(default)s).")
    parser.add_argument("--methods", default=",".join(DEFAULT_METHODS),
                        help=f"Comma-separated subset of: {', '.join(ALL_METHODS)} (default: %(default)s).")
    parser.add_argument("--seeds", type=int, default=20, help="Number of paired seeds (default: %(default)s).")
    parser.add_argument("--base-seed", type=int, default=0, help="First seed value (default: %(default)s).")
    parser.add_argument("--initial-size", type=int, default=25,
                        help="Shared starting variants per seed (default: %(default)s).")
    parser.add_argument("--rounds", type=int, default=4,
                        help="Selection rounds after the start (default: %(default)s).")
    parser.add_argument("--batch-size", type=int, default=25,
                        help="New variants per round (default: %(default)s).")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_ROOT),
                        help="Parent folder for run folders (default: %(default)s).")
    parser.add_argument("--provenance-mode", choices=("require", "acknowledge-unflagged"), default="require",
                        help="'require' (default) refuses to run unless measured records can be "
                             "distinguished from the published KNN-imputed ones.")
    parser.add_argument("--imputed-ids", default=None,
                        help="File listing the imputed (non-measured) variants to exclude.")
    parser.add_argument("--max-variants", type=int, default=None,
                        help="Optionally subsample the candidate universe (seeded) for a faster run.")
    parser.add_argument("--subsample-seed", type=int, default=0, help="Seed for --max-variants.")
    parser.add_argument("--synthetic-size", type=int, default=0,
                        help="Variant count for the synthetic smoke-test landscape; 0 (the default) "
                             "covers the complete 20^4-1 combo space so the production adapter's "
                             "selections are in-universe.")
    parser.add_argument("--alde-path", default=None,
                        help="Path to an existing clone of https://github.com/jsunn-y/ALDE.")
    parser.add_argument("--alde-python", default=None,
                        help="Python executable of ALDE's own isolated environment. Defaults to "
                             "<alde-path>/../.venv-alde/Scripts/python.exe (or bin/python).")
    parser.add_argument("--alde-workspace", default=None,
                        help="Folder for ALDE's generated data, tensors and logs. Defaults to "
                             "<alde-path>/../benchmark_runs. Never inside this repository.")
    parser.add_argument("--alde-timeout", type=int, default=7200,
                        help="Per-seed timeout in seconds for the ALDE subprocess (default: %(default)s).")
    parser.add_argument("--alde-train-iter", type=int, default=300,
                        help="ALDE train_iter; 300 is the authors' default (default: %(default)s).")
    parser.add_argument("--enable-live-muse", action="store_true",
                        help="Explicitly permit live Meta/Muse API calls. Off by default; without it the "
                             "live comparison is reported as SKIPPED rather than silently mocked.")
    parser.add_argument("--max-model-calls", type=int, default=50,
                        help="Hard cap on live model calls across the whole run (default: %(default)s).")
    parser.add_argument("--verbose", action="store_true", help="Debug-level console logging.")
    return parser


def _load_repo_dotenv() -> None:
    """
    Loads REPO_ROOT/.env before any credential is read.

    `controller.orchestrator` calls load_dotenv() at import time, but it is imported
    lazily inside PanoptesSelector._planner(), which runs long after the live-Muse
    availability gate has already consulted os.environ. Without this, a key present
    in .env makes the gate report panoptes_muse as UNAVAILABLE and the live
    comparison is silently skipped. Deliberate and explicit rather than relying on
    another module's import side effect.
    """
    env_path = REPO_ROOT / ".env"
    if not env_path.is_file():
        return
    try:
        from dotenv import load_dotenv
    except ImportError:
        logger.debug("python-dotenv unavailable; relying on the ambient environment")
        return
    load_dotenv(env_path, override=False)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    _load_repo_dotenv()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    logging.getLogger("workers").setLevel(logging.WARNING)
    logging.getLogger("controller").setLevel(logging.WARNING)

    if args.smoke_test:
        # A bounded CPU run: small landscape, few seeds, short rounds, no paid calls.
        args.seeds = min(args.seeds, 2)
        args.initial_size = min(args.initial_size, 10)
        args.rounds = min(args.rounds, 2)
        args.batch_size = min(args.batch_size, 10)
        args.enable_live_muse = False

    for name, value in (("--seeds", args.seeds), ("--rounds", args.rounds),
                        ("--batch-size", args.batch_size), ("--initial-size", args.initial_size)):
        if value < 1:
            parser.error(f"{name} must be >= 1")

    try:
        return execute(args)
    except ProvenanceError as exc:
        print(f"\nPROVENANCE ERROR\n{'=' * 78}\n{exc}\n", file=sys.stderr)
        return 3
    except BenchmarkError as exc:
        print(f"\nERROR: {exc}\n", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

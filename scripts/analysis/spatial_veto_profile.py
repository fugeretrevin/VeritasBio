"""How much of each landscape can the spatial cohesion veto actually reach?

WHAT THIS MEASURES
------------------
`workers/spatial_filter.py::filter_spatially_viable_variants` is applied to the
engine's candidate pool inside `workers/math_worker.py::select_next_micro_batch`
(~line 2214). It rejects any multi-mutant whose mutated positions have a **mean
pairwise C-alpha distance > 16.0 A or a max pairwise distance > 24.0 A**
(`evaluate_mutation_cluster_cohesion`, ~line 349), and separately hard-rejects
variants that touch an annotated catalytic position or that collapse the buried
core (`check_steric_and_core_clash`, ~line 358).

**No CLI flag disables it.** `--ablate-components spatial` only sets
`spatial_map = None` (a feature cross-term) and `--ablate-components prior` only
removes the structural prior; the filter keeps running in both cases.

The cohesion veto is monotone in mutation depth: the more positions you mutate,
the likelier they span more than 16 A. So in `--production-candidate-source
enumerate` mode - where the engine is supposed to rank the *whole* untested
universe - the veto can silently replace the depth-4 generator ceiling with a
*spatial* ceiling. An experiment that concluded "enumerate did not help, so the
candidate pool was never the bottleneck" would then be measuring the filter, not
the pool.

This script answers, per library, over its FULL eligible universe:

  1. overall pass rate through `filter_spatially_viable_variants`;
  2. pass rate per mutation depth;
  3. **the fraction of that landscape's top-1% hits that survive the veto** -
     the decisive number, because a hit that cannot pass the filter is
     unreachable no matter how large the enumerated pool is;
  4. the mean / max pairwise C-alpha distance distribution per depth, so the
     mechanism is visible rather than inferred.

ONE NUANCE WORTH KNOWING
------------------------
The call site guards with `if filtered_candidates:` (math_worker.py:2225), so a
batch in which the filter rejects *everything* falls back to the unfiltered pool.
The veto therefore cannot empty the ballot - but as long as a single candidate
survives, every other candidate is removed. Over a whole enumerated universe
hundreds to thousands survive on every library measured here, so the fallback
never fires in `enumerate` mode and the pass rates below are the live ceiling.

HOW IT STAYS HONEST
-------------------
* Depths and the hit set come from `veritas.calibration.landscapes.load_landscape`,
  never from a new genotype parser (HANDOFF rule 2.1 - there are already five
  canonicalisations in this repo and one of them inverted a published result).
  Variant notation is built with `benchmark_panoptes_isolated._canonical_variant_id`.
* The pass/fail verdict is produced by calling the REAL
  `filter_spatially_viable_variants`, not a re-implementation. The only
  intervention is that `SpatialContactMap` is memoised for the duration of the
  run: the production function constructs a fresh one on **every call**
  (spatial_filter.py:704) rather than using the cached `get_spatial_contact_map`
  (:661), which would re-parse the PDB once per batch. The constructor arguments
  are recorded and asserted to be the ones the production call site uses
  (`pdb_path` only; `chain_id=None`, `reference_sequence=None`), so the cached
  map is the same object the production path would have built.
* A 100% pass rate has two completely different meanings - "the veto found
  nothing to reject" and "the PDB never loaded, so the filter returned the input
  unchanged" (spatial_filter.py:705-707). Every library therefore carries an
  explicit STRUCTURE verdict, and a third case is called out separately: the
  structure parses but the library's assayed positions are absent from it, which
  also makes the cohesion veto inert.
* The per-depth distance distribution is computed independently from a
  precomputed position-by-position C-alpha submatrix, and its cohesion verdict is
  cross-checked against the count of cohesion prunes the real filter reported. A
  disagreement is printed as a DISAGREEMENT line rather than hidden.

USAGE
-----
    .venv\\Scripts\\python.exe -X utf8 scripts/analysis/spatial_veto_profile.py
    .venv\\Scripts\\python.exe -X utf8 scripts/analysis/spatial_veto_profile.py --library creilov
    ... --max-variants 20000      # deterministic subsample, announced loudly
    ... --json out.json           # machine-readable copy of every table

Read-only: it loads prepared CSVs and PDBs and writes nothing unless --json is
given. It runs no campaign.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np  # noqa: E402

import benchmark_panoptes_isolated as bench  # noqa: E402
from workers import spatial_filter as sf  # noqa: E402
from veritas.calibration.landscapes import REGISTRY, load_landscape  # noqa: E402

# The thresholds under test, read from the module rather than retyped, so this
# report cannot drift from the code if someone retunes them.
MEAN_THRESHOLD_A = 16.0
MAX_THRESHOLD_A = 24.0

#: Universe size above which a deterministic subsample is taken by default.
#: No current library reaches it; it exists so the script stays runnable if one does.
DEFAULT_SAMPLE_TRIGGER = 400_000
DEFAULT_SAMPLE_SIZE = 20_000

COHESION_PREFIX = "Spatially uncoupled"
CATALYTIC_MARKERS = ("Catalytic Invariant Violation", "Catalytic Core Invariant Violation")

#: HANDOFF rule 2.2: reproduce both of these before trusting any new number out of
#: genotype code. {landscape: (hits, variants) in the shallow shell 1 <= depth <= 2}.
SHALLOW_SHELL_ANCHORS: Dict[str, Tuple[int, int]] = {"trpb": (168, 2_242), "gb1": (72, 2_167)}


def check_shallow_shell_anchors() -> bool:
    """Prints the two rule-2.2 anchors and whether they reproduce. Never skipped."""
    print()
    print("GENOTYPE VALIDATION (HANDOFF rule 2.2 - shallow shell 1 <= depth <= 2)")
    ok = True
    for key, (want_hits, want_n) in SHALLOW_SHELL_ANCHORS.items():
        try:
            shell = load_landscape(key).shell(2)
        except Exception as exc:
            print(f"  {key:<6} COULD NOT CHECK: {type(exc).__name__}: {exc}")
            ok = False
            continue
        got = (shell.n_hits, shell.n_variants)
        good = got == (want_hits, want_n)
        ok = ok and good
        print(f"  {key:<6} {shell.n_hits:>5,} / {shell.n_variants:>6,} = "
              f"{100.0 * shell.rate:7.4f}%   "
              f"{'PASS' if good else f'FAIL (expected {want_hits}/{want_n})'}")
    if not ok:
        print("  *** The depth/hit pipeline does not reproduce its anchors. Nothing below "
              "is trustworthy. ***")
    return ok


# ==============================================================================
# Structure health
# ==============================================================================
@dataclass
class StructureReport:
    """Everything needed to tell 'no veto' apart from 'no structure'."""

    pdb_path: Optional[Path]
    file_exists: bool
    file_bytes: Optional[int]
    n_ca: int
    chain_id: Optional[str]
    n_cofactor_atoms: int
    declared_positions: Tuple[int, ...]
    mapped_positions: Tuple[int, ...]
    missing_positions: Tuple[int, ...]
    identity_mismatches: Tuple[Tuple[int, str, str], ...]
    verdict: str
    detail: str

    @property
    def veto_active(self) -> bool:
        return self.verdict == "LIVE"

    def to_json(self) -> Dict[str, Any]:
        return {
            "pdb_path": str(self.pdb_path) if self.pdb_path else None,
            "file_exists": self.file_exists,
            "file_bytes": self.file_bytes,
            "n_ca_residues": self.n_ca,
            "chain_id": self.chain_id,
            "n_cofactor_atoms": self.n_cofactor_atoms,
            "n_declared_positions": len(self.declared_positions),
            "n_mapped_positions": len(self.mapped_positions),
            "missing_positions": list(self.missing_positions),
            "identity_mismatches": [
                {"position": p, "structure": s, "library_parent": l}
                for p, s, l in self.identity_mismatches
            ],
            "verdict": self.verdict,
            "detail": self.detail,
        }


_MAP_CACHE: Dict[Tuple[Optional[str], Optional[str], Optional[str]], sf.SpatialContactMap] = {}
_CTOR_CALLS: List[Tuple[Any, Optional[str], Optional[str]]] = []
_REAL_CONTACT_MAP = sf.SpatialContactMap


def _cached_contact_map(pdb_path=None, chain_id=None, reference_sequence=None):
    """Memoised stand-in for `SpatialContactMap`, installed into spatial_filter.

    `filter_spatially_viable_variants` builds a fresh map on every call, which
    would re-parse (and, for a missing file, re-download) the PDB once per batch.
    Every construction argument is recorded in `_CTOR_CALLS` so the caller can
    assert the production argument shape was never deviated from.
    """
    _CTOR_CALLS.append((pdb_path, chain_id, reference_sequence))
    key = (str(pdb_path) if pdb_path is not None else None, chain_id, reference_sequence)
    if key not in _MAP_CACHE:
        _MAP_CACHE[key] = _REAL_CONTACT_MAP(
            pdb_path, chain_id=chain_id, reference_sequence=reference_sequence
        )
    return _MAP_CACHE[key]


def inspect_structure(library: "bench.LibrarySpec") -> StructureReport:
    """Loads the library's structure exactly as production would and grades it."""
    pdb_path = library.structure_pdb
    exists = bool(pdb_path and Path(pdb_path).is_file())
    size = Path(pdb_path).stat().st_size if exists else None
    declared = tuple(library.positions)

    if not exists:
        return StructureReport(
            pdb_path=pdb_path, file_exists=False, file_bytes=None, n_ca=0, chain_id=None,
            n_cofactor_atoms=0, declared_positions=declared, mapped_positions=(),
            missing_positions=declared, identity_mismatches=(),
            verdict="NO_FILE",
            detail=f"structure_pdb not on disk: {pdb_path}",
        )

    # Production passes pdb_path ONLY (math_worker.py:2218-2222): no chain_id and no
    # reference_sequence, so the residue-identity guard at spatial_filter.py:221 never
    # fires in a benchmark run. Mirrored here deliberately.
    cmap = _cached_contact_map(pdb_path)
    mapped = tuple(p for p in declared if p in cmap.ca_coordinates)
    missing = tuple(p for p in declared if p not in cmap.ca_coordinates)
    mismatches: List[Tuple[int, str, str]] = []
    for slot, pos in enumerate(declared):
        entry = cmap.ca_coordinates.get(pos)
        if entry is None:
            continue
        one = sf.THREE_TO_ONE.get(entry[0], "?")
        want = library.parent_residues[slot]
        if one != want:
            mismatches.append((pos, one, want))

    n_cof = sum(len(v) for v in cmap.cofactor_coordinates.values())

    if not cmap.ca_coordinates:
        verdict, detail = "PARSE_FAILED", (
            "file present but zero C-alpha coordinates parsed; "
            "filter_spatially_viable_variants returns its input UNFILTERED "
            "(spatial_filter.py:705-707), so any pass rate below is vacuous"
        )
    elif len(mapped) < 2:
        verdict, detail = "POSITIONS_UNMAPPED", (
            f"structure parsed ({len(cmap.ca_coordinates)} C-alphas) but only "
            f"{len(mapped)} of {len(declared)} assayed positions exist in it, so no "
            f"pairwise distance is measurable and the cohesion veto is inert"
        )
    else:
        verdict = "LIVE"
        detail = f"{len(mapped)}/{len(declared)} assayed positions mapped"
        if missing:
            detail += (
                f"; {len(missing)} missing {list(missing)} - pairs touching those are "
                f"dropped from the cohesion mean, which WEAKENS the veto"
            )
    return StructureReport(
        pdb_path=Path(pdb_path), file_exists=True, file_bytes=size,
        n_ca=len(cmap.ca_coordinates), chain_id=cmap.chain_id, n_cofactor_atoms=n_cof,
        declared_positions=declared, mapped_positions=mapped, missing_positions=missing,
        identity_mismatches=tuple(mismatches), verdict=verdict, detail=detail,
    )


# ==============================================================================
# Distance mechanism (independent of the filter, cross-checked against it)
# ==============================================================================
class SlotGeometry:
    """Mean / max pairwise C-alpha distance for any subset of library slots.

    Built from a slot-by-slot submatrix of the contact map so the mechanism table
    is computed without a second pass through the filter. A pair whose residue is
    absent from the structure is NaN and is dropped from the mean, which is what
    `evaluate_mutation_cluster_cohesion` does with a `None` distance.
    """

    def __init__(self, cmap: sf.SpatialContactMap, positions: Sequence[int]) -> None:
        n = len(positions)
        self.n = n
        self.d = np.full((n, n), np.nan, dtype=np.float64)
        for i in range(n):
            for j in range(i + 1, n):
                dist = cmap.get_distance(positions[i], positions[j])
                if dist is not None:
                    self.d[i, j] = self.d[j, i] = float(dist)
        self._cache: Dict[Tuple[int, ...], Tuple[Optional[float], Optional[float]]] = {}

    def stats(self, slots: Tuple[int, ...]) -> Tuple[Optional[float], Optional[float]]:
        """(mean, max) over measurable pairs, or (None, None) when there are none."""
        if len(slots) <= 1:
            return None, None
        hit = self._cache.get(slots)
        if hit is not None:
            return hit
        idx = np.asarray(slots, dtype=np.intp)
        sub = self.d[np.ix_(idx, idx)]
        vals = sub[np.triu_indices(len(slots), k=1)]
        vals = vals[~np.isnan(vals)]
        out: Tuple[Optional[float], Optional[float]]
        if vals.size == 0:
            out = (None, None)
        else:
            out = (float(vals.mean()), float(vals.max()))
        self._cache[slots] = out
        return out


# ==============================================================================
# Per-library measurement
# ==============================================================================
@dataclass
class DepthRow:
    depth: int
    n: int = 0
    n_pass: int = 0
    n_hits: int = 0
    n_hits_pass: int = 0
    n_cohesion_pruned: int = 0
    n_catalytic_pruned: int = 0
    n_core_pruned: int = 0
    n_unmeasurable: int = 0
    mean_d: List[float] = field(default_factory=list)
    max_d: List[float] = field(default_factory=list)
    n_mean_over: int = 0
    n_max_over: int = 0

    @property
    def pass_pct(self) -> float:
        return 100.0 * self.n_pass / self.n if self.n else float("nan")

    @property
    def hit_pass_pct(self) -> float:
        return 100.0 * self.n_hits_pass / self.n_hits if self.n_hits else float("nan")

    def to_json(self) -> Dict[str, Any]:
        q = _quantiles(self.mean_d)
        r = _quantiles(self.max_d)
        return {
            "depth": self.depth, "n": self.n, "n_pass": self.n_pass,
            "pass_pct": _r(self.pass_pct), "n_hits": self.n_hits,
            "n_hits_pass": self.n_hits_pass, "hit_pass_pct": _r(self.hit_pass_pct),
            "pruned_cohesion": self.n_cohesion_pruned,
            "pruned_catalytic": self.n_catalytic_pruned,
            "pruned_core_or_steric": self.n_core_pruned,
            "n_distance_unmeasurable": self.n_unmeasurable,
            "mean_pairwise_A": q, "max_pairwise_A": r,
            "pct_mean_over_16": _r(100.0 * self.n_mean_over / self.n if self.n else float("nan")),
            "pct_max_over_24": _r(100.0 * self.n_max_over / self.n if self.n else float("nan")),
        }


def _quantiles(values: Sequence[float]) -> Dict[str, Optional[float]]:
    if not values:
        return {"p50": None, "p90": None, "max": None}
    arr = np.asarray(values, dtype=np.float64)
    return {
        "p50": _r(float(np.percentile(arr, 50))),
        "p90": _r(float(np.percentile(arr, 90))),
        "max": _r(float(arr.max())),
    }


def _r(x: Optional[float]) -> Optional[float]:
    if x is None:
        return None
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return round(float(x), 2)


@dataclass
class LibraryResult:
    key: str
    label: str
    structure: StructureReport
    catalytic_positions: List[int]
    catalytic_in_library: List[int]
    catalytic_note: str
    n_universe: int
    n_measured: int
    sampled: bool
    threshold: float
    threshold_source: str
    n_hits: int
    n_pass: int
    n_hits_pass: int
    rows: List[DepthRow]
    reason_counts: Dict[str, int]
    reason_examples: Dict[str, str]
    cross_check: str
    seconds: float

    @property
    def pass_pct(self) -> float:
        return 100.0 * self.n_pass / self.n_measured if self.n_measured else float("nan")

    @property
    def hit_pass_pct(self) -> float:
        return 100.0 * self.n_hits_pass / self.n_hits if self.n_hits else float("nan")

    def first_biting_depth(self) -> Optional[int]:
        """Shallowest depth at which more than 5% of variants are vetoed, any reason."""
        for row in self.rows:
            if row.n and (100.0 - row.pass_pct) > 5.0:
                return row.depth
        return None

    def first_cohesion_depth(self) -> Optional[int]:
        """Shallowest depth at which the COHESION veto alone removes more than 5%.

        Reported separately because cohesion is the depth-monotone veto; the core
        and catalytic checks can fire at depth 1 and say nothing about a ceiling.
        """
        for row in self.rows:
            if row.n and 100.0 * row.n_cohesion_pruned / row.n > 5.0:
                return row.depth
        return None

    def to_json(self) -> Dict[str, Any]:
        return {
            "library": self.key, "label": self.label,
            "structure": self.structure.to_json(),
            "catalytic_positions": self.catalytic_positions,
            "catalytic_positions_in_library": self.catalytic_in_library,
            "catalytic_note": self.catalytic_note,
            "n_universe": self.n_universe, "n_measured": self.n_measured,
            "sampled": self.sampled,
            "threshold": self.threshold, "threshold_source": self.threshold_source,
            "n_pass": self.n_pass, "pass_pct": _r(self.pass_pct),
            "n_hits": self.n_hits, "n_hits_pass": self.n_hits_pass,
            "hit_pass_pct": _r(self.hit_pass_pct),
            "first_biting_depth": self.first_biting_depth(),
            "first_cohesion_biting_depth": self.first_cohesion_depth(),
            "by_depth": [r.to_json() for r in self.rows],
            "prune_reasons": self.reason_counts,
            "prune_reason_examples": self.reason_examples,
            "cross_check": self.cross_check,
            "seconds": round(self.seconds, 1),
        }


def _classify(reason: str) -> str:
    if reason.startswith(COHESION_PREFIX):
        return "cohesion"
    if any(m in reason for m in CATALYTIC_MARKERS):
        return "catalytic"
    if reason.startswith("Cofactor"):
        return "cofactor"
    if reason.startswith("Steric clash"):
        return "steric"
    return "buried_core"


def _resolve_catalytic(library: "bench.LibrarySpec") -> Tuple[List[int], str]:
    """`metadata.catalytic_positions`, exactly as the benchmark resolves it."""
    try:
        metadata = bench.resolve_library_metadata(library)
    except Exception as exc:  # pragma: no cover - reported, never guessed
        return [], f"UNDETERMINED: resolve_library_metadata raised {type(exc).__name__}: {exc}"
    positions = sorted({int(p) for p in (getattr(metadata, "catalytic_positions", None) or [])})
    if not positions:
        if not library.accession:
            return [], ("none - this library declares no accession, so "
                        "resolve_library_metadata returns catalytic_positions=() "
                        "(benchmark_panoptes_isolated.py:5651)")
        return [], f"none annotated for accession {library.accession}"
    return positions, f"from accession {library.accession}"


def profile_library(
    key: str,
    *,
    max_variants: int,
    sample_trigger: int,
    sample_size: int,
    seed: int,
) -> LibraryResult:
    started = time.monotonic()
    library = bench.LIBRARIES[key]
    structure = inspect_structure(library)
    catalytic, catalytic_note = _resolve_catalytic(library)
    catalytic_in_library = sorted(set(catalytic) & set(library.positions))

    landscape = load_landscape(key)
    combos = landscape.combos
    depths = landscape.depths
    hit_mask = landscape.hit_mask
    n_universe = len(combos)

    # Deterministic subsample only if the universe is impractically large.
    cap = max_variants if max_variants > 0 else (
        sample_size if n_universe > sample_trigger else 0
    )
    if cap and cap < n_universe:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(n_universe, size=cap, replace=False))
        sampled = True
    else:
        idx = np.arange(n_universe)
        sampled = False

    notations = [bench._canonical_variant_id(combos[i], library) for i in range(len(combos))]

    # The real filter, on the real code path. catalytic_invariants is exactly what
    # the benchmark hands select_next_micro_batch (benchmark_panoptes_isolated.py:3429),
    # which forwards it unchanged (math_worker.py:2221).
    subset = [notations[i] for i in idx]
    approved, pruned = sf.filter_spatially_viable_variants(
        subset,
        pdb_path=library.structure_pdb,
        catalytic_invariants=set(catalytic) or None,
    )
    approved_set = set(approved)
    prune_reason: Dict[str, str] = dict(pruned)

    # Geometry, independent of the filter.
    geom = SlotGeometry(_cached_contact_map(library.structure_pdb), library.positions) \
        if structure.file_exists else None
    parent = library.parent_residues

    rows: Dict[int, DepthRow] = {}
    reason_counts: Dict[str, int] = {}
    reason_examples: Dict[str, str] = {}
    n_pass = n_hits = n_hits_pass = 0
    cohesion_from_geometry = 0

    for i in idx:
        combo = combos[i]
        d = int(depths[i])
        row = rows.setdefault(d, DepthRow(depth=d))
        row.n += 1
        is_hit = bool(hit_mask[i])
        if is_hit:
            n_hits += 1
            row.n_hits += 1

        passed = notations[i] in approved_set
        if passed:
            n_pass += 1
            row.n_pass += 1
            if is_hit:
                n_hits_pass += 1
                row.n_hits_pass += 1
        else:
            reason = prune_reason.get(notations[i], "<unrecorded>")
            bucket = _classify(reason)
            reason_counts[bucket] = reason_counts.get(bucket, 0) + 1
            reason_examples.setdefault(bucket, f"{notations[i]} -> {reason}")
            if bucket == "cohesion":
                row.n_cohesion_pruned += 1
            elif bucket == "catalytic":
                row.n_catalytic_pruned += 1
            else:
                row.n_core_pruned += 1

        if geom is not None and d >= 2:
            # Depth 1 has no pair by definition and is never cohesion-tested
            # (spatial_filter.py:731), so it is excluded from the distance columns
            # rather than counted as an unmeasurable pair.
            slots = tuple(s for s, ch in enumerate(combo) if ch != parent[s])
            mean_d, max_d = geom.stats(slots)
            if mean_d is None:
                # Every pair touches a position absent from the structure, so the
                # cohesion test cannot fire and the variant passes it vacuously.
                row.n_unmeasurable += 1
            else:
                row.mean_d.append(mean_d)
                row.max_d.append(max_d or 0.0)
                over_mean = mean_d > MEAN_THRESHOLD_A
                over_max = (max_d or 0.0) > MAX_THRESHOLD_A
                if over_mean:
                    row.n_mean_over += 1
                if over_max:
                    row.n_max_over += 1
                if over_mean or over_max:
                    cohesion_from_geometry += 1

    observed_cohesion = reason_counts.get("cohesion", 0)
    if geom is None:
        cross_check = "not run: no structure on disk"
    elif observed_cohesion == cohesion_from_geometry:
        cross_check = (
            f"OK - independent geometry predicts {cohesion_from_geometry} cohesion "
            f"rejections, the filter reported {observed_cohesion}"
        )
    else:
        cross_check = (
            f"DISAGREEMENT - independent geometry predicts {cohesion_from_geometry} "
            f"cohesion rejections, the filter reported {observed_cohesion}. Some "
            f"variants are rejected by an earlier check (catalytic / core) before the "
            f"cohesion test is reached, which accounts for a filter count that is LOWER."
        )

    return LibraryResult(
        key=key, label=library.label, structure=structure,
        catalytic_positions=catalytic, catalytic_in_library=catalytic_in_library,
        catalytic_note=catalytic_note,
        n_universe=n_universe, n_measured=len(idx), sampled=sampled,
        threshold=landscape.threshold, threshold_source=landscape.threshold_source,
        n_hits=n_hits, n_pass=n_pass, n_hits_pass=n_hits_pass,
        rows=[rows[d] for d in sorted(rows)],
        reason_counts=reason_counts, reason_examples=reason_examples,
        cross_check=cross_check, seconds=time.monotonic() - started,
    )


# ==============================================================================
# Reporting
# ==============================================================================
RULE = "=" * 100
THIN = "-" * 100


def _pct(x: float) -> str:
    return "     -" if not math.isfinite(x) else f"{x:6.2f}"


def _num(x: Optional[float]) -> str:
    return "    -" if x is None else f"{x:5.1f}"


def print_library(result: LibraryResult) -> None:
    s = result.structure
    print()
    print(RULE)
    print(f"LIBRARY  {result.key}   {result.label}")
    print(RULE)

    print("STRUCTURE")
    print(f"  pdb                 : {s.pdb_path}")
    print(f"  file                : {'present' if s.file_exists else 'MISSING'}"
          + (f", {s.file_bytes:,} bytes" if s.file_bytes else ""))
    print(f"  parsed              : {s.n_ca:,} C-alpha residues, chain "
          f"{s.chain_id or '<blank>'}, {s.n_cofactor_atoms} cofactor atoms")
    print(f"  assayed positions   : {len(s.mapped_positions)}/{len(s.declared_positions)} "
          f"present in the structure"
          + (f"; MISSING {list(s.missing_positions)}" if s.missing_positions else ""))
    if s.identity_mismatches:
        shown = ", ".join(f"{p}:{st}!={lb}" for p, st, lb in s.identity_mismatches[:6])
        print(f"  residue identity    : {len(s.identity_mismatches)} mismatch(es) vs the "
              f"library parent ({shown})")
        print("                        (production passes no reference_sequence, so this is "
              "NOT checked at runtime - distances are taken from whatever residue is there)")
    else:
        print("  residue identity    : all mapped positions match the library parent residue")
    print(f"  VETO STATUS         : {s.verdict} - {s.detail}")

    print(f"  catalytic_positions : {result.catalytic_positions or '[]'} "
          f"({result.catalytic_note})")
    print(f"    intersect library : {result.catalytic_in_library or 'NONE'}"
          + ("  <-- second, independent hard veto is ACTIVE"
             if result.catalytic_in_library else "  (second hard veto cannot fire)"))

    print()
    scope = (f"deterministic subsample of {result.n_measured:,}"
             if result.sampled else "FULL eligible universe")
    print(f"UNIVERSE              : {result.n_universe:,} eligible variants ({scope})")
    print(f"  top-1% threshold    : {result.threshold:.6g}  [{result.threshold_source}]")
    print(f"  hits in scope       : {result.n_hits:,}")
    print(f"  OVERALL PASS RATE   : {result.n_pass:,}/{result.n_measured:,} "
          f"= {_pct(result.pass_pct)}%")
    print(f"  HITS THAT CAN PASS  : {result.n_hits_pass:,}/{result.n_hits:,} "
          f"= {_pct(result.hit_pass_pct)}%   <<< DECISIVE NUMBER")

    print()
    print("PASS RATE BY MUTATION DEPTH")
    print(f"  {'d':>2} {'n':>9} {'pass':>9} {'pass%':>7} | "
          f"{'hits':>6} {'hitpass':>8} {'hitpass%':>9} | "
          f"{'meanD50':>7} {'meanD90':>7} {'maxD50':>7} {'maxD90':>7} | "
          f"{'%>16A':>7} {'%>24A':>7} | {'nodist':>7}")
    print("  " + THIN)
    for row in result.rows:
        mq, rq = _quantiles(row.mean_d), _quantiles(row.max_d)
        print(f"  {row.depth:>2} {row.n:>9,} {row.n_pass:>9,} {_pct(row.pass_pct):>7} | "
              f"{row.n_hits:>6,} {row.n_hits_pass:>8,} {_pct(row.hit_pass_pct):>9} | "
              f"{_num(mq['p50']):>7} {_num(mq['p90']):>7} "
              f"{_num(rq['p50']):>7} {_num(rq['p90']):>7} | "
              f"{_pct(100.0 * row.n_mean_over / row.n if row.n else float('nan')):>7} "
              f"{_pct(100.0 * row.n_max_over / row.n if row.n else float('nan')):>7} | "
              f"{row.n_unmeasurable:>7,}")

    if result.reason_counts:
        print()
        print("WHY VARIANTS WERE REJECTED")
        total = sum(result.reason_counts.values())
        for bucket, count in sorted(result.reason_counts.items(), key=lambda kv: -kv[1]):
            print(f"  {bucket:<14} {count:>9,}  ({100.0 * count / result.n_measured:5.2f}% "
                  f"of scope, {100.0 * count / total:5.1f}% of rejections)")
            print(f"                 e.g. {result.reason_examples.get(bucket, '')[:120]}")

    print()
    print(f"CROSS-CHECK           : {result.cross_check}")
    bite = result.first_biting_depth()
    coh = result.first_cohesion_depth()
    print(f"FIRST BITING DEPTH    : "
          + (f"depth {bite} (first depth losing >5% of its variants, any reason)"
             if bite is not None
             else "none - the filter never removes >5% at any depth"))
    print(f"  cohesion veto alone : "
          + (f"depth {coh} (first depth where the >16A/>24A rule alone removes >5%)"
             if coh is not None
             else "none - the cohesion rule never removes >5% at any depth"))
    print(f"(measured in {result.seconds:.1f}s)")


def print_summary(results: List[LibraryResult]) -> None:
    print()
    print(RULE)
    print("SUMMARY - can an enumerated pool reach the hits at all?")
    print(RULE)
    print(f"  {'library':<14} {'veto':<18} {'universe':>10} {'pass%':>7} "
          f"{'hits':>7} {'hits pass':>10} {'hit pass%':>10} {'any':>5} {'cohesion':>9}")
    print("  " + THIN)
    for r in results:
        bite = r.first_biting_depth()
        coh = r.first_cohesion_depth()
        print(f"  {r.key:<14} {r.structure.verdict:<18} {r.n_measured:>10,} "
              f"{_pct(r.pass_pct):>7} {r.n_hits:>7,} {r.n_hits_pass:>10,} "
              f"{_pct(r.hit_pass_pct):>10} "
              f"{('d' + str(bite)) if bite is not None else '-':>5} "
              f"{('d' + str(coh)) if coh is not None else '-':>9}")
    print()
    print("  veto column: LIVE = structure parsed and assayed positions mapped, so the")
    print("  veto really ran. PARSE_FAILED / NO_FILE / POSITIONS_UNMAPPED mean a 100.00%")
    print("  pass rate is the filter returning its input untouched, NOT an absence of")
    print("  spatial rejection. Do not read those rows as 'the veto is harmless'.")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("USAGE")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--library", action="append", default=None, metavar="KEY",
        help="limit to one library (repeatable). Default: every library in "
             "benchmark_panoptes_isolated.LIBRARIES that has a landscape and a structure.",
    )
    parser.add_argument(
        "--include-missing-structure", action="store_true",
        help="also profile libraries whose structure_pdb is absent, to show the "
             "unfiltered baseline explicitly.",
    )
    parser.add_argument(
        "--max-variants", type=int, default=0, metavar="N",
        help="cap the scope at N variants per library with a deterministic sample "
             "(0 = full universe, the default).",
    )
    parser.add_argument("--sample-trigger", type=int, default=DEFAULT_SAMPLE_TRIGGER,
                        help=f"auto-subsample above this universe size "
                             f"(default {DEFAULT_SAMPLE_TRIGGER:,}).")
    parser.add_argument("--sample-size", type=int, default=DEFAULT_SAMPLE_SIZE,
                        help=f"size of that auto-subsample (default {DEFAULT_SAMPLE_SIZE:,}).")
    parser.add_argument("--seed", type=int, default=0, help="subsample seed (default 0).")
    parser.add_argument("--json", type=Path, default=None, metavar="PATH",
                        help="also write every table as JSON.")
    args = parser.parse_args(argv)

    # Memoise the contact map for the duration of the run. The production function
    # constructs one per call (spatial_filter.py:704); over a whole universe that is
    # one PDB re-parse per batch. Nothing else about the code path is changed.
    sf.SpatialContactMap = _cached_contact_map
    try:
        requested = args.library or [k for k in bench.LIBRARIES if k in REGISTRY]
        unknown = [k for k in requested if k not in bench.LIBRARIES]
        if unknown:
            parser.error(f"unknown library key(s) {unknown}; "
                         f"choose from {sorted(bench.LIBRARIES)}")
        no_landscape = [k for k in requested if k not in REGISTRY]
        if no_landscape:
            print(f"NOTE: no landscape loader for {no_landscape}; skipping "
                  f"(veritas.calibration.landscapes.REGISTRY has "
                  f"{sorted(REGISTRY)}).")
            requested = [k for k in requested if k in REGISTRY]

        print(RULE)
        print("SPATIAL VETO PROFILE")
        print(RULE)
        print(f"  cohesion veto under test : mean pairwise C-alpha > {MEAN_THRESHOLD_A} A "
              f"OR max pairwise > {MAX_THRESHOLD_A} A  ->  REJECT")
        print("  code under test          : workers/spatial_filter.py::"
              "filter_spatially_viable_variants (called unconditionally from "
              "workers/math_worker.py::select_next_micro_batch)")
        print("  depths and hit sets      : veritas.calibration.landscapes.load_landscape "
              "(HANDOFF rule 2.1 - no new genotype parser)")
        print(f"  libraries in scope       : {', '.join(requested)}")
        anchors_ok = check_shallow_shell_anchors()

        results: List[LibraryResult] = []
        for key in requested:
            library = bench.LIBRARIES[key]
            if library.structure_pdb is None or not Path(library.structure_pdb).is_file():
                if not args.include_missing_structure:
                    print()
                    print(f"SKIP {key}: structure_pdb is "
                          f"{library.structure_pdb} which is not on disk. The filter would "
                          f"return every candidate unfiltered. Re-run with "
                          f"--include-missing-structure to profile it anyway.")
                    continue
            try:
                results.append(profile_library(
                    key,
                    max_variants=args.max_variants,
                    sample_trigger=args.sample_trigger,
                    sample_size=args.sample_size,
                    seed=args.seed,
                ))
            except Exception as exc:
                print()
                print(f"FAILED {key}: {type(exc).__name__}: {exc}")
                continue
            print_library(results[-1])

        if results:
            print_summary(results)
        if not anchors_ok:
            print()
            print("  *** rule-2.2 anchors did NOT reproduce; treat every table above as "
                  "unverified. ***")

        # Assert the production argument shape was never deviated from.
        shapes = {(c[1], c[2]) for c in _CTOR_CALLS}
        print()
        print(f"CONTACT-MAP CONSTRUCTIONS : {len(_CTOR_CALLS)} call(s), "
              f"{len(_MAP_CACHE)} distinct structure(s) parsed; "
              f"(chain_id, reference_sequence) arg shapes seen: {sorted(shapes)}")
        if shapes - {(None, None)}:
            print("  WARNING: a construction used arguments the production call site does "
                  "not pass; the cached map may not be the one production would build.")

        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            args.json.write_text(json.dumps(
                {
                    "cohesion_mean_threshold_A": MEAN_THRESHOLD_A,
                    "cohesion_max_threshold_A": MAX_THRESHOLD_A,
                    "libraries": [r.to_json() for r in results],
                }, indent=2), encoding="utf-8")
            print(f"JSON written to {args.json}")
    finally:
        sf.SpatialContactMap = _REAL_CONTACT_MAP
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

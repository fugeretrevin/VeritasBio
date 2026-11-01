"""Label-only loader for the nine landscapes of the cross-protein study.

Model-free by construction: it reads the prepared CSVs under `scripts/data/`, the
library definition, and - where a run exists - that run's own `manifest.json`. No
surrogate, no campaign, no selector.

Three rules are deliberate and are what make the figures reproduce:

1. **The threshold comes from the run's manifest where one exists.** `top1.threshold`
   is the 99th percentile over the eligible set as the benchmark computed it, so
   analysis scores variants with exactly the threshold the benchmark used
   (`scripts/analysis/README.md`, "A note on thresholds"). Where no run exists - the six
   CombinGym landscapes have never been run - the same definition is recomputed and the
   source is recorded as `computed`, never silently mixed with a manifest value.

2. **The wild type is not a library member.** A row canonicalising to the all-parent
   genotype is dropped, which is the benchmark's own rule and what makes GB1's universe
   149,360 rather than 149,361. Every eligible variant therefore has depth >= 1, and
   the shells partition the mutants.

3. **HER2's unlabelled split never enters the universe.** Its 25,114 `stage != train`
   rows carry `DMS_score` exactly 0.0, which is an absent label, not a measurement of
   no binding. Each discarded row is CHECKED against that placeholder, so a real
   measurement cannot be thrown away as an absent label.
"""
from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from veritas.calibration.endpoint import Hit, diversity_ceiling
from veritas.calibration.genotypes import Combo, canonical_combo, depth

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "scripts" / "data"
RUNS_DIR = REPO_ROOT / "benchmark_comparison_results"


@dataclass(frozen=True)
class LandscapeSpec:
    """One landscape's library definition and eligibility rules.

    `positions` and `parent_residues` for TrpB, GB1 and HER2 are the benchmark's own
    (`benchmark_panoptes_isolated.py:LIBRARIES`). For the six CombinGym landscapes they
    are the values cross-validated in `docs/dataset-inventory.md` §1 against the
    registry, the varying sequence indices and the wild-type FASTA; `load_landscape`
    re-checks both against the data on every load, so a silent mis-mapping cannot pass.
    """

    key: str
    label: str
    csv_name: str
    positions: Tuple[int, ...]
    parent_residues: str
    protein: str
    labelled_stage: Optional[str] = None
    unlabelled_placeholder: Optional[float] = None
    #: Deepened tier where the protocol's tier ladder (§9.1) deepens it; 0.01 otherwise.
    tier_fraction: float = 0.01

    @property
    def n_sites(self) -> int:
        return len(self.positions)

    @property
    def path(self) -> Path:
        return DATA_DIR / self.csv_name


REGISTRY: Dict[str, LandscapeSpec] = {
    "trpb": LandscapeSpec(
        key="trpb", label="TrpB four-site (Johnston et al. 2024)",
        csv_name="TRPB_THEMA_Johnston_2024.csv",
        positions=(183, 184, 227, 228), parent_residues="VFVS", protein="TrpB",
    ),
    "gb1": LandscapeSpec(
        key="gb1", label="GB1 four-site (Wu et al. 2016)",
        csv_name="GB1_STRSG_Wu_2016.csv",
        positions=(39, 40, 41, 54), parent_residues="VDGV", protein="GB1",
    ),
    "her2": LandscapeSpec(
        key="her2", label="Trastuzumab CDR-H3 ten-site, HER2 binding (Mason et al. 2021)",
        csv_name="HER2_TRASTUZUMAB_Mason_2021.csv",
        positions=(99, 100, 101, 102, 103, 104, 105, 106, 107, 108),
        parent_residues="WGGDGFYAMD", protein="trastuzumab",
        labelled_stage="train", unlabelled_placeholder=0.0,
    ),
    "eqfp611_blue": LandscapeSpec(
        key="eqfp611_blue", label="eqFP611/mTagBFP2 13-site, blue channel (Poelwijk et al. 2019)",
        csv_name="EQFP611_BLUE_Poelwijk_2019.csv",
        positions=(20, 45, 63, 127, 143, 158, 168, 172, 174, 197, 206, 207, 231),
        parent_residues="DVLTFNSAAYNNK", protein="eqFP611",
    ),
    "eqfp611_red": LandscapeSpec(
        key="eqfp611_red", label="eqFP611/mTagBFP2 13-site, red channel (Poelwijk et al. 2019)",
        csv_name="EQFP611_RED_Poelwijk_2019.csv",
        positions=(20, 45, 63, 127, 143, 158, 168, 172, 174, 197, 206, 207, 231),
        parent_residues="DVLTFNSAAYNNK", protein="eqFP611",
    ),
    "creilov": LandscapeSpec(
        key="creilov", label="CreiLOV 15-site fluorescence (Chen et al. 2023)",
        csv_name="CREILOV_CHLRE_Chen_2023.csv",
        positions=(3, 4, 5, 7, 29, 34, 47, 60, 61, 92, 96, 98, 107, 109, 113),
        parent_residues="GLRTAGQRDIDRVVT", protein="CreiLOV",
    ),
    "cr9114_h1": LandscapeSpec(
        key="cr9114_h1", label="CR9114 16-site, H1 binding (Phillips et al. 2021)",
        csv_name="CR9114_H1_Phillips_2021.csv",
        positions=(29, 30, 31, 52, 57, 58, 59, 71, 74, 75, 76, 77, 84, 87, 95, 106),
        parent_residues="FSSITANTKSTSSRYY", protein="CR9114",
    ),
    "cr6261_h1": LandscapeSpec(
        key="cr6261_h1", label="CR6261 11-site, H1 binding (Phillips et al. 2021)",
        csv_name="CR6261_H1_Phillips_2021.csv",
        positions=(28, 30, 58, 59, 62, 74, 75, 76, 77, 79, 104),
        parent_residues="TSANQKSTSAL", protein="CR6261",
        # §4.5: the protocol's tier ladder deepens CR6261 to 2%.
        tier_fraction=0.02,
    ),
    "cr6261_h9": LandscapeSpec(
        key="cr6261_h9", label="CR6261 11-site, H9 binding (Phillips et al. 2021)",
        csv_name="CR6261_H9_Phillips_2021.csv",
        positions=(28, 30, 58, 59, 62, 74, 75, 76, 77, 79, 104),
        parent_residues="TSANQKSTSAL", protein="CR6261",
        tier_fraction=0.02,
    ),
}


@dataclass(frozen=True)
class Shell:
    """A depth shell's size, hit count and hit rate."""

    n_variants: int
    n_hits: int
    depth_cut: int

    @property
    def rate(self) -> float:
        return self.n_hits / self.n_variants if self.n_variants else 0.0


@dataclass
class Landscape:
    """One landscape's eligible universe, in label space only."""

    spec: LandscapeSpec
    combos: List[Combo]
    fitness: np.ndarray
    depths: np.ndarray
    threshold: float
    threshold_source: str
    computed_threshold: float
    n_rows_read: int = 0
    n_dropped_unparseable: int = 0
    n_dropped_unlabelled: int = 0
    n_dropped_wild_type: int = 0
    n_dropped_nonfinite: int = 0
    n_dropped_duplicate: int = 0
    _ceiling_cache: Dict[int, int] = field(default_factory=dict, repr=False)

    # ---- universe ------------------------------------------------------------------
    @property
    def key(self) -> str:
        return self.spec.key

    @property
    def n_universe(self) -> int:
        return len(self.combos)

    @property
    def hit_mask(self) -> np.ndarray:
        """At-or-above semantics: every variant at the threshold counts as a hit."""
        return self.fitness >= self.threshold

    @property
    def n_hits(self) -> int:
        return int(np.count_nonzero(self.hit_mask))

    @property
    def hits(self) -> List[Hit]:
        """(fitness, combo) for every hit, fitness descending."""
        idx = np.flatnonzero(self.hit_mask)
        pairs = [(float(self.fitness[i]), self.combos[i]) for i in idx]
        pairs.sort(key=lambda p: -p[0])
        return pairs

    @property
    def hit_combos(self) -> List[Combo]:
        return [c for _f, c in self.hits]

    # ---- shells --------------------------------------------------------------------
    def shell(self, depth_cut: int) -> Shell:
        """`1 <= depth <= depth_cut`. Depth 0 is not in the universe at all (rule 2)."""
        sel = (self.depths >= 1) & (self.depths <= depth_cut)
        return Shell(
            n_variants=int(np.count_nonzero(sel)),
            n_hits=int(np.count_nonzero(sel & self.hit_mask)),
            depth_cut=depth_cut,
        )

    def deep_shell(self, depth_cut: int) -> Shell:
        sel = self.depths > depth_cut
        return Shell(
            n_variants=int(np.count_nonzero(sel)),
            n_hits=int(np.count_nonzero(sel & self.hit_mask)),
            depth_cut=depth_cut,
        )

    @property
    def primary_depth_cut(self) -> int:
        """§6.3's pre-specified shell cut: `floor(S/2)`."""
        return self.spec.n_sites // 2

    def depth_histogram(self) -> Dict[int, Tuple[int, int]]:
        """`{depth: (hits, variants)}`."""
        out: Dict[int, Tuple[int, int]] = {}
        for d in range(int(self.depths.min()), int(self.depths.max()) + 1):
            sel = self.depths == d
            n = int(np.count_nonzero(sel))
            if n:
                out[d] = (int(np.count_nonzero(sel & self.hit_mask)), n)
        return out

    # ---- endpoint ------------------------------------------------------------------
    def ceilings(self, distances: Sequence[int] = (2, 3, 4)) -> Dict[int, int]:
        """Diversity-adjusted ceiling at each distance, over the whole hit set."""
        wanted = [int(d) for d in distances if int(d) not in self._ceiling_cache]
        if wanted:
            self._ceiling_cache.update(diversity_ceiling(self.hits, wanted))
        return {int(d): self._ceiling_cache[int(d)] for d in distances}

    @property
    def tier_size(self) -> int:
        """`ceil(tier_fraction * N)` - the hit set the tier ladder implies, for contrast."""
        return int(math.ceil(self.spec.tier_fraction * self.n_universe))


#: Where a run manifest records the top-1% threshold. Current schema first.
#:
#: Both shapes are read because this function used to know only the legacy one, and the
#: consequence was silent and forbidden: on no match the caller RECOMPUTES the 99th
#: percentile, which rule 2.3 prohibits outright ("Read `top1.threshold` from each run's
#: own `manifest.json`. Do not recompute it."). Today's manifests are
#: `manifest_schema_version` 2.x and record `results.top1.threshold_value`, so every
#: modern run was being skipped. It went unnoticed only because 160 legacy run
#: directories survive in the main checkout and satisfied the legacy path; it would have
#: begun failing unconditionally the moment they were pruned, and it already fails inside
#: a git worktree, where the gitignored results root is not visible.
MANIFEST_THRESHOLD_PATHS: Tuple[Tuple[str, ...], ...] = (
    ("results", "top1", "threshold_value"),
    ("top1", "threshold"),
)

#: Where a run manifest records the library it scored against. Current schema first.
MANIFEST_LIBRARY_PATHS: Tuple[Tuple[str, ...], ...] = (
    ("configuration", "dataset", "library"),
    ("dataset",),
)

#: Where the synthetic-dataset flag lives. Current schema first.
MANIFEST_SYNTHETIC_PATHS: Tuple[Tuple[str, ...], ...] = (
    ("configuration", "dataset", "synthetic"),
    ("dataset", "synthetic"),
)


def _manifest_dig(mapping: Any, path: Sequence[str]) -> Any:
    node: Any = mapping
    for key in path:
        if not isinstance(node, Mapping) or key not in node:
            return None
        node = node[key]
    return node


def _manifest_first(mapping: Any, paths: Sequence[Sequence[str]]) -> Any:
    for path in paths:
        found = _manifest_dig(mapping, path)
        if found is not None:
            return found
    return None


def _describes_landscape(library: Mapping[str, Any], spec: LandscapeSpec) -> bool:
    """Whether a manifest's library block is this landscape's.

    Matched on positions AND parent residues rather than on the key, because that works
    identically in both manifest schemas and is strictly stronger - it catches a library
    renamed under the same geometry, and a key reused for a different geometry. The key
    is checked too when present, so a rename cannot pass silently either.
    """
    positions = library.get("positions", library.get("active_site_positions"))
    if tuple(int(p) for p in (positions or ())) != spec.positions:
        return False
    if str(library.get("parent_residues") or "") != spec.parent_residues:
        return False
    key = library.get("key")
    return key is None or str(key) == spec.key


def _manifest_threshold(spec: LandscapeSpec) -> Optional[Tuple[float, str]]:
    """`top1.threshold` from this landscape's own runs, or None when none exist.

    Every non-synthetic manifest whose library matches is read and they must AGREE:
    two runs of one landscape scoring against different thresholds would mean the
    eligible universe moved between them, which is a finding, not something to average.

    Returning None sends the caller to a recomputed 99th percentile, which rule 2.3
    forbids - so the two reasons for None are kept apart and reported. "No run exists
    yet" is legitimate and is exactly the state a newly reserved landscape is in. "Runs
    exist for this landscape but their threshold could not be read" is a defect, and it
    RAISES rather than falling through, because that is the case where a real threshold
    is silently replaced by an invented one.
    """
    if not RUNS_DIR.is_dir():
        return None

    found: Dict[float, List[str]] = {}
    n_read = 0
    n_synthetic = 0
    n_other_landscape = 0
    matched_without_threshold: List[str] = []

    for manifest_path in sorted(RUNS_DIR.glob("*/manifest.json")):
        try:
            m = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        n_read += 1
        if _manifest_first(m, MANIFEST_SYNTHETIC_PATHS):
            n_synthetic += 1
            continue
        library = _manifest_first(m, MANIFEST_LIBRARY_PATHS)
        if not isinstance(library, Mapping) or not _describes_landscape(library, spec):
            n_other_landscape += 1
            continue
        threshold = _manifest_first(m, MANIFEST_THRESHOLD_PATHS)
        if threshold is None or isinstance(threshold, Mapping):
            matched_without_threshold.append(manifest_path.parent.name)
            continue
        found.setdefault(float(threshold), []).append(manifest_path.parent.name)

    if matched_without_threshold and not found:
        raise ValueError(
            f"{spec.key}: {len(matched_without_threshold)} run manifest(s) describe this "
            f"landscape but record no threshold at any known path "
            f"({', '.join('.'.join(p) for p in MANIFEST_THRESHOLD_PATHS)}) - e.g. "
            f"{matched_without_threshold[0]}. Rule 2.3 forbids recomputing it, so this "
            f"cannot fall back to a percentile. Add the manifest's actual path to "
            f"MANIFEST_THRESHOLD_PATHS, or establish why a scored run recorded none."
        )

    if not found:
        # Legitimate: no run has scored this landscape yet. The caller computes a
        # bootstrap percentile and the source string says so, naming what was read -
        # a reader must be able to tell "nothing exists" from "nothing was recognised".
        return None

    if len(found) > 1:
        detail = "; ".join(
            f"{t!r} in {len(runs)} run(s) e.g. {runs[0]}" for t, runs in sorted(found.items())
        )
        raise ValueError(
            f"{spec.key}: run manifests disagree on top1.threshold - {detail}. "
            f"Resolve which eligible universe is correct before scoring anything."
        )
    threshold, runs = next(iter(found.items()))
    return threshold, f"manifest ({len(runs)} run(s), e.g. {runs[0]})"


def _no_manifest_reason(spec: LandscapeSpec) -> str:
    """Why no manifest threshold was available, for the `threshold_source` string.

    The old string was the single literal "computed (no run manifest exists)", which was
    FALSE whenever manifests existed but none were recognised - the exact condition the
    schema drift created. A source field that lies about its own provenance is worse than
    one that is merely unhelpful.
    """
    if not RUNS_DIR.is_dir():
        return f"computed (no results directory at {RUNS_DIR})"
    n_manifests = sum(1 for _ in RUNS_DIR.glob("*/manifest.json"))
    if not n_manifests:
        return "computed (no run manifest exists)"
    return (f"computed ({n_manifests} manifest(s) read, none scored {spec.key}) - "
            f"NOT a measured threshold")


def load_landscape(key: str, *, strict_rows: bool = True) -> "Landscape":
    """Load one landscape's eligible universe, verifying the library against the data.

    `strict_rows` raises when a row does not map onto the declared design, because a
    dataset whose rows disagree with its library definition is a configuration error.
    """
    return _load_cached(key, strict_rows)


@lru_cache(maxsize=None)
def _load_cached(key: str, strict_rows: bool) -> "Landscape":
    try:
        spec = REGISTRY[key]
    except KeyError:
        raise KeyError(f"unknown landscape {key!r}; known: {sorted(REGISTRY)}") from None
    if not spec.path.is_file():
        raise FileNotFoundError(f"{spec.key}: prepared dataset not found at {spec.path}")

    combos: List[Combo] = []
    values: List[float] = []
    seen: Dict[Combo, int] = {}
    n_rows = n_unparseable = n_unlabelled = n_wt = n_nonfinite = n_dup = 0
    observed_positions: Dict[int, str] = {}

    with spec.path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        columns = reader.fieldnames or []
        id_col = _pick(columns, ("mutant",))
        value_col = _pick(columns, ("dms_score", "fitness", "score"))
        stage_col = _pick(columns, ("stage",))
        if id_col is None or value_col is None:
            raise ValueError(f"{spec.key}: no variant-id / value column in {columns}")
        for row in reader:
            n_rows += 1
            variant_id = (row.get(id_col) or "").strip()

            if spec.labelled_stage is not None:
                stage = (row.get(stage_col) or "").strip() if stage_col else ""
                if stage != spec.labelled_stage:
                    # Verify the discarded row really is unlabelled before discarding it.
                    if spec.unlabelled_placeholder is not None:
                        try:
                            held = float(row[value_col])
                        except (TypeError, ValueError, KeyError):
                            held = float("nan")
                        if held != spec.unlabelled_placeholder:
                            raise ValueError(
                                f"{spec.key}: row {variant_id!r} in split {stage!r} carries "
                                f"{held!r}, not the declared unlabelled placeholder "
                                f"{spec.unlabelled_placeholder!r} - that is a measurement "
                                f"about to be discarded as an absent label"
                            )
                    n_unlabelled += 1
                    continue

            combo = canonical_combo(variant_id, spec.positions, spec.parent_residues)
            if combo is None:
                if strict_rows and variant_id and variant_id.upper() != "WT":
                    raise ValueError(
                        f"{spec.key}: row {variant_id!r} does not map onto the declared "
                        f"library {spec.positions} / {spec.parent_residues!r}"
                    )
                n_unparseable += 1
                continue
            for pos, _wt, _mut in _tokens_for_check(variant_id):
                observed_positions.setdefault(pos, _wt)

            if depth(combo, spec.parent_residues) == 0:
                # The parent is not a library member (benchmark rule); GB1 149,360.
                n_wt += 1
                continue
            try:
                value = float(row[value_col])
            except (TypeError, ValueError, KeyError):
                n_nonfinite += 1
                continue
            if not math.isfinite(value):
                n_nonfinite += 1
                continue
            if combo in seen:
                n_dup += 1
                continue
            seen[combo] = len(combos)
            combos.append(combo)
            values.append(value)

    if not combos:
        raise ValueError(f"{spec.key}: no eligible variants parsed from {spec.path}")

    # Registry cross-validation: the positions that actually vary must be exactly the
    # declared ones, and the parent residue at each must agree with the data's own.
    if set(observed_positions) != set(spec.positions):
        raise ValueError(
            f"{spec.key}: data varies positions {sorted(observed_positions)} but the "
            f"library declares {list(spec.positions)}"
        )
    spelled = "".join(observed_positions[p] for p in spec.positions)
    if spelled != spec.parent_residues:
        raise ValueError(
            f"{spec.key}: the data's wild-type residues spell {spelled!r}, the library "
            f"declares {spec.parent_residues!r}"
        )

    fitness = np.asarray(values, dtype=np.float64)
    # The percentile must follow the landscape's DECLARED tier, not a hard-coded 99.0.
    # This was `np.percentile(fitness, 99.0)` for every landscape, while the `tier_size`
    # property used `tier_fraction` - so any landscape declaring a tier other than 1%
    # had a `threshold` and a `tier_size` describing different tiers. It bit exactly two:
    # CR6261-H1 and CR6261-H9 both declare `tier_fraction=0.02` under §4.5's tier ladder
    # and were handed a 1% bar, giving 19 hits where their own tier implies 38 and 37.
    # Neither has ever been scored - 0 of 170 run manifests mention them - so no measured
    # threshold was overridden and nothing downstream of a manifest moves.
    computed = float(np.percentile(fitness, 100.0 * (1.0 - spec.tier_fraction)))
    manifest = _manifest_threshold(spec)
    threshold, source = (manifest if manifest is not None
                         else (computed, _no_manifest_reason(spec)))
    depths = np.asarray(
        [depth(c, spec.parent_residues) for c in combos], dtype=np.int16
    )
    return Landscape(
        spec=spec, combos=combos, fitness=fitness, depths=depths,
        threshold=float(threshold), threshold_source=source, computed_threshold=computed,
        n_rows_read=n_rows, n_dropped_unparseable=n_unparseable,
        n_dropped_unlabelled=n_unlabelled, n_dropped_wild_type=n_wt,
        n_dropped_nonfinite=n_nonfinite, n_dropped_duplicate=n_dup,
    )


def _pick(columns: Sequence[str], wanted: Sequence[str]) -> Optional[str]:
    lowered = {c.lower(): c for c in columns}
    for name in wanted:
        if name in lowered:
            return lowered[name]
    return None


def _tokens_for_check(variant_id: str) -> List[Tuple[int, str, str]]:
    from veritas.calibration.genotypes import parse_mutations

    return parse_mutations(variant_id) or []

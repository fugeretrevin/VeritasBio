"""The pre-flight gate: assertions on the artefacts that already recorded the defect.

Five separate defects in this project invalidated or nearly invalidated a result, and
every one of them was **already written down in an artefact we produced ourselves**:

| defect | the artefact that already held the number |
|---|---|
| ALDE a null arm on 7 of 9 landscapes, 26 unique measurements of a 125 budget | `round_metrics.csv` |
| our own arms measuring 120-124 where baselines measured 125 | `round_metrics.csv` |
| CR9114 missing 25% of its variants with depth-correlated attrition | `docs/dataset-inventory.md` |
| the primary endpoint degenerate on HER2 (d>=2 ceiling 113 == raw ceiling 113) | both numbers published |
| the qualification bar gating nothing on 8 of 9 landscapes | every HER2 seed logged "qualified before round 1" |

So this is not a detection problem, it is an **assertion** problem, and this module is
the assertions. Six checks:

1. ``corpus``    - every file a landscape reads is hashed in the corpus lock, and a
                   derived "clean" table is flagged when a raw one exists at the same pin.
2. ``encoding``  - no constant column in any encoded design space. A constant one-hot
                   column is a residue never observed at that position, which became
                   ``0/0 = NaN`` inside ALDE and silently degenerated it.
3. ``budget``    - distinct measurements equal the declared budget, and a shortfall is
                   classified DEGENERATE (could not propose a distinct candidate;
                   unrecoverable) or RECOVERABLE (proposals ineligible, backfilled).
4. ``endpoint``  - the diversity-adjusted ceiling differs materially from the raw count,
                   and the budget is smaller than the hit tier.
5. ``threshold`` - the qualification bar is a quantile resolved from the same values
                   array that defines ``top1.threshold``, and the manifest records the
                   quantile, the resolved absolute and the source.
6. ``row_order`` - shuffling the input table's row order and re-running one seed gives
                   identical results. This one generalises furthest: our prepared tables
                   are sorted by fitness descending, so CR9114's row 0 IS the global
                   optimum, and an arm whose surrogate is entirely NaN would have posted
                   a perfect best-fitness there and beaten a working one.

Design rules, each of them a defect we already paid for:

* **No check returns a bare boolean.** Every check returns a :class:`CheckResult` with a
  status of PASS / FAIL / NOT-APPLICABLE / SKIPPED, the measured value, the expectation,
  and - for anything that is not a pass - a reason. ``bool(result)`` raises, because
  ``if check(...):`` collapsing "not applicable" into "fine" is the whole failure mode.
* **Skips are recorded, never silent.** A check that did not run says so, and says why.
* **Nothing is reimplemented.** The degeneracy rule comes from
  :mod:`veritas.calibration.endpoint`, the quantile conversion from
  :mod:`veritas.calibration.qualification`, genotype canonicalisation from
  :mod:`veritas.calibration.genotypes`, and corpus hashing from
  ``scripts/freeze_cross_protein_corpus.py``, which is *called*, not copied.
* **Mutation strings omit wild-type positions.** ``V183A:F184A:V227A`` is a four-site
  genotype whose fourth site is wild type. Everything here works on canonical
  full-length combos from :mod:`veritas.calibration.genotypes` and asserts equal length;
  nothing is zipped or truncated.

Two of the six are not purely static: ``budget`` reads a run's ``round_metrics.csv`` and
``row_order`` needs a runner it can call twice. Both say so rather than passing silently
when the data is absent.

Usage::

    .venv\\Scripts\\python.exe -X utf8 -m veritas.preflight --library her2 --arms alde,ridge
    .venv\\Scripts\\python.exe -X utf8 -m veritas.preflight --library gb1 --arms alde \\
        --run-dir benchmark_comparison_results/run_20261031_132813_ff1a3423
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import importlib
import importlib.util
import io
import json
import random
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import (Any, Callable, Dict, Iterable, List, Mapping, Optional,
                    Sequence, Set, Tuple)

import numpy as np

from veritas import paths
from veritas.calibration.endpoint import (DEGENERATE, INFORMATIVENESS_FLOOR,
                                          NEAR_DEGENERATE, USABLE, endpoint_verdict)
from veritas.calibration.qualification import (absolute_for_quantile,
                                               equivalent_quantile,
                                               p_qualify_in_draw, qualifying_count)

# ======================================================================================
# Statuses and check identifiers
# ======================================================================================
PASS = "PASS"
FAIL = "FAIL"
NOT_APPLICABLE = "NOT-APPLICABLE"
SKIPPED = "SKIPPED"

CHECK_CORPUS = "corpus"
CHECK_ENCODING = "encoding"
CHECK_BUDGET = "budget"
CHECK_ENDPOINT = "endpoint"
CHECK_THRESHOLD = "threshold"
CHECK_ROW_ORDER = "row_order"

#: Every check, in the order the gate runs them.
ALL_CHECKS: Tuple[str, ...] = (
    CHECK_CORPUS, CHECK_ENCODING, CHECK_BUDGET, CHECK_ENDPOINT, CHECK_THRESHOLD,
    CHECK_ROW_ORDER,
)

#: Checks that cannot be answered from the repository alone. Stated here so a report can
#: say which part of the gate a static invocation did not cover.
CHECKS_REQUIRING_A_RUN: Tuple[str, ...] = (CHECK_BUDGET, CHECK_ROW_ORDER)

# The 20 amino acids, in ALDE's own order (`src/encoding_utils.py:ALL_AAS`). Its one-hot
# encoder allocates a column per amino acid per site REGARDLESS of what the library
# contains, which is why a 2-residue site yields 18 structurally-zero columns.
ALL_AAS: Tuple[str, ...] = tuple("ACDEFGHIKLMNPQRSTVWY")

#: Arms handed an encoded design space whose pipeline normalises it by column range, and
#: therefore divide by zero on a constant column. The ALDE family is the known case
#: (`docs/alde-compatibility.md` §3); extend it with `encoded_arms=` for a new arm.
ENCODED_DESIGN_SPACE_ARMS: frozenset = frozenset({"alde", "alde_constant_column_guard"})

#: Constant one-hot columns per landscape, as measured against ALDE's own code in
#: `docs/alde-compatibility.md` §1. Reference figures: `encode_landscape` plus
#: `constant_columns` must reproduce them, or this check cannot catch the defect.
REFERENCE_CONSTANT_COLUMNS: Dict[str, Tuple[int, int]] = {
    "trpb": (0, 80),
    "gb1": (0, 80),
    "her2": (55, 200),
    "creilov": (265, 300),
    "cr9114_h1": (288, 320),
    "eqfp611_blue": (234, 260),
    "eqfp611_red": (234, 260),
    "cr6261_h1": (198, 220),
    "cr6261_h9": (198, 220),
}

#: sha256 of a prepared table under `scripts/data/`, where one has been published. The
#: corpus lock hashes the UPSTREAM corpus only, so a prepared table - the file an arm
#: actually reads - has no hash coverage unless it is pinned here.
#: cr9114_h1: `docs/cr9114-reprep.md` §4.
PREPARED_TABLE_DIGESTS: Dict[str, str] = {
    "cr9114_h1": "1d9a57745636f27d5fdea5f8798a852f4852b9dc7fa2e60884d49bd930972d83",
}

# Budget shortfall classifications. They are different failures and must not collapse.
BUDGET_DEGENERATE = "DEGENERATE"
BUDGET_RECOVERABLE = "RECOVERABLE"
BUDGET_UNCLASSIFIED = "UNCLASSIFIED"

_SEVERITY = {BUDGET_DEGENERATE: 3, BUDGET_UNCLASSIFIED: 2, BUDGET_RECOVERABLE: 1}

# A rejection reason naming any of these means the arm re-proposed something it had
# already measured or already had in the batch: it could not propose a distinct
# candidate, which is the degenerate case. Anything else is an eligibility rejection.
DEGENERACY_REASON_MARKERS: Tuple[str, ...] = (
    "duplicate", "already_observed", "already_measured", "repeat", "over_budget",
)

#: Fixed shuffle seed, recorded in the result so a row-order failure reproduces.
DEFAULT_SHUFFLE_SEED = 20261031

DEFAULT_BUDGET = 125
DEFAULT_INITIAL_SIZE = 25
DEFAULT_DISTANCE = 2
TOP1_QUANTILE = 99.0

#: A qualification bar that the 25-variant initial draw already clears with at least this
#: probability does not gate: the stage is satisfied before round 1 on more than half the
#: seeds. On HER2 the absolute 0.70 default qualifies 1,027 of 11,276 variants, so the
#: exact hypergeometric probability is **0.908** - which is why every HER2 seed logged
#: "qualified before round 1". A strict `p == 1.0` test would have missed that, and this
#: is a judgement bar rather than a derived one, so it is a named, overridable constant.
GATING_P_MAX = 0.5


# ======================================================================================
# The result type
# ======================================================================================
@dataclass(frozen=True)
class CheckResult:
    """One check's structured outcome. Never a bare boolean.

    `measured` and `expected` are what make a failure actionable: a report that says
    "budget check failed" sends the next reader to the code, while one that says
    "26 distinct measurements against a declared 125" sends them to the arm.
    """

    check: str
    status: str
    measured: Any
    expected: Any
    summary: str
    reason: str = ""
    detail: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.status not in (PASS, FAIL, NOT_APPLICABLE, SKIPPED):
            raise ValueError(f"unknown check status {self.status!r}")
        if self.status in (NOT_APPLICABLE, SKIPPED) and not str(self.reason).strip():
            raise ValueError(
                f"{self.check}: a {self.status} result must carry a reason - a check that "
                f"did not run is recorded, never silent"
            )

    def __bool__(self) -> bool:  # pragma: no cover - the raise IS the behaviour
        raise TypeError(
            f"{self.check}: a CheckResult is not a bare boolean. Collapsing "
            f"{NOT_APPLICABLE} or {SKIPPED} into 'fine' is the failure mode this gate "
            f"exists to prevent - read .status, .passed or .failed."
        )

    @property
    def passed(self) -> bool:
        return self.status == PASS

    @property
    def failed(self) -> bool:
        return self.status == FAIL

    @property
    def applicable(self) -> bool:
        return self.status in (PASS, FAIL)

    def render(self) -> str:
        line = f"{self.status:<14} {self.check:<42} {self.summary}"
        if self.applicable:
            line += f"\n{'':<14} measured {self.measured!r}; expected {self.expected!r}"
        if self.reason:
            line += f"\n{'':<14} reason: {self.reason}"
        return line

    def to_dict(self) -> Dict[str, Any]:
        return {
            "check": self.check,
            "status": self.status,
            "measured": _jsonable(self.measured),
            "expected": _jsonable(self.expected),
            "summary": self.summary,
            "reason": self.reason,
            "detail": _jsonable(dict(self.detail)),
        }


@dataclass(frozen=True)
class PreflightReport:
    """Every check's result for one (library, arm set), and the gate's exit code."""

    library: str
    arms: Tuple[str, ...]
    results: Tuple[CheckResult, ...]

    @property
    def failures(self) -> Tuple[CheckResult, ...]:
        return tuple(r for r in self.results if r.failed)

    @property
    def skipped(self) -> Tuple[CheckResult, ...]:
        return tuple(r for r in self.results if r.status == SKIPPED)

    @property
    def not_applicable(self) -> Tuple[CheckResult, ...]:
        return tuple(r for r in self.results if r.status == NOT_APPLICABLE)

    @property
    def exit_code(self) -> int:
        return 1 if self.failures else 0

    def render(self) -> str:
        head = (
            f"pre-flight gate: library={self.library} "
            f"arms={','.join(self.arms) if self.arms else '(none)'}"
        )
        lines = [head, "=" * len(head)]
        lines.extend(r.render() for r in self.results)
        counts = {s: sum(1 for r in self.results if r.status == s)
                  for s in (PASS, FAIL, NOT_APPLICABLE, SKIPPED)}
        lines.append("")
        lines.append(
            f"{counts[PASS]} passed, {counts[FAIL]} FAILED, "
            f"{counts[NOT_APPLICABLE]} not applicable, {counts[SKIPPED]} skipped"
        )
        if self.failures:
            lines.append("FAILED: " + ", ".join(r.check for r in self.failures))
        return "\n".join(lines)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema": "veritasbio/preflight-report/1",
            "library": self.library,
            "arms": list(self.arms),
            "exit_code": self.exit_code,
            "results": [r.to_dict() for r in self.results],
        }


def result_for(results: Iterable[CheckResult], check: str) -> CheckResult:
    """The one result with this check id. Raises rather than returning None."""
    materialised = list(results)
    for result in materialised:
        if result.check == check:
            return result
    raise KeyError(f"no result for {check!r}; have {[r.check for r in materialised]}")


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return [_jsonable(v) for v in value.tolist()]
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def _skipped(check: str, reason: str) -> CheckResult:
    return CheckResult(check=check, status=SKIPPED, measured=None, expected=None,
                       summary="not run", reason=reason)


# ======================================================================================
# 1. Corpus integrity
# ======================================================================================
def _freeze_module():
    """`scripts/freeze_cross_protein_corpus.py`, imported rather than reimplemented."""
    script = paths.repo_path("scripts/freeze_cross_protein_corpus.py")
    spec = importlib.util.spec_from_file_location("_veritas_freeze_corpus", script)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        raise ImportError(f"cannot load the freeze script at {script}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _normalise_table_token(name: str) -> str:
    """`bnAbs_CR6261_H1_clean.xlsx` -> `bnabscr6261h1`; `bnAbs_CR9114_raw.csv` -> `bnabscr9114`."""
    stem = Path(name).stem.lower()
    for suffix in ("_clean", "_raw"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
    return "".join(c for c in stem if c.isalnum())


def check_corpus_integrity(
    *,
    corpus_dir: Optional[Path] = None,
    lock_path: Optional[Path] = None,
    landscape_files: Optional[Mapping[str, str]] = None,
    acknowledged_derived: Optional[Mapping[str, str]] = None,
    prepared_digests: Optional[Mapping[str, str]] = None,
) -> List[CheckResult]:
    """Three results: the lock, the derived-table flag, and the prepared-table pin.

    `corpus.lock` is `scripts/freeze_cross_protein_corpus.py verify` - called, not
    reimplemented - which already checks that every pinned file is present, that every
    file a landscape reads is in the lock, and that both hashes match.

    The two results after it are what `verify` does **not** cover:

    * `corpus.derived_table` - a landscape reading a derived "clean" table while a raw
      table exists at the same pinned commit. CombinGym's clean CR9114 workbook silently
      excluded every genotype carrying both F29S and Y106S: 16,384 of them, a whole
      quadrant, with depth-correlated attrition that reached 0% at depth 16 and so
      removed the matured antibody the landscape is about (`docs/cr9114-reprep.md`).
      `verify` cannot see this, because the clean file it hashes is exactly the file that
      was read. Acknowledge a reviewed case with `acknowledged_derived={key: why}`; the
      acknowledgement is recorded in the result.
    * `corpus.prepared` - the sha256 of a prepared table under `scripts/data/`. The lock
      pins the upstream corpus; the arms read these derived CSVs, which it never hashes.
    """
    freeze = _freeze_module()
    corpus_dir = Path(corpus_dir) if corpus_dir is not None else Path(freeze.DEFAULT_CORPUS_DIR)
    lock_path = Path(lock_path) if lock_path is not None else Path(freeze.DEFAULT_LOCK)
    results: List[CheckResult] = []

    # ---- corpus.lock: the freeze script's own verify ----------------------------------
    buffer = io.StringIO()
    try:
        with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
            exit_code = int(freeze.main(
                ["--corpus-dir", str(corpus_dir), "verify", "--lock", str(lock_path)]))
    except BaseException as exc:  # noqa: BLE001 - a crash is a failed check, not a crash
        exit_code = 2
        buffer.write(f"\n{type(exc).__name__}: {exc}")
    output = buffer.getvalue().strip()
    results.append(CheckResult(
        check="corpus.lock",
        status=PASS if exit_code == 0 else FAIL,
        measured=exit_code,
        expected=0,
        summary=("every pinned file present and both hashes match"
                 if exit_code == 0 else
                 "freeze_cross_protein_corpus.py verify reported a corpus that is not "
                 "the one the protocol was frozen against"),
        detail={
            "verify_exit_code": exit_code,
            "verify_output": output,
            "corpus_dir": corpus_dir,
            "lock_path": lock_path,
            "delegated_to": "scripts/freeze_cross_protein_corpus.py verify",
        },
    ))

    # ---- corpus.derived_table ---------------------------------------------------------
    # Default scope: the landscapes this evaluation actually reads. A `reserve` landscape
    # is not part of the two-level validation assignment, so flagging its table would
    # bury the finding that matters under rows nothing scores against.
    files = dict(landscape_files) if landscape_files is not None else {
        key: spec["file"] for key, spec in freeze.LANDSCAPES.items()
        if spec.get("level") != "reserve"
    }
    acknowledged = dict(acknowledged_derived or {})
    known_raw = [rel for rel, _blob, _size in freeze.PINNED_FILES if "/Raw/" in rel]
    if corpus_dir.is_dir():
        raw_dir = corpus_dir / "Data" / "DMS" / "Raw"
        if raw_dir.is_dir():
            known_raw.extend(
                f"Data/DMS/Raw/{p.name}" for p in sorted(raw_dir.iterdir()) if p.is_file()
            )
    raw_tokens = {_normalise_table_token(rel): rel for rel in sorted(set(known_raw))}

    flagged: List[Dict[str, str]] = []
    for key, rel in sorted(files.items()):
        if "/Clean/" not in str(rel):
            continue
        token = _normalise_table_token(rel)
        alternative = next(
            (raw for raw_token, raw in sorted(raw_tokens.items())
             if token.startswith(raw_token) or raw_token.startswith(token)),
            None,
        )
        if alternative is None or key in acknowledged:
            continue
        flagged.append({"landscape": key, "derived_table": str(rel),
                        "raw_alternative": alternative})
    results.append(CheckResult(
        check="corpus.derived_table",
        status=PASS if not flagged else FAIL,
        measured=len(flagged),
        expected=0,
        summary=("no landscape reads a derived table that has a raw counterpart at the pin"
                 if not flagged else
                 f"{len(flagged)} landscape(s) read a derived 'clean' table while a raw "
                 f"table exists at the same pinned commit - the CR9114 failure mode"),
        detail={
            "flagged": flagged,
            "acknowledged": acknowledged,
            "n_landscapes_examined": len(files),
            "landscapes_examined": sorted(files),
            "scope": ("landscapes supplied by the caller" if landscape_files is not None
                      else "every non-reserve landscape in the freeze script's registry"),
            "raw_tables_at_the_pin": sorted(raw_tokens.values()),
            "motivating_defect": "docs/cr9114-reprep.md",
        },
    ))

    # ---- corpus.prepared --------------------------------------------------------------
    pins = dict(PREPARED_TABLE_DIGESTS if prepared_digests is None else prepared_digests)
    if not pins:
        results.append(CheckResult(
            check="corpus.prepared", status=NOT_APPLICABLE, measured=None, expected=None,
            summary="no prepared-table digest is pinned",
            reason=("the corpus lock hashes the UPSTREAM corpus only, so the prepared "
                    "tables the arms actually read have no hash coverage; pin them in "
                    "PREPARED_TABLE_DIGESTS or pass prepared_digests="),
            detail={"lock_covers": "upstream corpus files only"},
        ))
        return results

    from veritas.calibration.landscapes import REGISTRY

    checked: Dict[str, str] = {}
    mismatched: Dict[str, Dict[str, str]] = {}
    missing: List[str] = []
    for key, pinned in sorted(pins.items()):
        spec = REGISTRY.get(key)
        table = Path(spec.path) if spec is not None else paths.dataset_dir() / f"{key}.csv"
        if not table.is_file():
            missing.append(f"{key}: {table}")
            continue
        observed = hashlib.sha256(table.read_bytes()).hexdigest()
        checked[key] = observed
        if observed != pinned:
            mismatched[key] = {"pinned": pinned, "observed": observed,
                               "path": str(table)}
    bad = bool(mismatched or missing)
    results.append(CheckResult(
        check="corpus.prepared",
        status=FAIL if bad else PASS,
        measured=len(mismatched) + len(missing),
        expected=0,
        summary=(f"{len(checked)} prepared table(s) match their pinned sha256"
                 if not bad else
                 f"{len(mismatched)} prepared table(s) do not match their pinned sha256 "
                 f"and {len(missing)} are absent"),
        detail={"checked": checked, "mismatched": mismatched, "missing": missing},
    ))
    return results


# ======================================================================================
# 2. Encoding compatibility
# ======================================================================================
def one_hot(combos: Sequence[Sequence[str]],
            alphabet: Sequence[str] = ALL_AAS) -> np.ndarray:
    """ALDE's encoding: `len(sites) x len(ALL_AAS)` columns, whatever the library holds.

    `src/encoding_utils.py:46` builds its index from the module-level 20-amino-acid
    tuple and allocates a column per amino acid per site regardless of the library's
    actual alphabet, so a binary site contributes 18 structurally-zero columns. That
    waste is what later becomes a `0/0`.

    Genotypes must already be canonical full-length combos over the library's positions
    (`veritas.calibration.genotypes.canonical_combo`). Unequal lengths raise: a mutation
    string omits wild-type positions, and zipping or truncating them is the bug that
    inverted a conclusion earlier in this project.
    """
    rows = [tuple(c) for c in combos]
    if not rows:
        raise ValueError("cannot encode an empty design space")
    widths = {len(r) for r in rows}
    if len(widths) != 1:
        raise ValueError(
            f"genotypes of unequal length {sorted(widths)}: canonicalise to full library "
            f"length first - never zip or truncate"
        )
    width = widths.pop()
    index = {residue: i for i, residue in enumerate(alphabet)}
    X = np.zeros((len(rows), width * len(alphabet)), dtype=np.float64)
    for row_at, combo in enumerate(rows):
        for site, residue in enumerate(combo):
            try:
                X[row_at, site * len(alphabet) + index[residue]] = 1.0
            except KeyError:
                raise ValueError(
                    f"residue {residue!r} at site {site} is outside the alphabet "
                    f"{''.join(alphabet)}"
                ) from None
    return X


def encode_landscape(key: str) -> np.ndarray:
    """The encoded design space an arm would be handed for this landscape.

    The benchmark hands its ALDE arm `view.combos` - one residue per library position -
    so this encodes exactly that, from the canonical full-length combos
    `veritas.calibration.landscapes` already builds.
    """
    from veritas.calibration.landscapes import load_landscape

    return one_hot(load_landscape(key).combos)


def constant_columns(X: np.ndarray) -> np.ndarray:
    """Indices where `X.max(0) == X.min(0)`: the columns that normalise to `0/0`.

    botorch's `normalize` is `(X - lower) / (upper - lower)` and ALDE's `domain_discrete`
    takes `lower, upper` as the column min and max, so a column that never varies
    produces NaN, NaN propagates through the first `Linear` layer, the "exclude already
    queried" mask writes NaN over NaN, and `np.argmax` on an all-NaN array returns 0 -
    the same index, forever (`docs/alde-compatibility.md` §3).

    `benchmark_panoptes_isolated.constant_encoding_columns` is the same one-line condition,
    enforced inside a run. This one is deliberately independent of it: the pre-flight gate
    has to answer BEFORE a run and without importing a 4,800-line harness that is under
    active edit. Consolidating the two is a follow-up, and the formula - not an
    implementation - is what both share.
    """
    A = np.asarray(X, dtype=np.float64)
    if A.ndim != 2:
        raise ValueError(f"expected a 2-D design matrix, got shape {A.shape}")
    return np.flatnonzero(A.max(axis=0) == A.min(axis=0))


def check_encoding_compatibility(
    *,
    arm: str,
    X: Optional[np.ndarray] = None,
    library: str = "",
    available: Optional[bool] = None,
    encoded_arms: Optional[Iterable[str]] = None,
) -> CheckResult:
    """`(X.max(0) == X.min(0)).any()` must be false for any arm handed an encoded space.

    **Availability is not compatibility.** `probe_alde` checks the checkout, the
    interpreter, the dependencies and the imports, and correctly reports
    `available: true` on all nine landscapes - including the seven where the arm returns
    one index repeated for its entire budget. The two concepts are therefore separate
    keys in `detail`: `available` is what the probe said, `compatible` is what this
    check measured, and they disagree on seven of nine.
    """
    arms = frozenset(encoded_arms) if encoded_arms is not None else ENCODED_DESIGN_SPACE_ARMS
    check = f"{CHECK_ENCODING}.{arm}"
    if arm not in arms:
        return CheckResult(
            check=check, status=NOT_APPLICABLE, measured=None, expected=None,
            summary=f"{arm} does not normalise an encoded design space by column range",
            reason=(f"the check applies to an arm handed an encoded design space whose "
                    f"pipeline normalises it by column range; {arm!r} is not one of "
                    f"{sorted(arms)} - pass encoded_arms= to extend the set"),
            detail={"available": available, "compatible": None, "library": library},
        )
    if X is None:
        if not library:
            raise ValueError("check_encoding_compatibility needs X= or library=")
        try:
            X = encode_landscape(library)
        except Exception as exc:  # noqa: BLE001 - an unencodable space is a failed check
            return CheckResult(
                check=check, status=FAIL, measured=f"{type(exc).__name__}", expected=0,
                summary=f"the design space for {library} could not be encoded: {exc}",
                detail={"available": available, "compatible": False, "library": library},
            )
    indices = constant_columns(X)
    n_columns = int(np.asarray(X).shape[1])
    n_constant = int(indices.size)
    return CheckResult(
        check=check,
        status=PASS if n_constant == 0 else FAIL,
        measured=n_constant,
        expected=0,
        summary=(f"no constant column in {n_columns} ({library or 'design space'}): every "
                 f"residue in the alphabet is observed at every position"
                 if n_constant == 0 else
                 f"{n_constant} of {n_columns} one-hot columns are constant on "
                 f"{library or 'this design space'} - each is a residue never observed at "
                 f"that position, and normalisation turns it into 0/0 = NaN"),
        detail={
            # Availability and compatibility are deliberately separate.
            "available": available,
            "compatible": n_constant == 0,
            "library": library,
            "n_columns": n_columns,
            "n_constant": n_constant,
            "constant_fraction": n_constant / n_columns if n_columns else 0.0,
            "first_constant_columns": indices[:16].tolist(),
            "reference": REFERENCE_CONSTANT_COLUMNS.get(library),
            "motivating_defect": "docs/alde-compatibility.md",
        },
    )


# ======================================================================================
# 3. Budget
# ======================================================================================
@dataclass(frozen=True)
class BudgetObservation:
    """One (method, seed)'s measurement accounting, as `round_metrics.csv` records it."""

    method: str
    seed: int
    n_unique_measurements: int
    rejected_duplicate: int = 0
    rejected_ineligible: int = 0
    n_distinct_selected: Optional[int] = None
    #: `requested`/`accepted`/`rejected` summed over the campaign's rounds. These separate
    #: "the arm's proposals were refused" from "the arm proposed fewer than the batch
    #: size", which look identical in the distinct-measurement count alone.
    requested: int = 0
    accepted: int = 0
    rejected: int = 0
    #: From the newer `round_metrics.csv` columns, where a run wrote them. These make the
    #: classification exact instead of inferred: `backfill_added` is the brief's
    #: "ineligible and backfilled" (RECOVERABLE) and `degenerate_proposals` is the
    #: benchmark's own determination that the arm's proposals collapsed (DEGENERATE).
    backfill_added: int = 0
    distinct_proposed: Optional[int] = None
    degenerate_proposals: int = 0

    @property
    def has_reason_breakdown(self) -> bool:
        return bool(self.rejected_duplicate or self.rejected_ineligible)

    @property
    def proposed_short_of_the_batch(self) -> bool:
        """The arm offered fewer candidates than were requested, with nothing refused."""
        return bool(self.requested and self.rejected == 0
                    and self.accepted < self.requested)


def classify_budget_shortfall(
    observation: BudgetObservation,
    *,
    declared_budget: int,
    initial_size: int = DEFAULT_INITIAL_SIZE,
) -> Tuple[Optional[str], str]:
    """`(classification, why)`. These are different failures and must not be collapsed.

    * DEGENERATE - the arm could not propose a distinct candidate. HER2 ALDE spent 26 of
      125: the 25 shared starting variants plus one variant measured 100 times. There is
      no budget left to recover and the trajectory is not a search.
    * RECOVERABLE - the arm's proposals were ineligible and the round ran short. GB1
      reliability averaged 123.23 (min 120) and TrpB frontier 120.23 (min 116) this way;
      the arm searched, and spent slightly less than the stated budget.
    * UNCLASSIFIED - a shortfall with no rejection-reason evidence. It still fails,
      because the two cases above require opposite responses.
    """
    shortfall = int(declared_budget) - int(observation.n_unique_measurements)
    if shortfall <= 0:
        return None, "distinct measurements equal the declared budget"
    if observation.degenerate_proposals:
        # The benchmark's own per-round determination, reused rather than second-guessed.
        return BUDGET_DEGENERATE, (
            f"{observation.degenerate_proposals} round(s) flagged degenerate_proposals in "
            f"round_metrics.csv"
            + (f" (distinct_proposed {observation.distinct_proposed} against "
               f"{observation.requested} requested)"
               if observation.distinct_proposed is not None else "")
        )
    if observation.n_distinct_selected is not None and observation.n_distinct_selected <= 1:
        return BUDGET_DEGENERATE, (
            f"the arm selected {observation.n_distinct_selected} distinct variant(s) over "
            f"the whole campaign"
        )
    if observation.n_unique_measurements <= initial_size + 1:
        return BUDGET_DEGENERATE, (
            f"{observation.n_unique_measurements} distinct measurements is the "
            f"{initial_size}-variant starting sample plus "
            f"{observation.n_unique_measurements - initial_size}: the arm proposed nothing "
            f"new it could measure"
        )
    if observation.has_reason_breakdown:
        if observation.rejected_duplicate >= observation.rejected_ineligible:
            return BUDGET_DEGENERATE, (
                f"{observation.rejected_duplicate} rejection(s) were duplicates of variants "
                f"already measured or already in the batch, against "
                f"{observation.rejected_ineligible} ineligible: the arm could not propose a "
                f"distinct candidate"
            )
        return BUDGET_RECOVERABLE, (
            f"{observation.rejected_ineligible} rejection(s) were ineligible proposals "
            f"against {observation.rejected_duplicate} duplicates: the shortfall is "
            f"eligibility-driven"
        )
    if observation.backfill_added:
        # The brief's RECOVERABLE case, stated by the run itself: the arm's proposals were
        # ineligible and the batch was filled from elsewhere.
        return BUDGET_RECOVERABLE, (
            f"{observation.backfill_added} measurement(s) were backfilled after ineligible "
            f"proposals: the shortfall is eligibility-driven, not a collapsed surrogate"
        )
    if observation.proposed_short_of_the_batch:
        # GB1 reliability (mean 123.23, min 120) and TrpB frontier (mean 120.23, min 116)
        # are this case: nothing was refused, the eligible pool simply ran out of
        # candidates the policy would offer. The arm searched.
        return BUDGET_RECOVERABLE, (
            f"the arm offered {observation.accepted} of {observation.requested} requested "
            f"candidates with 0 rejections: the eligible pool ran short, the arm did not "
            f"collapse onto one candidate"
        )
    return BUDGET_UNCLASSIFIED, (
        f"{shortfall} measurement(s) short with no rejection-reason evidence, so "
        f"'could not propose a distinct candidate' cannot be separated from 'proposals "
        f"were ineligible'"
    )


def check_budget(
    observations: Iterable[BudgetObservation],
    *,
    declared_budget: int,
    initial_size: int = DEFAULT_INITIAL_SIZE,
) -> List[CheckResult]:
    """One result per arm: distinct measurements must equal the declared budget.

    `round_metrics.csv` already recorded ALDE's 26 against 125, and already recorded our
    own arms' 120-124 against the baselines' exactly 125. Nothing asserted on either.
    """
    by_method: Dict[str, List[BudgetObservation]] = {}
    for observation in observations:
        by_method.setdefault(observation.method, []).append(observation)

    results: List[CheckResult] = []
    for method, group in sorted(by_method.items()):
        counts = [int(o.n_unique_measurements) for o in group]
        per_seed: List[Dict[str, Any]] = []
        worst: Optional[str] = None
        reasons: List[str] = []
        for observation in sorted(group, key=lambda o: o.seed):
            classification, why = classify_budget_shortfall(
                observation, declared_budget=declared_budget, initial_size=initial_size)
            per_seed.append({
                "seed": observation.seed,
                "n_unique_measurements": observation.n_unique_measurements,
                "classification": classification,
                "why": why,
                "rejected_duplicate": observation.rejected_duplicate,
                "rejected_ineligible": observation.rejected_ineligible,
                "requested": observation.requested,
                "accepted": observation.accepted,
                "rejected": observation.rejected,
            })
            if classification is None:
                continue
            reasons.append(f"seed {observation.seed}: {why}")
            if worst is None or _SEVERITY[classification] > _SEVERITY[worst]:
                worst = classification

        short = [c for c in counts if c < declared_budget]
        missing_evidence = ""
        if worst == BUDGET_UNCLASSIFIED:
            missing_evidence = (
                "no rejection-reason breakdown; selection_trace.csv carries the reason "
                "column that separates a duplicate re-proposal from an ineligible one"
            )
        if worst == BUDGET_DEGENERATE:
            summary = (
                f"{method}: {min(counts)} of {declared_budget} distinct measurements on the "
                f"worst of {len(counts)} seed(s) - DEGENERATE and unrecoverable: the arm "
                f"could not propose a distinct candidate"
            )
        elif worst == BUDGET_RECOVERABLE:
            summary = (
                f"{method}: {min(counts)} of {declared_budget} distinct measurements on the "
                f"worst of {len(counts)} seed(s) - RECOVERABLE: proposals were ineligible "
                f"and the round ran short"
            )
        elif worst == BUDGET_UNCLASSIFIED:
            summary = (
                f"{method}: {min(counts)} of {declared_budget} distinct measurements on the "
                f"worst of {len(counts)} seed(s), and the artefacts do not say which kind "
                f"of shortfall it is"
            )
        else:
            summary = (
                f"{method}: every one of {len(counts)} seed(s) spent all "
                f"{declared_budget} distinct measurements"
            )
        results.append(CheckResult(
            check=f"{CHECK_BUDGET}.{method}",
            status=PASS if worst is None else FAIL,
            measured=min(counts),
            expected=int(declared_budget),
            summary=summary,
            detail={
                "classification": worst,
                "n_seeds": len(counts),
                "min": min(counts),
                "max": max(counts),
                "mean": sum(counts) / len(counts),
                "n_seeds_short": len(short),
                "per_seed": per_seed,
                "reasons": reasons,
                "missing_evidence": missing_evidence,
                "initial_size": initial_size,
            },
        ))
    return results


def budget_observations_from_run(run_dir: Path) -> List[BudgetObservation]:
    """Read `round_metrics.csv` (and `selection_trace.csv` where present) from a run.

    `n_unique_measurements` is taken from each (method, seed)'s LAST round, which is the
    cumulative distinct count the run finished with. The rejection-reason breakdown comes
    from the trace's `reason` column, because `round_metrics.csv` records how many
    selections were rejected but not why - and why is exactly what separates a degenerate
    arm from an eligibility-limited one.
    """
    run_dir = Path(run_dir)
    rounds_path = run_dir / "round_metrics.csv"
    if not rounds_path.is_file():
        raise FileNotFoundError(
            f"no round_metrics.csv in {run_dir}. Distinct-measurement counts are written "
            f"by a run; checked: {rounds_path}"
        )
    last: Dict[Tuple[str, int], Tuple[int, int]] = {}
    totals: Dict[Tuple[str, int], List[int]] = {}
    # The newer degeneracy and fill columns, summed where the run wrote them. A run from
    # before they existed simply has no key here, and the classifier falls back.
    extra: Dict[Tuple[str, int], Dict[str, int]] = {}
    extra_columns = ("backfill_added", "distinct_proposed", "degenerate_proposals")
    with rounds_path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            method = (row.get("method") or "").strip()
            if not method:
                continue
            seed = int(float(row.get("seed") or 0))
            round_index = int(float(row.get("round") or 0))
            n_unique = int(float(row.get("n_unique_measurements") or 0))
            key = (method, seed)
            if key not in last or round_index >= last[key][0]:
                last[key] = (round_index, n_unique)
            running = totals.setdefault(key, [0, 0, 0])
            for at, column in enumerate(("requested", "accepted", "rejected")):
                running[at] += int(float(row.get(column) or 0))
            for column in extra_columns:
                if row.get(column) not in (None, ""):
                    bucket = extra.setdefault(key, {})
                    bucket[column] = bucket.get(column, 0) + int(float(row[column]))

    duplicate: Dict[Tuple[str, int], int] = {}
    ineligible: Dict[Tuple[str, int], int] = {}
    distinct_selected: Dict[Tuple[str, int], Set[str]] = {}
    trace_path = run_dir / "selection_trace.csv"
    if trace_path.is_file():
        with trace_path.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                method = (row.get("method") or "").strip()
                if not method:
                    continue
                key = (method, int(float(row.get("seed") or 0)))
                status = (row.get("status") or "").strip().lower()
                if status == "accepted":
                    distinct_selected.setdefault(key, set()).add(
                        (row.get("variant_id") or "").strip())
                    continue
                if status != "rejected":
                    continue
                reason = (row.get("reason") or "").strip().lower()
                if any(marker in reason for marker in DEGENERACY_REASON_MARKERS):
                    duplicate[key] = duplicate.get(key, 0) + 1
                else:
                    ineligible[key] = ineligible.get(key, 0) + 1

    return [
        BudgetObservation(
            method=method,
            seed=seed,
            n_unique_measurements=n_unique,
            rejected_duplicate=duplicate.get((method, seed), 0),
            rejected_ineligible=ineligible.get((method, seed), 0),
            n_distinct_selected=(len(distinct_selected[(method, seed)])
                                 if (method, seed) in distinct_selected else None),
            requested=totals.get((method, seed), (0, 0, 0))[0],
            accepted=totals.get((method, seed), (0, 0, 0))[1],
            rejected=totals.get((method, seed), (0, 0, 0))[2],
            backfill_added=extra.get((method, seed), {}).get("backfill_added", 0),
            distinct_proposed=extra.get((method, seed), {}).get("distinct_proposed"),
            degenerate_proposals=extra.get((method, seed), {}).get(
                "degenerate_proposals", 0),
        )
        for (method, seed), (_round, n_unique) in sorted(last.items())
    ]


# ======================================================================================
# 4. Endpoint viability
# ======================================================================================
def check_endpoint_viability(
    *,
    library: str,
    raw: int,
    ceiling: int,
    distance: int = DEFAULT_DISTANCE,
    budget: int = DEFAULT_BUDGET,
    tier_size: Optional[int] = None,
    floor: float = INFORMATIVENESS_FLOOR,
) -> CheckResult:
    """Two conditions, both pre-registered, both already published as numbers.

    1. The diversity-adjusted ceiling at the intended `d` must differ materially from the
       raw top-1% count. The rule is `veritas.calibration.endpoint.endpoint_verdict`,
       imported rather than duplicated: DEGENERATE when the ceiling equals the raw count,
       NEAR-DEGENERATE under 5% removed, USABLE otherwise. HER2 is the only DEGENERATE
       target in this corpus - 113 == 113 - so the adjustment is a no-op there and the
       primary endpoint silently becomes the secondary one.
    2. The measurement budget must be smaller than the hit tier. At 125 measurements both
       eqFP611 channels (tier 82) and both CR6261 landscapes (tier 19) are saturated: an
       arm can reach the whole tier, so the endpoint cannot separate two arms.
    """
    tier = int(raw if tier_size is None else tier_size)
    verdict = endpoint_verdict(raw=int(raw), ceiling=int(ceiling), floor=floor)
    saturated = int(budget) >= tier
    violations: List[str] = []
    if not verdict.usable_as_primary:
        violations.append(
            f"the d>={distance} ceiling is {ceiling} against a raw count of {raw} "
            f"({verdict.informativeness:.2%} removed, floor {floor:.0%}): {verdict.verdict}"
        )
    if saturated:
        violations.append(
            f"a budget of {budget} measurements is not smaller than the {tier}-variant hit "
            f"tier: an arm can reach the whole tier, so the endpoint cannot discriminate"
        )
    return CheckResult(
        check=f"{CHECK_ENDPOINT}.{library}" if library else CHECK_ENDPOINT,
        status=PASS if not violations else FAIL,
        measured={"raw": int(raw), "ceiling": int(ceiling),
                  "informativeness": verdict.informativeness, "budget": int(budget),
                  "tier_size": tier},
        expected={"verdict": USABLE,
                  "informativeness": f">= {floor:.0%} of the hit set removed",
                  "budget": f"< the {tier}-variant hit tier"},
        summary=("; ".join(violations) if violations else
                 f"{library or 'target'}: d>={distance} ceiling {ceiling} of {raw} "
                 f"({verdict.informativeness:.1%} removed) and a budget of {budget} under a "
                 f"{tier}-variant tier - usable as the primary endpoint"),
        detail={
            "verdict": verdict.verdict,
            "distance": int(distance),
            "raw": int(raw),
            "ceiling": int(ceiling),
            "informativeness": verdict.informativeness,
            "floor": floor,
            "budget": int(budget),
            "tier_size": tier,
            "budget_saturates_tier": bool(saturated),
            "violations": violations,
            "rule": "veritas.calibration.endpoint.endpoint_verdict",
            "verdict_vocabulary": [DEGENERATE, NEAR_DEGENERATE, USABLE],
        },
    )


# ======================================================================================
# 5. Threshold transfer
# ======================================================================================
#: Where the qualification bar is recorded, best shape first. The first is what a manifest
#: SHOULD carry - a quantile, the absolute it resolved to, and the array it came from. The
#: rest are the shapes real manifests carry today: a bare absolute and nothing else.
QUALIFICATION_MANIFEST_PATHS: Tuple[Tuple[str, ...], ...] = (
    ("qualification",),
    ("configuration", "qualification"),
    ("configuration", "components", "qualify_threshold"),
    ("invocation", "resolved_args", "qualify_threshold"),
)

#: Where `top1.threshold` is recorded, best shape first.
TOP1_MANIFEST_PATHS: Tuple[Tuple[str, ...], ...] = (
    ("top1", "threshold"),
    ("results", "top1", "threshold_value"),
    ("configuration", "top1", "threshold"),
)


def _dig(mapping: Mapping[str, Any], path: Sequence[str]) -> Any:
    node: Any = mapping
    for key in path:
        if not isinstance(node, Mapping) or key not in node:
            return None
        node = node[key]
    return node


def qualification_from_manifest(
    manifest: Mapping[str, Any]
) -> Tuple[Dict[str, Any], str, List[str]]:
    """`(block, read_from, paths_checked)` for the qualification bar.

    A bare number at any of these paths is an absolute constant with no quantile and no
    source, which is exactly the defect: `--qualify-threshold` defaults to 0.70 and a run
    records that number and nothing else. It is normalised into a block so the check
    reports what IS recorded rather than crashing on the shape.

    Every path checked is returned, so a failure never says only "not found".
    """
    checked = [".".join(path) for path in QUALIFICATION_MANIFEST_PATHS]
    for path in QUALIFICATION_MANIFEST_PATHS:
        found = _dig(manifest, path)
        if found is None:
            continue
        where = ".".join(path)
        if isinstance(found, Mapping):
            return dict(found), where, checked
        return {"absolute": float(found), "quantile": None, "source": None}, where, checked
    return {}, "", checked


def top1_threshold_from_manifest(manifest: Mapping[str, Any]) -> Tuple[Optional[float], str]:
    for path in TOP1_MANIFEST_PATHS:
        found = _dig(manifest, path)
        if found is not None and not isinstance(found, Mapping):
            return float(found), ".".join(path)
    return None, ""


def check_threshold_transfer(
    *,
    library: str,
    values: Sequence[float],
    manifest: Mapping[str, Any],
    n_initial: int = DEFAULT_INITIAL_SIZE,
    legacy_absolute: float = 0.70,
    top1_quantile: float = TOP1_QUANTILE,
    tolerance: float = 1e-9,
    gating_p_max: float = GATING_P_MAX,
) -> CheckResult:
    """The qualification bar must be a quantile, resolved from the array top1 uses.

    An absolute constant does not transfer. 0.70 is TrpB's 99.96th percentile - a bar in
    the extreme upper tail - and lies below the global minimum on four landscapes, where
    it qualifies every variant instantly. Every HER2 seed logged "qualified before
    round 1" and the frontier policy's qualification stage gated nothing on 8 of 9
    landscapes.

    So the manifest must record three things, and all three are checked against the data
    rather than believed: the `quantile`, the `absolute` it resolves to on this
    landscape, and the `source` array it was resolved from. The conversion is
    `veritas.calibration.qualification`, reused, not reimplemented.
    """
    v = np.asarray(values, dtype=np.float64)
    block, read_from, paths_checked = qualification_from_manifest(manifest)
    violations: List[str] = []
    detail: Dict[str, Any] = {
        "library": library,
        "n_universe": int(v.size),
        "fitness_min": float(v.min()) if v.size else None,
        "fitness_max": float(v.max()) if v.size else None,
        "recorded": block,
        "read_from": read_from,
        "manifest_paths_checked": paths_checked,
    }
    if not block:
        violations.append(
            "the manifest records no qualification bar at all; checked "
            + ", ".join(paths_checked)
        )

    quantile = block.get("quantile")
    absolute = block.get("absolute")
    source = str(block.get("source") or "").strip()

    if quantile is None:
        violations.append(
            "qualification.quantile is not recorded: the bar is an absolute constant, and "
            "an absolute does not transfer - 0.70 is TrpB's 99.96th percentile but lies "
            "below the global minimum on four landscapes"
        )
    if absolute is None:
        violations.append("qualification.absolute is not recorded")
    if not source:
        violations.append(
            "qualification.source is not recorded, so the array the bar was resolved from "
            "is unknown"
        )

    if quantile is not None:
        resolved = absolute_for_quantile(v, float(quantile))
        detail["resolved_absolute"] = resolved
        if absolute is not None:
            scale = max(1.0, abs(resolved), abs(float(absolute)))
            if abs(float(absolute) - resolved) > tolerance * scale:
                violations.append(
                    f"the recorded absolute {float(absolute)!r} does not resolve from this "
                    f"landscape's values at quantile {float(quantile)!r}, which gives "
                    f"{resolved!r}: the bar was resolved from a different values array"
                )

    recorded_top1, top1_read_from = top1_threshold_from_manifest(manifest)
    detail["top1_read_from"] = top1_read_from
    if recorded_top1 is not None and v.size:
        recomputed_top1 = float(np.percentile(v, top1_quantile))
        detail["top1_threshold_recorded"] = float(recorded_top1)
        detail["top1_threshold_recomputed"] = recomputed_top1
        scale = max(1.0, abs(recomputed_top1), abs(float(recorded_top1)))
        if abs(float(recorded_top1) - recomputed_top1) > tolerance * scale:
            violations.append(
                f"the manifest records top1.threshold {float(recorded_top1)!r} but this "
                f"values array gives {recomputed_top1!r}: the qualification bar and "
                f"top1.threshold are not expressed on one scale, so they cannot be read "
                f"against each other"
            )

    effective = absolute if absolute is not None else detail.get("resolved_absolute")
    if effective is not None and v.size:
        effective = float(effective)
        n_qualifying = qualifying_count(v, effective)
        p_qualify = p_qualify_in_draw(
            n_universe=int(v.size), n_qualifying=n_qualifying, n_draws=n_initial)
        detail.update({
            "effective_absolute": effective,
            "n_qualifying": n_qualifying,
            "qualifying_fraction": n_qualifying / v.size,
            "p_qualify_in_initial_draw": p_qualify,
            "equivalent_quantile_of_legacy_absolute": equivalent_quantile(v, legacy_absolute),
            "legacy_absolute": float(legacy_absolute),
        })
        if effective <= float(v.min()):
            violations.append(
                f"the bar {effective!r} lies at or below the global minimum "
                f"{float(v.min())!r}: it gates nothing, every variant qualifies"
            )
        elif p_qualify >= gating_p_max:
            violations.append(
                f"the {n_initial}-variant initial draw already qualifies with probability "
                f"{p_qualify:.6f} (>= {gating_p_max}): the bar gates nothing. This is the "
                f"HER2 case at the absolute 0.70 default (p = 0.908) - every seed logged "
                f"'qualified before round 1', so the frontier policy's qualification stage "
                f"never ran"
            )

    detail["violations"] = violations
    return CheckResult(
        check=f"{CHECK_THRESHOLD}.{library}" if library else CHECK_THRESHOLD,
        status=PASS if not violations else FAIL,
        measured={"quantile": quantile, "absolute": absolute, "source": source or None},
        expected={"quantile": "recorded, and resolving to the recorded absolute on this "
                              "landscape's own values",
                  "absolute": "recorded",
                  "source": "recorded, and the same array that defines top1.threshold"},
        summary=("; ".join(violations) if violations else
                 f"{library or 'target'}: the bar is quantile {quantile}, resolving to "
                 f"{detail.get('effective_absolute')!r} on this landscape's own values, "
                 f"from the array that defines top1.threshold"),
        detail=detail,
    )


# ======================================================================================
# 6. Row-order invariance
# ======================================================================================
def shuffled_order(n: int, shuffle_seed: int = DEFAULT_SHUFFLE_SEED) -> List[int]:
    """A reproducible permutation of `range(n)`. Same seed, same permutation, always."""
    order = list(range(int(n)))
    random.Random(int(shuffle_seed)).shuffle(order)
    return order


def _fingerprint(value: Any) -> str:
    try:
        payload = json.dumps(_jsonable(value), sort_keys=True, default=str)
    except (TypeError, ValueError):  # pragma: no cover - defensive
        payload = repr(value)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def row_zero_hazard(fitness: Sequence[float]) -> Dict[str, Any]:
    """How much a row-0 index bug would be rewarded on this table.

    Our prepared tables are sorted by fitness descending, so CR9114's row 0 IS the global
    optimum: an arm with an all-NaN surrogate posts a perfect `best_measured_fitness`
    there and beats a working one. The hazard is a property of the table, so it is
    reported even when the invariance run itself cannot be made.
    """
    f = np.asarray(fitness, dtype=np.float64)
    if f.size == 0:
        return {"row_zero_is_global_optimum": False, "percentile": None, "n_rows": 0}
    row_zero = float(f[0])
    maximum = float(f.max())
    percentile = (100.0 * float(np.count_nonzero(f < row_zero)) / (f.size - 1)
                  if f.size > 1 else 100.0)
    return {
        "row_zero_is_global_optimum": bool(row_zero >= maximum),
        "row_zero_fitness": row_zero,
        "max_fitness": maximum,
        "percentile": percentile,
        "n_rows": int(f.size),
    }


def check_row_order_invariance(
    *,
    rows: Sequence[Any],
    runner: Optional[Callable[[Sequence[Any], int], Any]] = None,
    seed: int = 0,
    shuffle_seed: int = DEFAULT_SHUFFLE_SEED,
    fitness: Optional[Sequence[float]] = None,
    library: str = "",
) -> CheckResult:
    """Shuffle the input table's row order, re-run one seed, assert identical results.

    Any correct method is invariant to input row order. ALDE was not: with an all-NaN
    surrogate `np.argmax` returned index 0 on every query of every round, so its result
    was a property of how the corpus happened to be exported. This check catches that
    entire class of index bug, including instances we have not hit.

    Cheap by construction - one seed, one extra run - and the shuffle seed is recorded in
    the result so a failure reproduces exactly. `runner(rows, seed)` may return anything
    JSON-serialisable or `repr`-able; the two returns are compared by fingerprint.

    **Not static.** Without a `runner` the check is NOT-APPLICABLE and says so; the
    row-0 hazard is still reported, because that part is a property of the table.
    """
    check = f"{CHECK_ROW_ORDER}.{library}" if library else CHECK_ROW_ORDER
    detail: Dict[str, Any] = {
        "shuffle_seed": int(shuffle_seed),
        "run_seed": int(seed),
        "n_rows": len(rows),
    }
    if fitness is not None:
        if len(fitness) != len(rows):
            raise ValueError(
                f"fitness has {len(fitness)} entries for {len(rows)} rows: refusing to "
                f"zip unequal-length inputs"
            )
        detail["row_zero_hazard"] = row_zero_hazard(fitness)

    if len(rows) < 2:
        return CheckResult(
            check=check, status=NOT_APPLICABLE, measured=None, expected=None,
            summary="a table of fewer than two rows has only one order",
            reason=f"{len(rows)} row(s): there is no second order to compare against",
            detail=detail,
        )
    if runner is None:
        return CheckResult(
            check=check, status=NOT_APPLICABLE, measured=None, expected=None,
            summary="row-order invariance was not exercised",
            reason=("this check is not static: it needs an actual campaign run, so pass "
                    "runner=(rows, seed) -> result (CLI: --invariance-runner module:function)"),
            detail=detail,
        )

    order = shuffled_order(len(rows), shuffle_seed)
    detail["permutation_head"] = order[:8]
    baseline = runner(list(rows), int(seed))
    shuffled = runner([rows[i] for i in order], int(seed))
    baseline_fingerprint = _fingerprint(baseline)
    shuffled_fingerprint = _fingerprint(shuffled)
    detail["baseline_fingerprint"] = baseline_fingerprint
    detail["shuffled_fingerprint"] = shuffled_fingerprint
    identical = baseline_fingerprint == shuffled_fingerprint
    return CheckResult(
        check=check,
        status=PASS if identical else FAIL,
        measured=shuffled_fingerprint,
        expected=baseline_fingerprint,
        summary=(f"seed {seed} returns identical results under shuffle seed {shuffle_seed}"
                 if identical else
                 f"seed {seed} returns DIFFERENT results when the input table's rows are "
                 f"shuffled (shuffle seed {shuffle_seed}): the method depends on input row "
                 f"order, which is an index bug, not a result"),
        detail=detail,
    )


# ======================================================================================
# The entry point
# ======================================================================================
def _selected(only: Optional[Iterable[str]], skip: Optional[Iterable[str]]) -> Dict[str, str]:
    """`{check: ''}` to run, `{check: reason}` to skip. Unknown names are an error."""
    only_set = set(only) if only is not None else None
    skip_set = set(skip) if skip is not None else set()
    unknown = {name for name in ((only_set or set()) | skip_set) if name not in ALL_CHECKS}
    if unknown:
        raise ValueError(f"unknown check(s) {sorted(unknown)}; known: {list(ALL_CHECKS)}")
    plan: Dict[str, str] = {}
    for check in ALL_CHECKS:
        if only_set is not None and check not in only_set:
            plan[check] = f"not in --only {sorted(only_set)}"
        elif check in skip_set:
            plan[check] = "excluded by --skip"
        else:
            plan[check] = ""
    return plan


def run_preflight(
    *,
    library: str,
    arms: Sequence[str],
    only: Optional[Iterable[str]] = None,
    skip: Optional[Iterable[str]] = None,
    run_dir: Optional[Path] = None,
    manifest: Optional[Mapping[str, Any]] = None,
    budget: int = DEFAULT_BUDGET,
    initial_size: int = DEFAULT_INITIAL_SIZE,
    distance: int = DEFAULT_DISTANCE,
    invariance_runner: Optional[Callable[[Sequence[Any], int], Any]] = None,
    invariance_seed: int = 0,
    shuffle_seed: int = DEFAULT_SHUFFLE_SEED,
    corpus_dir: Optional[Path] = None,
    lock_path: Optional[Path] = None,
    acknowledged_derived: Optional[Mapping[str, str]] = None,
    encoded_arms: Optional[Iterable[str]] = None,
) -> PreflightReport:
    """Run every applicable check for one library and arm set, and return the report.

    Each check is individually skippable through `only`/`skip`, and a skip becomes a
    SKIPPED result carrying its reason rather than a silent absence. A check whose input
    only a run can produce (`budget`, `row_order`) returns NOT-APPLICABLE naming what is
    missing; nothing here reports "fine" for a check it could not make.
    """
    plan = _selected(only, skip)
    arms = tuple(arms)
    results: List[CheckResult] = []
    landscape = None

    def _landscape():
        nonlocal landscape
        if landscape is None:
            from veritas.calibration.landscapes import load_landscape

            landscape = load_landscape(library)
        return landscape

    # ---- 1. corpus -------------------------------------------------------------------
    if plan[CHECK_CORPUS]:
        results.append(_skipped(CHECK_CORPUS, plan[CHECK_CORPUS]))
    else:
        results.extend(check_corpus_integrity(
            corpus_dir=corpus_dir, lock_path=lock_path,
            acknowledged_derived=acknowledged_derived))

    # ---- 2. encoding -----------------------------------------------------------------
    if plan[CHECK_ENCODING]:
        results.append(_skipped(CHECK_ENCODING, plan[CHECK_ENCODING]))
    elif not arms:
        results.append(CheckResult(
            check=CHECK_ENCODING, status=NOT_APPLICABLE, measured=None, expected=None,
            summary="no arm to check", reason="no arms were supplied"))
    else:
        for arm in arms:
            results.append(check_encoding_compatibility(
                arm=arm, library=library, encoded_arms=encoded_arms))

    # ---- 3. budget -------------------------------------------------------------------
    if plan[CHECK_BUDGET]:
        results.append(_skipped(CHECK_BUDGET, plan[CHECK_BUDGET]))
    elif run_dir is None:
        results.append(CheckResult(
            check=CHECK_BUDGET, status=NOT_APPLICABLE, measured=None, expected=None,
            summary="distinct-measurement counts were not available",
            reason=("this check is not static: distinct measurements per (method, seed) "
                    "live in a run's round_metrics.csv, so pass run_dir= (CLI: --run-dir)"),
            detail={"declared_budget": budget}))
    else:
        try:
            observations = budget_observations_from_run(Path(run_dir))
        except FileNotFoundError as exc:
            results.append(CheckResult(
                check=CHECK_BUDGET, status=FAIL, measured=str(exc), expected=budget,
                summary=f"the run directory does not carry a budget accounting: {exc}"))
        else:
            if arms:
                observations = [o for o in observations if o.method in set(arms)]
            if not observations:
                results.append(CheckResult(
                    check=CHECK_BUDGET, status=NOT_APPLICABLE, measured=None, expected=None,
                    summary="no matching (method, seed) rows",
                    reason=f"round_metrics.csv in {run_dir} holds no rows for arms {list(arms)}"))
            else:
                results.extend(check_budget(
                    observations, declared_budget=budget, initial_size=initial_size))

    # ---- 4. endpoint -----------------------------------------------------------------
    if plan[CHECK_ENDPOINT]:
        results.append(_skipped(CHECK_ENDPOINT, plan[CHECK_ENDPOINT]))
    else:
        try:
            ls = _landscape()
            results.append(check_endpoint_viability(
                library=library, raw=ls.n_hits, ceiling=ls.ceilings((distance,))[distance],
                distance=distance, budget=budget))
        except Exception as exc:  # noqa: BLE001 - an unloadable landscape is a failure
            results.append(CheckResult(
                check=f"{CHECK_ENDPOINT}.{library}", status=FAIL,
                measured=f"{type(exc).__name__}: {exc}", expected="a loadable landscape",
                summary=f"the endpoint could not be evaluated for {library}: {exc}"))

    # ---- 5. threshold ----------------------------------------------------------------
    if plan[CHECK_THRESHOLD]:
        results.append(_skipped(CHECK_THRESHOLD, plan[CHECK_THRESHOLD]))
    else:
        resolved_manifest = dict(manifest) if manifest is not None else None
        if resolved_manifest is None and run_dir is not None:
            candidate = Path(run_dir) / "manifest.json"
            if candidate.is_file():
                resolved_manifest = json.loads(candidate.read_text(encoding="utf-8"))
        if resolved_manifest is None:
            results.append(CheckResult(
                check=CHECK_THRESHOLD, status=NOT_APPLICABLE, measured=None, expected=None,
                summary="no manifest to read the qualification bar from",
                reason=("the quantile, the resolved absolute and the source are recorded in "
                        "a run manifest, so pass manifest= or run_dir= (CLI: --manifest)")))
        else:
            try:
                ls = _landscape()
                results.append(check_threshold_transfer(
                    library=library, values=ls.fitness, manifest=resolved_manifest,
                    n_initial=initial_size))
            except Exception as exc:  # noqa: BLE001
                results.append(CheckResult(
                    check=f"{CHECK_THRESHOLD}.{library}", status=FAIL,
                    measured=f"{type(exc).__name__}: {exc}", expected="a loadable landscape",
                    summary=f"the qualification bar could not be checked: {exc}"))

    # ---- 6. row order ----------------------------------------------------------------
    if plan[CHECK_ROW_ORDER]:
        results.append(_skipped(CHECK_ROW_ORDER, plan[CHECK_ROW_ORDER]))
    else:
        try:
            ls = _landscape()
            results.append(check_row_order_invariance(
                rows=ls.combos, fitness=ls.fitness, runner=invariance_runner,
                seed=invariance_seed, shuffle_seed=shuffle_seed, library=library))
        except Exception as exc:  # noqa: BLE001
            results.append(CheckResult(
                check=f"{CHECK_ROW_ORDER}.{library}", status=FAIL,
                measured=f"{type(exc).__name__}: {exc}", expected="a loadable landscape",
                summary=f"row-order invariance could not be checked: {exc}"))

    return PreflightReport(library=library, arms=arms, results=tuple(results))


# ======================================================================================
# CLI
# ======================================================================================
def _load_runner(spec: str) -> Callable[[Sequence[Any], int], Any]:
    """`package.module:function` -> the callable. An explicit spec is a constraint."""
    if ":" not in spec:
        raise ValueError(
            f"--invariance-runner must be 'module:function', got {spec!r}")
    module_name, _, attribute = spec.partition(":")
    module = importlib.import_module(module_name)
    runner = getattr(module, attribute, None)
    if not callable(runner):
        raise ValueError(f"{spec!r} does not name a callable")
    return runner


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m veritas.preflight",
        description=("The pre-flight gate: six assertions over artefacts that already "
                     "recorded a result-invalidating defect. Exits non-zero on failure."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--library", required=True,
                        help="Landscape key, e.g. her2, trpb, cr9114_h1.")
    parser.add_argument("--arms", default="",
                        help="Comma-separated arm names, e.g. alde,ridge,panoptes_default.")
    parser.add_argument("--only", default="",
                        help=f"Run only these checks (comma-separated): {', '.join(ALL_CHECKS)}.")
    parser.add_argument("--skip", default="",
                        help="Skip these checks; the skip is recorded in the report.")
    parser.add_argument("--run-dir", type=Path, default=None,
                        help="A run folder, for round_metrics.csv and manifest.json.")
    parser.add_argument("--manifest", type=Path, default=None,
                        help="A manifest.json, if not the one in --run-dir.")
    parser.add_argument("--budget", type=int, default=DEFAULT_BUDGET,
                        help="Declared distinct-measurement budget (default: %(default)s).")
    parser.add_argument("--initial-size", type=int, default=DEFAULT_INITIAL_SIZE)
    parser.add_argument("--distance", type=int, default=DEFAULT_DISTANCE,
                        help="Intended diversity distance d (default: %(default)s).")
    parser.add_argument("--invariance-runner", default="",
                        help="'module:function' taking (rows, seed); without it the "
                             "row-order check is reported NOT-APPLICABLE.")
    parser.add_argument("--invariance-seed", type=int, default=0)
    parser.add_argument("--shuffle-seed", type=int, default=DEFAULT_SHUFFLE_SEED,
                        help="Recorded in the result so a failure reproduces.")
    parser.add_argument("--corpus-dir", type=Path, default=None)
    parser.add_argument("--lock", type=Path, default=None)
    parser.add_argument("--acknowledge-derived", action="append", default=[],
                        metavar="KEY=WHY",
                        help="Accept a derived table for this landscape, with a reason; "
                             "the acknowledgement is recorded in the report.")
    parser.add_argument("--encoded-arms", default="",
                        help="Extend the set of arms that normalise an encoded design "
                             "space by column range.")
    parser.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    acknowledged: Dict[str, str] = {}
    for item in args.acknowledge_derived:
        key, _, why = str(item).partition("=")
        if not key.strip() or not why.strip():
            print(f"--acknowledge-derived must be KEY=WHY, got {item!r}", file=sys.stderr)
            return 2
        acknowledged[key.strip()] = why.strip()

    manifest = None
    if args.manifest is not None:
        manifest = json.loads(
            paths.explicit_path(args.manifest, what="--manifest").read_text(encoding="utf-8"))

    try:
        report = run_preflight(
            library=args.library,
            arms=tuple(a.strip() for a in args.arms.split(",") if a.strip()),
            only={c.strip() for c in args.only.split(",") if c.strip()} or None,
            skip={c.strip() for c in args.skip.split(",") if c.strip()} or None,
            run_dir=args.run_dir,
            manifest=manifest,
            budget=args.budget,
            initial_size=args.initial_size,
            distance=args.distance,
            invariance_runner=(_load_runner(args.invariance_runner)
                               if args.invariance_runner else None),
            invariance_seed=args.invariance_seed,
            shuffle_seed=args.shuffle_seed,
            corpus_dir=args.corpus_dir,
            lock_path=args.lock,
            acknowledged_derived=acknowledged,
            encoded_arms=({a.strip() for a in args.encoded_arms.split(",") if a.strip()}
                          | set(ENCODED_DESIGN_SPACE_ARMS)) if args.encoded_arms else None,
        )
    except (ValueError, KeyError, ImportError) as exc:
        print(f"pre-flight could not run: {exc}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(report.to_dict(), indent=2))
    else:
        print(report.render())
    return report.exit_code


if __name__ == "__main__":
    raise SystemExit(main())

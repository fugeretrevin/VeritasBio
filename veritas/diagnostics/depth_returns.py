"""
veritas/diagnostics/depth_returns.py - per-round, per-mutation-depth RETURNS.

Why this exists
---------------
The GB1-vs-TrpB reversal has the same mechanism as the TrpB `frontier` collapse:
mutation depth. The planned follow-up intervention is depth-adaptive allocation - a
policy that learns from its own measurements which depth is paying and allocates
later batches toward it. That hypothesis is only testable if per-depth returns were
recorded WHILE the campaign ran, round by round. If the cross-protein benchmark does
not record them, the hypothesis cannot be tested without repeating the whole study.

So this ledger records, per run / method / seed / round / depth:
    1. candidates PROPOSED at that depth
    2. candidates SELECTED at that depth
    3. top-1% HITS RETURNED at that depth
    4. mean and max revealed fitness of the selections at that depth
as tidy long rows, so the round-by-round trajectory can be reconstructed rather than
only the endpoint.

This module does NOT implement depth-adaptive allocation. That is deliberately gated
behind the cross-protein generalisation result: TrpB and GB1 are development data
now, and tuning a policy against them is how one builds something that wins twice and
never again. This is only the measurement that makes it testable later.

The outcome boundary, and why it is a write-only ledger
------------------------------------------------------
Per-depth hit returns are exactly the quantity it would be catastrophic to leak into
acquisition: a policy reading "depth 4 returns 0.9% and depth 2 returns 7.5% on this
landscape" is reading the answer key, not learning from experiments it paid for.

Two different things are mixed in these rows, and they are not equally safe:
  * the revealed fitness of SELECTED variants is a measurement the campaign paid for.
    The selector is already allowed to see it (it is its own evidence).
  * top-1% membership is NOT. The threshold is defined over the whole landscape and
    is evaluator-only, so a per-depth hit COUNT encodes information about variants
    nobody measured.

Because the rows carry both, the ledger is WRITE-ONLY while a run is in progress:
`record_round` returns None, and there is no accessor that reads anything back until
`seal()` has been called. A selector handed a reference to this object can write to
it and learn nothing from it. Sealing is irreversible and blocks all further
recording, so nothing read after sealing can reach a round that has yet to happen.

An online depth-allocation signal - the measured-only half of these rows, scoped to
one campaign - is a separate, reviewed change. It is not exposed here, on purpose.
"""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Hashable, Iterable, List, Optional, Sequence, Tuple

from veritas.diagnostics.genotype import GenotypeParseError, LibrarySpec


class OutcomeDisclosureError(RuntimeError):
    """Raised when outcome-bearing per-depth returns were read before the ledger was sealed."""


@dataclass(frozen=True, slots=True)
class MeasuredSelection:
    """
    One variant the campaign has ALREADY PAID TO MEASURE, with what came back.

    Constructing one is the explicit disclosure point: it must only ever be built
    from a batch the selector has already committed and the evaluator has already
    revealed. `is_top1` is an evaluator-side judgement about a landscape-wide
    threshold, which is why anything holding these objects is write-only until sealed.
    """

    variant_id: str
    revealed_fitness: float
    is_top1: bool

    def __post_init__(self) -> None:
        if not isinstance(self.variant_id, str) or not self.variant_id:
            raise ValueError("MeasuredSelection requires a variant identifier")
        float(self.revealed_fitness)  # raises on a non-numeric outcome


@dataclass(frozen=True, slots=True)
class DepthReturnRow:
    """One (run, method, seed, round, depth) cell. Tidy long format, for analysis."""

    run_id: str
    method: str
    seed: str
    round_index: int
    depth: int
    proposed: int
    selected: int
    top1_hits: int
    mean_selected_fitness: Optional[float]
    max_selected_fitness: Optional[float]
    cumulative_selected: int
    cumulative_top1_hits: int

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def field_names(cls) -> Tuple[str, ...]:
        return tuple(cls.__dataclass_fields__)


@dataclass(frozen=True, slots=True)
class DepthReturnTable:
    """
    The sealed result: every per-depth row, plus the roll-ups the depth-allocation
    hypothesis is judged on.
    """

    run_id: str
    library: LibrarySpec
    rows: Tuple[DepthReturnRow, ...]
    unparseable_selections: int

    # -- roll-ups --------------------------------------------------------------------
    def _filtered(self, method: Optional[str]) -> Tuple[DepthReturnRow, ...]:
        if method is None:
            return self.rows
        return tuple(row for row in self.rows if row.method == method)

    def methods(self) -> Tuple[str, ...]:
        return tuple(sorted({row.method for row in self.rows}))

    def mean_selected_depth(self, method: Optional[str] = None) -> float:
        """Mean mutation depth of the selections - the headline mechanism number."""
        rows = self._filtered(method)
        total = sum(row.selected for row in rows)
        if not total:
            return 0.0
        return sum(row.depth * row.selected for row in rows) / total

    def mean_proposed_depth(self, method: Optional[str] = None) -> float:
        rows = self._filtered(method)
        total = sum(row.proposed for row in rows)
        if not total:
            return 0.0
        return sum(row.depth * row.proposed for row in rows) / total

    def new_hits_per_round(self, method: Optional[str] = None) -> Dict[int, float]:
        """
        Mean NEW top-1% hits per campaign, per round - the per-round return trajectory
        (e.g. reliability on TrpB: 10.44 / 7.04 / 14.27 / 13.31).
        """
        rows = self._filtered(method)
        campaigns = {(row.method, row.seed) for row in rows}
        n_campaigns = max(len(campaigns), 1)
        per_round: Dict[int, int] = {}
        for row in rows:
            per_round[row.round_index] = per_round.get(row.round_index, 0) + row.top1_hits
        return {r: per_round[r] / n_campaigns for r in sorted(per_round)}

    def shell_returns(
        self, method: Optional[str] = None, shallow_max_depth: int = 2
    ) -> Dict[str, Dict[str, float]]:
        """
        Hit rate among the SELECTIONS in the shallow (depth <= k) and deep (depth > k)
        shells: top-1% hits divided by variants selected there.
        """
        rows = self._filtered(method)
        shells: Dict[str, Dict[str, float]] = {}
        for name, predicate in (
            ("shallow", lambda d: d <= shallow_max_depth),
            ("deep", lambda d: d > shallow_max_depth),
        ):
            selected = sum(row.selected for row in rows if predicate(row.depth))
            hits = sum(row.top1_hits for row in rows if predicate(row.depth))
            shells[name] = {
                "selected": float(selected),
                "top1_hits": float(hits),
                "hit_rate": (hits / selected) if selected else 0.0,
            }
        return shells

    def depth_profile(self, method: Optional[str] = None) -> Dict[int, Dict[str, float]]:
        rows = self._filtered(method)
        profile: Dict[int, Dict[str, float]] = {}
        for row in rows:
            cell = profile.setdefault(
                row.depth, {"proposed": 0.0, "selected": 0.0, "top1_hits": 0.0, "max": float("-inf")}
            )
            cell["proposed"] += row.proposed
            cell["selected"] += row.selected
            cell["top1_hits"] += row.top1_hits
            if row.max_selected_fitness is not None:
                cell["max"] = max(cell["max"], row.max_selected_fitness)
        for cell in profile.values():
            cell["hit_rate"] = (cell["top1_hits"] / cell["selected"]) if cell["selected"] else 0.0
            if cell["max"] == float("-inf"):
                cell["max"] = 0.0
        return dict(sorted(profile.items()))

    def to_dict(self) -> Dict[str, Any]:
        return {
            "run_id": self.run_id,
            "library": dict(self.library.to_dict()),
            "unparseable_selections": self.unparseable_selections,
            "methods": list(self.methods()),
            "mean_selected_depth": {m: self.mean_selected_depth(m) for m in self.methods()},
            "mean_proposed_depth": {m: self.mean_proposed_depth(m) for m in self.methods()},
            "new_hits_per_round": {
                m: {str(k): v for k, v in self.new_hits_per_round(m).items()}
                for m in self.methods()
            },
            "shell_returns": {m: self.shell_returns(m) for m in self.methods()},
            "rows": [row.to_dict() for row in self.rows],
        }


class DepthReturnLedger:
    """
    Write-only, per-round, per-depth return ledger.

    Live usage, after the evaluator has revealed a committed batch:

        ledger.record_round(
            method=arm, seed=seed, round_index=r,
            proposed=[...ids...],
            measured=[MeasuredSelection(id, fitness, is_top1), ...],
        )

    and once, when the run is over:

        table = ledger.seal()

    Nothing can be read back before `seal()`. That is the enforcement, not a comment:
    an acquisition path holding this object can contribute to it and learn nothing
    from it.
    """

    __slots__ = (
        "_library",
        "_run_id",
        "_cells",
        "_depth_cache",
        "_unparseable",
        "_sealed",
    )

    def __init__(self, library: LibrarySpec, run_id: str = "") -> None:
        if not isinstance(library, LibrarySpec):
            raise TypeError(f"library must be a LibrarySpec, got {type(library).__name__}")
        self._library = library
        self._run_id = str(run_id)
        # (method, seed, round_index, depth) -> tallies
        self._cells: Dict[Tuple[str, str, int, int], Dict[str, Any]] = {}
        self._depth_cache: Dict[str, int] = {}
        self._unparseable = 0
        self._sealed = False

    @property
    def sealed(self) -> bool:
        return self._sealed

    # -- the write-only guarantee ----------------------------------------------------
    @property
    def rows(self) -> Tuple[DepthReturnRow, ...]:
        raise OutcomeDisclosureError(
            "DepthReturnLedger is write-only while a run is in progress. Its rows carry "
            "top-1% membership, which is defined over the whole landscape and is "
            "evaluator-only, so they are readable only from the table returned by "
            "seal(). An online depth-allocation signal must be built from measured "
            "outcomes alone and is a separate, reviewed change."
        )

    def _depth(self, variant_id: str) -> int:
        cached = self._depth_cache.get(variant_id)
        if cached is None:
            cached = self._library.mutation_depth(variant_id)
            self._depth_cache[variant_id] = cached
        return cached

    def record_round(
        self,
        *,
        method: str,
        seed: Hashable,
        round_index: int,
        proposed: Sequence[str] = (),
        measured: Sequence[MeasuredSelection] = (),
    ) -> None:
        """
        Records one round. Returns None - by design, this call discloses nothing.

        `proposed` is the ballot (identifiers only); `measured` is the committed batch
        AFTER the evaluator revealed it.
        """
        if self._sealed:
            raise RuntimeError(
                "this DepthReturnLedger is sealed and can record no further round"
            )
        key_method, key_seed = str(method), str(seed)

        for variant_id in proposed or ():
            if not isinstance(variant_id, str):
                raise TypeError(
                    "proposed must contain bare variant-identifier strings, got "
                    f"{type(variant_id).__name__}"
                )
            try:
                depth = self._depth(variant_id)
            except GenotypeParseError:
                self._unparseable += 1
                continue
            self._cell(key_method, key_seed, round_index, depth)["proposed"] += 1

        for selection in measured or ():
            if not isinstance(selection, MeasuredSelection):
                raise TypeError(
                    "measured must contain MeasuredSelection records - constructing one "
                    "is the explicit point at which an already-paid-for outcome is "
                    f"disclosed. Got {type(selection).__name__}."
                )
            try:
                depth = self._depth(selection.variant_id)
            except GenotypeParseError:
                self._unparseable += 1
                continue
            cell = self._cell(key_method, key_seed, round_index, depth)
            cell["selected"] += 1
            cell["hits"] += 1 if selection.is_top1 else 0
            cell["values"].append(float(selection.revealed_fitness))

    def _cell(self, method: str, seed: str, round_index: int, depth: int) -> Dict[str, Any]:
        key = (method, seed, int(round_index), int(depth))
        cell = self._cells.get(key)
        if cell is None:
            cell = {"proposed": 0, "selected": 0, "hits": 0, "values": []}
            self._cells[key] = cell
        return cell

    def seal(self) -> DepthReturnTable:
        """Freezes the ledger and returns the readable table."""
        self._sealed = True
        rows: List[DepthReturnRow] = []
        running: Dict[Tuple[str, str, int], Tuple[int, int]] = {}
        for key in sorted(self._cells):
            method, seed, round_index, depth = key
            cell = self._cells[key]
            selected_before, hits_before = running.get((method, seed, depth), (0, 0))
            cumulative_selected = selected_before + cell["selected"]
            cumulative_hits = hits_before + cell["hits"]
            running[(method, seed, depth)] = (cumulative_selected, cumulative_hits)
            values: List[float] = cell["values"]
            rows.append(
                DepthReturnRow(
                    run_id=self._run_id,
                    method=method,
                    seed=seed,
                    round_index=round_index,
                    depth=depth,
                    proposed=cell["proposed"],
                    selected=cell["selected"],
                    top1_hits=cell["hits"],
                    mean_selected_fitness=(sum(values) / len(values)) if values else None,
                    max_selected_fitness=max(values) if values else None,
                    cumulative_selected=cumulative_selected,
                    cumulative_top1_hits=cumulative_hits,
                )
            )
        return DepthReturnTable(
            run_id=self._run_id,
            library=self._library,
            rows=tuple(rows),
            unparseable_selections=self._unparseable,
        )


def write_depth_returns_csv(
    path: Path,
    rows: Iterable[DepthReturnRow],
    stamp: str,
) -> int:
    """Tidy long format, one row per run/method/seed/round/depth. Returns the row count."""
    header = ("provenance",) + DepthReturnRow.field_names()
    written = 0
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow(
                (stamp,)
                + tuple(
                    "" if getattr(row, name) is None else getattr(row, name)
                    for name in DepthReturnRow.field_names()
                )
            )
            written += 1
    return written

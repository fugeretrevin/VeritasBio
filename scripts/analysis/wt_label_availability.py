"""Per-landscape: is the wild type's MEASURED label available, and what is it?

Why this exists
---------------
`workers/math_worker.py::select_next_micro_batch` injects one synthetic training row per
round for the wild type carrying the literal label `y = 0.0` at variance `1e-4`. Ablating
that row is worth +29.75 top-1% hits on CreiLOV (t=12.17, 12/12 seeds). Ablation deletes
the anchor; the principled repair is to anchor at the wild type's MEASURED label. That
repair is only possible where such a label exists, so this script establishes, per
landscape and before any fix is written, whether one does.

Three independent sources are consulted and cross-checked, because a single source would
be exactly the "check validated only against its own fixture" failure of HANDOFF rule
2.10:

  A. The PREPARED table (`scripts/data/*.csv`). The repo convention drops the parent -
     "the parent is not a library member" - so a hit here is the exception, not the rule.
     GB1 and HER2 keep a literal `WT` row; nobody else does.

  B. The UPSTREAM CombinGym workbook's `n_mut == 0` row, for the six CombinGym
     landscapes. Usable only once the upstream label column is PROVEN to be the same
     quantity as our prepared `DMS_score`: the script joins the two sorted label vectors
     over the mutants and refuses the value unless they agree to float precision.

  C. The prepared table's own `DMS_score_bin` column. `docs/dataset-inventory.md` §7
     records the rule `DMS_score_bin = 1 iff DMS_score >= the wild-type value`, so the
     column BRACKETS the wild-type label without reading any external file:
     max{score : bin == 0} < wt_label <= min{score : bin == 1}. An A/B value outside that
     bracket is a provenance error and is reported as one.

What the numbers mean
---------------------
Two different distances get confused, and the difference decides which landscapes are
negative controls:

  * `zero_sd_from_centre`  - |0 - mean| / sd. How much of an OUTLIER the fabricated row
    is among the labels. This is the 15.68 / 14.14 / 0.35 / 0.20 ordering in
    docs/HANDOFF.md.
  * `anchor_error_sd`      - |wt_label - 0| / sd. How WRONG the fabricated label is about
    the point it claims to pin. This is the quantity the repair corrects, and it is NOT
    the same ordering: TrpB's parent measures 0.408 on a label sd of 0.0599.

Run:
  .venv\\Scripts\\python.exe -X utf8 scripts/analysis/wt_label_availability.py
  .venv\\Scripts\\python.exe -X utf8 scripts/analysis/wt_label_availability.py --json
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from veritas.calibration.genotypes import canonical_combo, depth  # noqa: E402
from veritas.calibration.landscapes import REGISTRY, _pick, load_landscape  # noqa: E402
from veritas.paths import repository_root  # noqa: E402

#: Where the pinned CombinGym corpus lives. Gitignored, so it resolves against the REAL
#: repository root and never against a worktree - HANDOFF rule 2.8.
CORPUS_DIR = repository_root() / "downloads" / "combingym"

#: Upstream workbook, its parent-label column, and the column our prepared CSV was cut
#: from. `None` means this landscape has no CombinGym upstream.
UPSTREAM: Dict[str, Tuple[str, str]] = {
    "creilov": ("Data/DMS/Clean/CreiLOV_clean.xlsx", "mean_log"),
    "eqfp611_blue": ("Data/DMS/Clean/eqFP611_clean.xlsx", "blue"),
    "eqfp611_red": ("Data/DMS/Clean/eqFP611_clean.xlsx", "red"),
    "cr6261_h1": ("Data/DMS/Clean/bnAbs_CR6261_H1_clean.xlsx", "h1_mean"),
    "cr6261_h9": ("Data/DMS/Clean/bnAbs_CR6261_H9_clean.xlsx", "h9_mean"),
    "gb1": ("Data/DMS/Clean/GB1_clean.xlsx", "fitness"),
    # CR9114's prepared table was re-cut from the RAW 65,536-genotype file
    # (docs/cr9114-reprep.md), so the clean workbook is a 48,841-row SUBSET and its
    # label vector cannot join against ours. The germline value is pinned in-repo
    # instead; see D below.
}

#: Values already committed to this repository as the parent's measured fitness. Each is
#: an independent in-repo record, not a value this script derives.
IN_REPO: Dict[str, Tuple[float, str]] = {
    "cr9114_h1": (
        8.424622015254,
        "scripts/data/prepare_cr9114_h1.py::CR9114_GERMLINE_SCORE",
    ),
    "trpb": (
        0.408073925,
        "workers/data/TRPB_TM9D8.json::veritasbio_provenance.parent_measured_fitness",
    ),
}


@dataclass
class Row:
    key: str
    n_universe: int = 0
    label_min: float = math.nan
    label_max: float = math.nan
    label_mean: float = math.nan
    label_sd: float = math.nan
    threshold: float = math.nan
    #: A. prepared-table parent row
    prepared_wt_label: Optional[float] = None
    prepared_wt_rows: int = 0
    #: B. upstream n_mut == 0 row, with the join that licenses it
    upstream_wt_label: Optional[float] = None
    upstream_join: Optional[str] = None
    #: C. DMS_score_bin bracket
    bracket: Optional[Tuple[Optional[float], Optional[float]]] = None
    #: D. in-repo committed constant
    in_repo_label: Optional[float] = None
    in_repo_source: Optional[str] = None
    #: resolution
    wt_label: Optional[float] = None
    wt_source: str = "ABSENT"
    notes: List[str] = field(default_factory=list)

    # -- derived -------------------------------------------------------------------
    @property
    def zero_sd_from_centre(self) -> float:
        return abs(0.0 - self.label_mean) / self.label_sd if self.label_sd else math.nan

    @property
    def anchor_error_sd(self) -> Optional[float]:
        if self.wt_label is None or not self.label_sd:
            return None
        return abs(self.wt_label - 0.0) / self.label_sd

    @property
    def wt_percentile(self) -> Optional[float]:
        return self._pct

    _pct: Optional[float] = None


def _prepared_parent(key: str) -> Tuple[Optional[float], int]:
    """The parent row's label in our own prepared CSV, if the table kept one.

    Rule 2.1: depth comes from `canonical_combo` + `depth`, which fill wild type at every
    omitted position by construction. Nothing is zipped and no string is compared
    positionally.
    """
    spec = REGISTRY[key]
    labels: List[float] = []
    with spec.path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        id_col = _pick(cols, ("mutant",))
        value_col = _pick(cols, ("dms_score", "fitness", "score"))
        stage_col = _pick(cols, ("stage",))
        for row in reader:
            if spec.labelled_stage is not None and stage_col:
                if (row.get(stage_col) or "").strip() != spec.labelled_stage:
                    continue
            vid = (row.get(id_col) or "").strip()
            combo = canonical_combo(vid, spec.positions, spec.parent_residues)
            if combo is None:
                continue
            if depth(combo, spec.parent_residues) != 0:
                continue
            try:
                labels.append(float(row[value_col]))
            except (TypeError, ValueError, KeyError):
                continue
    if not labels:
        return None, 0
    return labels[0], len(labels)


def _bin_bracket(key: str) -> Tuple[Optional[float], Optional[float]]:
    """(max label with bin==0, min label with bin==1) - the wild type lies in between.

    `docs/dataset-inventory.md` §7: `DMS_score_bin = 1 iff DMS_score >= wild-type value`.
    Returns (None, None) where the column is absent or constant.
    """
    spec = REGISTRY[key]
    lo: Optional[float] = None
    hi: Optional[float] = None
    with spec.path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        value_col = _pick(cols, ("dms_score", "fitness", "score"))
        bin_col = _pick(cols, ("dms_score_bin",))
        stage_col = _pick(cols, ("stage",))
        id_col = _pick(cols, ("mutant",))
        if bin_col is None:
            return None, None
        for row in reader:
            if spec.labelled_stage is not None and stage_col:
                if (row.get(stage_col) or "").strip() != spec.labelled_stage:
                    continue
            # Exclude the parent's own row: it is on both sides of its own cutoff.
            vid = (row.get(id_col) or "").strip()
            combo = canonical_combo(vid, spec.positions, spec.parent_residues)
            if combo is not None and depth(combo, spec.parent_residues) == 0:
                continue
            try:
                v = float(row[value_col])
                b = int(float(row[bin_col]))
            except (TypeError, ValueError, KeyError):
                continue
            if b == 0:
                lo = v if lo is None else max(lo, v)
            else:
                hi = v if hi is None else min(hi, v)
    return lo, hi


def _upstream_parent(key: str) -> Tuple[Optional[float], Optional[str]]:
    """The upstream `n_mut == 0` label, and the join that proves the column is ours.

    Refuses the value unless the upstream mutant labels and our prepared labels are the
    same multiset to float precision. A column that is merely plausibly named is not
    evidence; HANDOFF rule 2.10.
    """
    entry = UPSTREAM.get(key)
    if entry is None:
        return None, None
    rel, col = entry
    path = CORPUS_DIR / rel
    if not path.is_file():
        return None, f"UPSTREAM MISSING: checked {path}"
    try:
        import pandas as pd
    except ImportError:  # pragma: no cover
        return None, "pandas unavailable; cannot read the upstream workbook"

    df = pd.read_excel(path)
    if "n_mut" not in df.columns or col not in df.columns:
        return None, f"{path.name}: no n_mut/{col} column (has {list(df.columns)})"
    parent = df[df["n_mut"] == 0]
    if len(parent) != 1:
        return None, f"{path.name}: {len(parent)} rows at n_mut == 0, expected exactly 1"

    up = df[df["n_mut"] > 0][col].dropna().to_numpy(dtype=float)
    ours = load_landscape(key).fitness
    if len(up) != len(ours):
        return None, (f"JOIN REFUSED: upstream has {len(up)} mutant labels, the prepared "
                      f"table has {len(ours)}")
    diff = float(np.abs(np.sort(up) - np.sort(ours)).max())
    if diff > 1e-9:
        return None, f"JOIN REFUSED: max |upstream - prepared| = {diff:.3e}"
    return float(parent[col].iloc[0]), f"join exact (max abs diff {diff:.1e}, n={len(up)})"


def audit(key: str) -> Row:
    r = Row(key=key)
    ls = load_landscape(key)
    f = ls.fitness
    r.n_universe = int(ls.n_universe)
    r.label_min, r.label_max = float(f.min()), float(f.max())
    r.label_mean, r.label_sd = float(f.mean()), float(f.std(ddof=0))
    r.threshold = float(ls.threshold)

    r.prepared_wt_label, r.prepared_wt_rows = _prepared_parent(key)
    r.upstream_wt_label, r.upstream_join = _upstream_parent(key)
    r.bracket = _bin_bracket(key)
    if key in IN_REPO:
        r.in_repo_label, r.in_repo_source = IN_REPO[key]

    # Resolution order: our own table, then a join-verified upstream, then an in-repo
    # committed constant. Every source that produced a value must AGREE with the others.
    candidates = [
        ("prepared_table", r.prepared_wt_label),
        ("upstream_combingym", r.upstream_wt_label),
        ("in_repo_constant", r.in_repo_label),
    ]
    found = [(s, v) for s, v in candidates if v is not None]
    if found:
        r.wt_source = "+".join(s for s, _ in found)
        r.wt_label = found[0][1]
        vals = [v for _, v in found]
        if max(vals) - min(vals) > 1e-6:
            r.notes.append(
                f"SOURCES DISAGREE: {dict(found)} - the label is NOT resolved")
            r.wt_label = None
            r.wt_source = "CONFLICTED"
    if r.upstream_join and r.upstream_wt_label is None:
        r.notes.append(r.upstream_join)

    # Bracket cross-check (source C against A/B/D).
    lo, hi = r.bracket
    if r.wt_label is not None and (lo is not None or hi is not None):
        ok_lo = lo is None or r.wt_label > lo
        ok_hi = hi is None or r.wt_label <= hi
        if ok_lo and ok_hi:
            r.notes.append(f"bracket OK: {lo!r} < {r.wt_label:.9g} <= {hi!r}")
        else:
            r.notes.append(
                f"BRACKET VIOLATED: DMS_score_bin implies ({lo!r}, {hi!r}] but the "
                f"resolved label is {r.wt_label:.9g}")
    if r.wt_label is not None:
        r._pct = float((f < r.wt_label).mean() * 100.0)
    return r


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    args = ap.parse_args()

    rows = [audit(k) for k in REGISTRY]

    if args.json:
        out = []
        for r in rows:
            d = {k: v for k, v in r.__dict__.items() if not k.startswith("_")}
            d["zero_sd_from_centre"] = r.zero_sd_from_centre
            d["anchor_error_sd"] = r.anchor_error_sd
            d["wt_percentile"] = r.wt_percentile
            d["wt_is_hit"] = (None if r.wt_label is None
                              else bool(r.wt_label >= r.threshold))
            out.append(d)
        json.dump(out, sys.stdout, indent=2, default=str)
        print()
        return 0

    print("WILD-TYPE MEASURED LABEL AVAILABILITY, per landscape")
    print(f"corpus dir checked: {CORPUS_DIR}  (exists={CORPUS_DIR.is_dir()})")
    print()
    head = (f"{'landscape':14s} {'in prep?':9s} {'WT label':>14s} {'source':34s} "
            f"{'0-SD':>7s} {'err-SD':>7s} {'pct':>6s} {'hit?':>5s}")
    print(head)
    print("-" * len(head))
    for r in rows:
        prep = "YES" if r.prepared_wt_rows else "no"
        lab = "ABSENT" if r.wt_label is None else f"{r.wt_label:.9g}"
        err = "-" if r.anchor_error_sd is None else f"{r.anchor_error_sd:7.2f}"
        pct = "-" if r.wt_percentile is None else f"{r.wt_percentile:6.2f}"
        hit = "-" if r.wt_label is None else ("HIT" if r.wt_label >= r.threshold else "no")
        print(f"{r.key:14s} {prep:9s} {lab:>14s} {r.wt_source:34s} "
              f"{r.zero_sd_from_centre:7.2f} {err:>7s} {pct:>6s} {hit:>5s}")
    print()
    print("0-SD   = |0 - label mean| / label sd      how much of an outlier y=0.0 is")
    print("err-SD = |wt_label - 0| / label sd        how WRONG y=0.0 is about the parent")
    print("pct    = percentile of the wild type within the eligible universe")
    print()
    for r in rows:
        print(f"{r.key}:")
        print(f"  labels [{r.label_min:.6g}, {r.label_max:.6g}] mean {r.label_mean:.6g} "
              f"sd {r.label_sd:.6g} threshold {r.threshold:.6g} n={r.n_universe}")
        if r.upstream_join:
            print(f"  upstream: {r.upstream_join}")
        if r.in_repo_source:
            print(f"  in-repo:  {r.in_repo_label!r} from {r.in_repo_source}")
        for n in r.notes:
            print(f"  note: {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

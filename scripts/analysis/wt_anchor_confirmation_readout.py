"""Paired read-out for the wt_anchor confirmation campaign.

Committed because HANDOFF section 1's criticism of the GB1 "before" figure was that it
was "a number quoted in a document, not a re-analysable artifact". Every figure in
docs/wt-anchor-confirmation-results.md comes out of this script.

RENAMED at the integration merge, 2026-11-01. This file was `wt_anchor_readout.py` on
branch `confirm-wt-anchor-48` and collided add/add with a DIFFERENT script of the same
name from `wt-anchor-measured-label`. They are not two versions of one tool: this one
resolves run directories BY NAME from a campaign log and reports the two-arm ablation
confirmation; the other takes a results root and reports the THREE anchor modes on both
endpoints. Each is the provenance for its own results document, so both are kept and
this one took the longer name.

Usage, from the checkout:

    <real-root>\\.venv\\Scripts\\python.exe -X utf8 \\
        scripts/analysis/wt_anchor_confirmation_readout.py \\
        wtanchor_20261101_103453.log

Run directories are resolved BY NAME from the campaign log the .ps1 wrote, never from
"the newest directory" - during a period of concurrent sweeps the newest directory
belonged to another campaign and TrpB was reported as CreiLOV for an hour
(docs/day1-reachability.md section 11, rule 3).

What it enforces before it will report anything:

  * a run with no manifest.json is a PARTIAL LAUNCH and not data - the shared results
    root carries orphans left by task-stops, which are not process kills;
  * every seed must have spent the full budget with status `completed`, or the arms are
    not budget-matched and the landscape is REFUSED;
  * `ridge` is independent of --ablate-components, so its full per-round sub-table must
    be bit-identical across the two arms, with matching seed_block_sha256 and top1
    threshold. This is a refusal, not a warning;
  * configuration.components.states.wt_anchor must read True in the default arm and
    False in the ablated one, so the manifest is confirmed to describe the arm that ran
    (rule 2.12: a record that reads one of two gates is worse than no record).

Statistics are paired over seeds (rule 2.11): a difference of means is not a verdict.
Every non-significant result is printed beside the minimum effect the design could have
detected, so a null cannot be mistaken for a shortage of power.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

#: initial + 4 x batch, from wt_anchor_confirmation.ps1's registered plan.
BUDGETS = {"creilov": 125, "trpb": 125, "gb1": 125, "cr9114_h1": 125,
           "eqfp611_red": 81, "eqfp611_blue": 81}

#: Wall-clock columns must differ between two runs; everything describing the RESULT
#: must not.
_VOLATILE = {"round_seconds", "cumulative_seconds", "token_usage"}

_DONE = re.compile(r"^## (\S+) \[(\w+)\] done .* run_dir=(run_\d{8}_\d{6}_[0-9a-f]+)")


def results_dir() -> Path:
    """The checkout's own gitignored results root.

    The benchmark writes run folders beside the sources it ran, so from a worktree they
    land in the worktree - while the interpreter and the corpus come from the real root.
    Rule 2.8's point is to keep those questions apart.
    """
    return Path(__file__).resolve().parent.parent.parent / "benchmark_comparison_results"


def parse_log(log_path: Path) -> dict[str, dict[str, str]]:
    if not log_path.exists():
        raise SystemExit(f"no campaign log. Checked: {log_path}")
    out: dict[str, dict[str, str]] = {}
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = _DONE.match(line.strip())
        if m:
            out.setdefault(m.group(1), {})[m.group(2)] = m.group(3)
    return out


def _round_rows(run_dir: Path) -> list[dict[str, str]]:
    path = run_dir / "round_metrics.csv"
    if not path.exists():
        raise SystemExit(f"no round_metrics.csv. Checked: {path}")
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def final_round_hits(run_dir: Path) -> dict[str, dict[int, float]]:
    """{method: {seed: cum_top1_hits_excluding_initial at that seed's LAST round}}"""
    best: dict[str, dict[int, tuple[int, float]]] = {}
    for row in _round_rows(run_dir):
        m, s, r = row["method"], int(row["seed"]), int(row["round"])
        v = float(row["cum_top1_hits_excluding_initial"])
        cur = best.setdefault(m, {}).get(s)
        if cur is None or r > cur[0]:
            best[m][s] = (r, v)
    return {m: {s: v for s, (_r, v) in d.items()} for m, d in best.items()}


def ridge_rows(run_dir: Path) -> dict[tuple[int, int], tuple]:
    rows: dict[tuple[int, int], tuple] = {}
    with (run_dir / "round_metrics.csv").open(newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        cols = [c for c in (rd.fieldnames or []) if c not in _VOLATILE]
        for row in rd:
            if row["method"] == "ridge":
                rows[(int(row["seed"]), int(row["round"]))] = tuple(row[c] for c in cols)
    return rows


def manifest(run_dir: Path) -> dict:
    mp = run_dir / "manifest.json"
    if not mp.exists():
        raise SystemExit(f"no manifest.json. Checked: {mp}")
    man = json.loads(mp.read_text(encoding="utf-8"))

    def dig(*keys):
        cur = man
        for k in keys:
            if not isinstance(cur, dict) or k not in cur:
                return None
            cur = cur[k]
        return cur

    return {
        "seed_block": dig("configuration", "seeding", "seed_block_sha256"),
        "fingerprint": dig("configuration_fingerprint"),
        "disabled": dig("configuration", "acquisition_policy", "components_disabled"),
        "wt_anchor_state": dig("configuration", "components", "states", "wt_anchor"),
        "threshold": dig("results", "top1", "threshold_value"),
        "library": dig("configuration", "dataset", "library", "key"),
    }


def integrity(run_dir: Path, budget: int | None, seeds: int) -> tuple[bool, list[str]]:
    if not run_dir.is_dir():
        return False, [f"run directory absent. Checked: {run_dir}"]
    if not (run_dir / "manifest.json").exists():
        return False, [f"NO MANIFEST - partial launch, not data. "
                       f"Checked: {run_dir / 'manifest.json'}"]
    if budget is None:
        return False, [f"no registered budget for this landscape; add it to BUDGETS"]
    rows = _round_rows(run_dir)
    notes: list[str] = []
    ok = True
    for meth in sorted({r["method"] for r in rows}):
        per: dict[int, tuple[int, str, int]] = {}
        for r in (x for x in rows if x["method"] == meth):
            s, rd = int(r["seed"]), int(r["round"])
            if s not in per or rd > per[s][0]:
                per[s] = (rd, r["status"], int(r["n_unique_measurements"]))
        problems = []
        if len(per) != seeds:
            problems.append(f"{len(per)} seeds, expected {seeds}")
        short = {s: v[2] for s, v in per.items() if v[2] != budget}
        if short:
            problems.append(f"{len(short)} seed(s) short of budget {budget}: "
                            f"{dict(list(short.items())[:5])}")
        bad = {s: v[1] for s, v in per.items() if v[1] != "completed"}
        if bad:
            problems.append(f"{len(bad)} not 'completed': {dict(list(bad.items())[:5])}")
        if problems:
            ok = False
            notes.extend(f"{meth}: {p}" for p in problems)
        else:
            notes.append(f"{meth}: {len(per)}/{seeds} seeds, all at budget {budget}, "
                         f"all completed")
    return ok, notes


def paired(a: dict[int, float], b: dict[int, float], n_boot: int = 10000, rs: int = 12345):
    """a - b over the seeds present in both, with a whole-seed bootstrap."""
    seeds = sorted(set(a) & set(b))
    d = np.array([a[s] - b[s] for s in seeds], dtype=float)
    n = len(d)
    mean = float(d.mean())
    sd = float(d.std(ddof=1)) if n > 1 else 0.0
    if sd > 0:
        t = mean / (sd / math.sqrt(n))
    else:
        t = math.copysign(math.inf, mean) if mean else 0.0
    rng = np.random.default_rng(rs)
    boots = d[rng.integers(0, n, size=(n_boot, n))].mean(axis=1)
    lo, hi = (float(x) for x in np.percentile(boots, [2.5, 97.5]))
    return {"n": n, "mean": mean, "sd": sd, "t": t, "lo": lo, "hi": hi,
            "mde": 2.0 * sd / math.sqrt(n) if n else float("nan"),
            "better": int((d > 0).sum()), "worse": int((d < 0).sum()),
            "tied": int((d == 0).sum())}


def line(label: str, r: dict) -> str:
    tt = f"{r['t']:+7.2f}" if math.isfinite(r["t"]) else ("    inf" if r["t"] > 0 else "   -inf")
    return (f"  {label:<40} d={r['mean']:+8.3f}  t={tt}  "
            f"B/W/T={r['better']:2d}/{r['worse']:2d}/{r['tied']:2d}  "
            f"boot95=[{r['lo']:+8.3f},{r['hi']:+8.3f}]  MDE@t2={r['mde']:6.3f}  n={r['n']}")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("log", help="campaign log file name inside benchmark_comparison_results")
    ap.add_argument("--seeds", type=int, default=48)
    args = ap.parse_args(argv)

    root = results_dir()
    plan = parse_log(root / args.log)
    if not plan:
        raise SystemExit(f"no completed arms parsed from {root / args.log}")

    rc = 0
    for lib, arms in plan.items():
        print("=" * 100)
        print(f"LANDSCAPE: {lib}")
        if set(arms) != {"default", "ablated"}:
            print(f"  INCOMPLETE: arms present = {sorted(arms)} - REFUSED")
            rc = 1
            continue

        integ_ok = True
        for tag in ("default", "ablated"):
            good, notes = integrity(root / arms[tag], BUDGETS.get(lib), args.seeds)
            for n in notes:
                print(f"  integrity[{tag}] {n}")
            integ_ok = integ_ok and good
        if not integ_ok:
            print("  >>> INTEGRITY REFUSED - not quotable")
            rc = 1
            continue

        hd, ha = final_round_hits(root / arms["default"]), final_round_hits(root / arms["ablated"])
        md, ma = manifest(root / arms["default"]), manifest(root / arms["ablated"])
        print(f"  default run_dir : {arms['default']}")
        print(f"  ablated run_dir : {arms['ablated']}")
        print(f"  top1 threshold  : {md['threshold']} / {ma['threshold']}")
        print(f"  components off  : {md['disabled']} / {ma['disabled']}")
        print(f"  seed_block      : {md['seed_block']}")

        q0, q1 = ridge_rows(root / arms["default"]), ridge_rows(root / arms["ablated"])
        checks = {
            "ridge_bit_identical": bool(q0 and q0 == q1),
            "seed_block_match": md["seed_block"] == ma["seed_block"],
            "threshold_match": md["threshold"] == ma["threshold"],
            "wt_anchor_state_flipped": md["wt_anchor_state"] is True
                                       and ma["wt_anchor_state"] is False,
        }
        print(f"  ridge rows compared: {len(q0)} vs {len(q1)} (all columns bar wall-clock)")
        print("  PAIRING: " + " ".join(f"{k}={v}" for k, v in checks.items()))
        if not checks["ridge_bit_identical"]:
            diff = [k for k in sorted(set(q0) | set(q1)) if q0.get(k) != q1.get(k)]
            print(f"    ridge MISMATCH on {len(diff)} (seed,round) cell(s); "
                  f"first 5: {diff[:5]}")
        if not all(checks.values()):
            print("  >>> PAIRING REFUSED - not reportable as a paired result")
            rc = 1
            continue
        print("  >>> PAIRING VERIFIED")

        pd_, pa_, rg = hd["panoptes_default"], ha["panoptes_default"], hd["ridge"]
        for nm, per in (("ridge", rg), ("panoptes default", pd_), ("panoptes ABLATED", pa_)):
            arr = np.array([per[s] for s in sorted(per)], dtype=float)
            print(f"  mean {nm:<18} = {arr.mean():8.3f}  (sd {arr.std(ddof=1):6.3f}, "
                  f"n {len(arr)})")
        print("  -- paired contrasts (rule 2.11) --")
        print(line("PRIMARY  ablated - default", paired(pa_, pd_)))
        print(line("         default - ridge", paired(pd_, rg)))
        print(line("         ablated - ridge", paired(pa_, rg)))
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

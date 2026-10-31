"""
alde_compat_probe.py - can the published ALDE implementation run on the new corpus?

Probes the six newly committed cross-protein landscapes (plus trpb/gb1/her2 as references)
against the AUTHORS' OWN code at https://github.com/jsunn-y/ALDE, commit d0b3593dd179,
run in ALDE's own interpreter. Nothing in the ALDE checkout is modified or reimplemented.

Two passes:

  encode  - build the per-site combo strings from `mutant` + the parent sequence, hand them to
            src.encoding_utils.generate_onehot unmodified, and report the encoding geometry
            plus the state of ALDE's own design-space normalisation
            (botorch.utils.transforms.normalize over utils.domain_discrete bounds), which is
            the first place a non-saturated library can break.

  run     - run src.optimize.BayesianOptimization.run end to end for N seeds with exactly the
            configuration benchmark_panoptes_isolated.py uses for its `alde` arm:
            DNN_ENSEMBLE, architecture [d, 30, 30, 1], TS, lrelu, lr 1e-3, train_iter 300,
            dropout 0, 25 shared initial observations, budget 100, batch_size 25.
            Reports whether it completed, and whether the trajectory it produced is a real
            search or a degenerate one.

MUST be run with ALDE's interpreter, not the repo venv:
  C:\\Users\\stoof\\OneDrive\\Documents\\GitHub\\.venv-alde\\Scripts\\python.exe

Usage:
  python alde_compat_probe.py encode --out <json>
  python alde_compat_probe.py run --landscape <key> --seeds 3 --out <json>
"""
from __future__ import annotations

import argparse
import json
import os
import random
import re
import sys
import time
import traceback
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import torch

ALDE_CHECKOUT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\ALDE")
REPO = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio")
DATA = REPO / "scripts" / "data"
WT_JSON = DATA / "library_wt_sequences.json"

sys.path.insert(0, str(ALDE_CHECKOUT))

# The six new landscapes, plus the three already-registered libraries as controls.
LANDSCAPES = {
    # key: (csv name, labelled stage filter or None)
    "creilov":    ("CREILOV_CHLRE_Chen_2023.csv", None),
    "cr9114_h1":  ("CR9114_H1_Phillips_2021.csv", None),
    "eqfp611_blue": ("EQFP611_BLUE_Poelwijk_2019.csv", None),
    "eqfp611_red":  ("EQFP611_RED_Poelwijk_2019.csv", None),
    "cr6261_h1":  ("CR6261_H1_Phillips_2021.csv", None),
    "cr6261_h9":  ("CR6261_H9_Phillips_2021.csv", None),
    # controls: the two landscapes ALDE was reported strongest on, and HER2.
    "trpb":       ("TRPB_THEMA_Johnston_2024.csv", None),
    "gb1":        ("GB1_STRSG_Wu_2016.csv", None),
    "her2":       ("HER2_TRASTUZUMAB_Mason_2021.csv", "train"),
}

# Parent residues for the three control libraries, which are not in library_wt_sequences.json.
CONTROL_PARENTS = {
    "trpb": ([183, 184, 227, 228], "VFVS"),
    "gb1":  ([39, 40, 41, 54], "VDGV"),
    "her2": ([99, 100, 101, 102, 103, 104, 105, 106, 107, 108], "WGGDGFYAMD"),
}

MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def load_landscape(key: str):
    """Returns (combos, fitness, positions, parent_residues) for one landscape.

    The genotype is one residue per LIBRARY POSITION, which is what ALDE's own
    generate_onehot consumes and what benchmark_panoptes_isolated.py hands its alde arm
    (`view.combos`). The full-length parent is used only to fill the silent positions.
    """
    csv_name, stage = LANDSCAPES[key]
    frame = pd.read_csv(DATA / csv_name, usecols=lambda c: c in ("mutant", "DMS_score", "stage"))
    if stage is not None:
        frame = frame[frame["stage"].astype(str) == stage]
    frame = frame[frame["mutant"].notna()]

    records = []
    positions_seen = {}
    for mutant, score in zip(frame["mutant"].astype(str), frame["DMS_score"].astype(float)):
        subs = []
        ok = True
        for token in mutant.split(":"):
            token = token.strip()
            if not token:
                continue
            match = MUT.match(token)
            if not match:
                ok = False
                break
            wt, pos, mt = match.group(1), int(match.group(2)), match.group(3)
            subs.append((wt, pos, mt))
            positions_seen.setdefault(pos, wt)
        if ok:
            records.append((subs, score))

    if key in CONTROL_PARENTS:
        positions, parent = CONTROL_PARENTS[key]
    else:
        wt_table = json.loads(WT_JSON.read_text(encoding="utf-8"))
        entry = wt_table[csv_name]
        sequence = entry["wild_type_sequence"]
        offset = int(entry["position_offset"])
        positions = sorted(positions_seen)
        parent = "".join(sequence[p - offset] for p in positions)
        # cross-check the registry residue against the parent sequence
        for i, p in enumerate(positions):
            if positions_seen[p] != parent[i]:
                raise SystemExit(
                    f"{key}: parent mismatch at position {p}: mutant strings say "
                    f"{positions_seen[p]}, parent sequence says {parent[i]}"
                )

    index_of = {p: i for i, p in enumerate(positions)}
    combos, fitness, dropped = [], [], 0
    for subs, score in records:
        residues = list(parent)
        bad = False
        for wt, pos, mt in subs:
            slot = index_of.get(pos)
            if slot is None or residues[slot] != wt:
                bad = True
                break
            residues[slot] = mt
        combo = "".join(residues)
        if bad or combo == parent:  # the parent itself is not a library member
            dropped += 1
            continue
        combos.append(combo)
        fitness.append(score)

    return combos, np.asarray(fitness, dtype=np.float64), positions, parent, dropped


def encode_report(key: str):
    from src.encoding_utils import generate_onehot, ALL_AAS
    import botorch
    import src.utils as alde_utils

    out = {"landscape": key, "csv": LANDSCAPES[key][0]}
    try:
        combos, fitness, positions, parent, dropped = load_landscape(key)
    except Exception as exc:
        out["status"] = "FAILS"
        out["stage"] = "genotype reconstruction"
        out["error"] = f"{type(exc).__name__}: {exc}"
        return out

    out["variants"] = len(combos)
    out["dropped_rows"] = dropped
    out["sites"] = len(positions)
    out["parent"] = parent
    out["parent_length_full_protein"] = None
    if key not in CONTROL_PARENTS:
        wt_table = json.loads(WT_JSON.read_text(encoding="utf-8"))
        out["parent_length_full_protein"] = wt_table[LANDSCAPES[key][0]]["length"]

    alphabets = [sorted({c[i] for c in combos} | {parent[i]}) for i in range(len(positions))]
    out["alphabet_sizes"] = [len(a) for a in alphabets]
    out["max_alphabet"] = max(out["alphabet_sizes"])
    out["is_binary_library"] = set(out["alphabet_sizes"]) == {2}
    out["characters_outside_ALL_AAS"] = sorted(
        {c for combo in combos for c in combo} - set(ALL_AAS)
    )

    # --- ALDE's own one-hot, unmodified -------------------------------------------------
    try:
        started = time.time()
        encoded = generate_onehot(combos)
        encoded = torch.reshape(encoded, (encoded.shape[0], -1))
        out["onehot_seconds"] = round(time.time() - started, 2)
        out["generate_onehot"] = "OK"
        out["encoding_dim"] = int(encoded.shape[1])
        out["encoding_dim_formula"] = f"{len(positions)} sites x {len(ALL_AAS)} (ALL_AAS) = {encoded.shape[1]}"
        # every row must be a valid one-hot over each site
        per_site = encoded.reshape(encoded.shape[0], len(positions), len(ALL_AAS))
        out["rows_are_valid_onehot"] = bool(torch.all(per_site.sum(dim=2) == 1).item())
    except Exception as exc:
        out["status"] = "FAILS"
        out["stage"] = "src.encoding_utils.generate_onehot"
        out["error"] = f"{type(exc).__name__}: {exc}"
        out["traceback"] = traceback.format_exc()
        return out

    # --- ALDE's design-space normalisation ----------------------------------------------
    # BayesianOptimization.__init__ does:
    #   self.disc_X = botorch.utils.transforms.normalize(disc_X, domain)
    # with domain = utils.domain_discrete(X) = (per-column min, per-column max).
    # A one-hot column that is CONSTANT across the design space gives upper == lower, so the
    # normalisation divides by zero.
    disc_X = encoded.double()
    lower, upper = alde_utils.domain_discrete(disc_X)
    constant = (upper == lower)
    out["constant_onehot_columns"] = int(constant.sum().item())
    out["constant_column_fraction"] = round(float(constant.float().mean().item()), 4)
    domain = (lower, upper)
    normalised = botorch.utils.transforms.normalize(disc_X, torch.stack(domain))
    out["nan_columns_after_normalize"] = int(torch.isnan(normalised).any(dim=0).sum().item())
    out["nan_cells_after_normalize"] = int(torch.isnan(normalised).sum().item())
    out["design_space_is_all_nan_rows"] = bool(
        torch.isnan(normalised).any(dim=1).all().item()
    ) if out["nan_cells_after_normalize"] else False

    out["status"] = "FAILS" if out["nan_columns_after_normalize"] else "RUNS"
    return out


def run_alde(key: str, seeds: int, initial: int, rounds: int, batch: int,
             train_iter: int, work_root: Path, guard_constant: bool = False):
    """Runs the authors' loop end to end, exactly as the benchmark's alde arm configures it.

    `guard_constant` is a DIAGNOSTIC DEVIATION from the published code, not a configuration
    option and not something the benchmark does. It widens the reported domain of every
    one-hot column that is constant across the design space, from (0, 0) to (0, 1), so that
    BayesianOptimization.__init__'s `normalize(disc_X, domain)` divides by 1 rather than by 0
    and leaves the column at its constant value instead of NaN. It changes no other number:
    on a library with no constant columns it is a no-op. It exists only to establish whether
    the 0/0 in the normalisation is the SOLE cause of the degenerate trajectory, and its
    results must be reported as a deviation from the authors' implementation.
    """
    from src.encoding_utils import generate_onehot
    from src.optimize import BayesianOptimization, BO_ARGS
    import src.objectives as objectives

    combos, fitness, positions, parent, _ = load_landscape(key)
    protein = f"{key}_probe"

    work_root.mkdir(parents=True, exist_ok=True)
    os.chdir(work_root)
    data_dir = Path("data") / protein
    data_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"Combo": combos, "fitness": fitness}).to_csv(
        data_dir / "fitness.csv", index=False)
    if not (data_dir / "onehot_x.pt").exists():
        encoded = generate_onehot(combos)
        encoded = torch.reshape(encoded, (encoded.shape[0], -1))
        torch.save(encoded, data_dir / "onehot_x.pt")

    obj = objectives.Combo(protein, "onehot")
    domain = obj.get_domain()
    disc_X = obj.get_points()[0]
    disc_y = obj.get_points()[1]

    n_constant = int((domain[1] == domain[0]).sum().item())
    if guard_constant and n_constant:
        lower, upper = domain[0].clone(), domain[1].clone()
        flat = upper == lower
        upper[flat] = lower[flat] + 1.0
        domain = (lower, upper)

    # top-1% hit set, for an honest read of the trajectory (evaluator-side only)
    threshold = float(np.percentile(fitness, 99))
    hit_set = set(np.flatnonzero(fitness >= threshold).tolist())

    results = []
    for seed in range(seeds):
        torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
        rng = np.random.default_rng(seed)
        start_indices = torch.tensor(
            np.sort(rng.choice(len(combos), size=initial, replace=False)), dtype=torch.long)
        start_x = disc_X[start_indices]
        start_y = disc_y[start_indices]
        arc = [domain[0].size(-1), 30, 30, 1]
        tag = f"{key}_guarded_seed{seed}" if guard_constant else f"{key}_seed{seed}"
        savedir = str(work_root / "runs" / tag / "run_")
        Path(savedir).parent.mkdir(parents=True, exist_ok=True)

        args = BO_ARGS(
            mtype="DNN_ENSEMBLE", kernel="RBF", acq_fn="TS", xi=4,
            architecture=arc, activation="lrelu", min_noise=1e-6, trainlr=1e-3,
            train_iter=train_iter, dropout=0, mcdropout=0, verbose=0,
            bb_fn=obj.objective, domain=domain, disc_X=disc_X, disc_y=disc_y,
            noise_std=0, n_rand_init=0, budget=rounds * batch, query_cost=1,
            queries_x=start_x, queries_y=start_y, indices=start_indices,
            savedir=savedir, batch_size=batch,
        )
        record = {"seed": seed}
        started = time.time()
        try:
            BayesianOptimization.run(args, seed)
            selected = torch.load(savedir + "indices.pt").reshape(-1).tolist()
            queried = [int(i) for i in selected[initial:]]
            record["status"] = "completed"
            record["n_queried"] = len(queried)
            record["n_unique_queried"] = len(set(queried))
            record["n_distinct_in_whole_trajectory"] = len(set(int(i) for i in selected))
            record["duplicate_picks"] = len(queried) - len(set(queried))
            record["first_ten_picks"] = queried[:10]
            record["hits_in_start"] = len(set(start_indices.tolist()) & hit_set)
            record["hits_total"] = len(set(int(i) for i in selected) & hit_set)
            record["best_fitness_found"] = float(
                max(fitness[int(i)] for i in selected))
        except Exception as exc:
            record["status"] = "failed"
            record["error"] = f"{type(exc).__name__}: {exc}"
            record["traceback"] = traceback.format_exc()[-2500:]
        record["seconds"] = round(time.time() - started, 1)
        results.append(record)
        print(f"  seed {seed}: {record['status']} "
              f"unique={record.get('n_unique_queried')} hits={record.get('hits_total')} "
              f"best={record.get('best_fitness_found')} ({record['seconds']}s)", flush=True)

    return {
        "landscape": key,
        "variants": len(combos),
        "sites": len(positions),
        "top1_threshold": threshold,
        "top1_hit_set_size": len(hit_set),
        "constant_onehot_columns": n_constant,
        "published_code_unmodified": not (guard_constant and n_constant),
        "deviation": ("domain of constant one-hot columns widened from (0,0) to (0,1) so that "
                      "normalize() does not divide by zero") if (guard_constant and n_constant) else None,
        "config": {"initial": initial, "rounds": rounds, "batch_size": batch,
                   "budget": rounds * batch, "train_iter": train_iter,
                   "model": "DNN_ENSEMBLE [d,30,30,1]", "acq": "TS"},
        "seeds": results,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("encode", "run"))
    parser.add_argument("--landscape", default=None)
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--initial", type=int, default=25)
    parser.add_argument("--rounds", type=int, default=4)
    parser.add_argument("--batch", type=int, default=25)
    parser.add_argument("--train-iter", type=int, default=300)
    parser.add_argument("--work", default=None)
    parser.add_argument("--guard-constant-columns", action="store_true",
                        help="DIAGNOSTIC DEVIATION from the published code: widen the domain of "
                             "constant one-hot columns so normalize() does not divide by zero. "
                             "Results obtained with this flag are not ALDE as published.")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if args.mode == "encode":
        keys = [args.landscape] if args.landscape else list(LANDSCAPES)
        reports = []
        for key in keys:
            print(f"== {key}", flush=True)
            report = encode_report(key)
            reports.append(report)
            print(json.dumps(report, indent=2), flush=True)
        out_path.write_text(json.dumps(reports, indent=2), encoding="utf-8")
    else:
        if not args.landscape:
            raise SystemExit("--landscape is required for mode=run")
        work = Path(args.work) if args.work else Path.cwd() / "alde_probe_work"
        report = run_alde(args.landscape, args.seeds, args.initial, args.rounds,
                          args.batch, args.train_iter, work,
                          guard_constant=args.guard_constant_columns)
        out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps({k: v for k, v in report.items() if k != "seeds"}, indent=2))
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()

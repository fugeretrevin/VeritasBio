#!/usr/bin/env bash
# Generality sweep: does library-aware candidate generation help where reachability was
# the bottleneck, and leave the saturated libraries alone?
#
# Pre-registered prediction (docs/day1-reachability.md section 7a):
#   improvement on the reachability-limited landscapes, and NO CHANGE on GB1, which is
#   100% reachable at depth <= 4 and so has nothing for these fixes to add.
#
# HARD GUARD: refuses to start if a benchmark is already running. An earlier attempt had
# two sweeps execute concurrently for three hours after a task-stop ended the task but
# left its shell alive; that voided ten run directories through CPU contention, a mixed
# code state and a mid-run source change (section 11). Waiting was the wrong behaviour -
# it let a second sweep queue up behind a stale one nobody knew was alive.

set -u
# Resolve the repository root from this script's own location rather than a
# hardcoded path, so the sweep runs from any checkout.
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)" || exit 1
PY=".venv/Scripts/python.exe"
COMMON="--seeds 48 --base-seed 0 --rounds 4 --provenance-mode acknowledge-unflagged"
METHODS="random,ridge,panoptes_default"

running=$(powershell.exe -NoProfile -NonInteractive -Command \
  "@(Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | Where-Object { \$_.CommandLine -match 'benchmark_panoptes' }).Count" \
  2>/dev/null | tr -cd '0-9')
[ -z "$running" ] && running=0
if [ "$running" -ne 0 ]; then
  echo "REFUSING TO START: $running benchmark process(es) already running."
  echo "Tree-kill them first:  taskkill //F //T //PID <shell pid>"
  exit 1
fi
echo "no benchmark process running; starting sweep at $(date +%H:%M:%S) local"

# corpus.derived_table flags seven landscapes reading a CombinGym "clean" workbook while a
# raw table exists at the pin. Reviewed by measuring observed-vs-expected genotype coverage
# PER DEPTH on all nine landscapes (scripts/analysis/depth_coverage.py): no depth is at zero
# anywhere, eqFP611 and TrpB are complete, CR6261-H9's deep half is slightly higher than
# its shallow half, and the largest gradient is CreiLOV at 4.5 points - none shows the
# CR9114 pattern of attrition reaching 100% at maximum depth. See section 10.
WHY="reviewed 2026-10-31 by measuring observed-vs-expected coverage per depth; no depth at zero, deep-half coverage within 4.5 points of shallow, unlike the CR9114 clean workbook which fell to zero at depth 16"
# An ARRAY, not a string: the reason contains spaces and an unquoted string word-splits.
ACK=()
for k in cr6261_h1 cr6261_h9 creilov eqfp611_blue eqfp611_combined eqfp611_red gb1_combingym; do
  ACK+=(--acknowledge-derived "${k}=${WHY}")
done

run_one () {
  local key="$1"; local initial="$2"; local batch="$3"; local why="$4"
  local budget=$((initial + 4 * batch))
  echo ""
  echo "################################################################"
  echo "## $key   initial=$initial batch=$batch budget=$budget   start $(date +%H:%M:%S)"
  echo "## $why"
  echo "################################################################"

  # Rule 2.7: the gate runs before any scored run, at the budget actually about to run.
  "$PY" -X utf8 -m veritas.preflight --library "$key" \
      --budget "$budget" --initial-size "$initial" "${ACK[@]}" 2>&1 | tail -4
  local gate=${PIPESTATUS[0]}
  if [ "$gate" -ne 0 ]; then
    echo "!! PREFLIGHT FAILED for $key (exit $gate) - DO NOT QUOTE THIS LANDSCAPE"
  else
    echo "   preflight OK for $key"
  fi

  "$PY" -X utf8 benchmark_panoptes_isolated.py \
      --library "$key" --methods "$METHODS" \
      --initial-size "$initial" --batch-size "$batch" \
      $COMMON 2>&1 | grep -E "Run state|Artifacts|^ +(random|ridge|panoptes)|DEGENERATE RUN|^[0-9]+ round"
  echo "## $key done $(date +%H:%M:%S) exit=${PIPESTATUS[0]}"
}

# Scope trimmed to what earns its runtime. Dropped: TrpB (already measured bit-identical
# on code carrying every generator fix, and it declares no per-site alphabet so the
# encoding change cannot touch it); CR6261 H1/H9 (19-variant top-1% tier, low power by
# construction, and one protein against two antigens so not independent of each other).
#
# Fast landscapes first, so signal arrives early, with CreiLOV third rather than last -
# it is both the most informative and the slowest, and should not be the thing lost if
# the sweep is interrupted.

# 1-2. eqFP611 red and blue: 13 sites binary, 100% genotype coverage, same geometry with
#    the parent axis varied. Red's prior 48-seed run found 0 top-1% hits on every seed.
run_one eqfp611_red  17 16 "13 sites binary, NON-functional parent, prior run found 0 hits"
run_one eqfp611_blue 17 16 "same geometry as red, functional parent - paired control"

# 3. CreiLOV: 15 sites, the ONLY non-uniform alphabet (2 x12, 3 x2, 5 x1), so it tests the
#    per-site path rather than a uniform binary one. Carries the 4.5-point coverage caveat.
run_one creilov 25 25 "15 sites, mixed alphabet - highest-value generality test"

# 4. GB1: 4 sites, saturated, 100% reachable at depth <= 4. NEGATIVE CONTROL. If this
#    moves, these are not geometry fixes and the diagnosis is wrong.
run_one gb1 25 25 "NEGATIVE CONTROL - saturated 4-site, must not change"

# 5. CR9114: re-run so `ridge` gets the corrected 32-column encoding instead of 320
#    columns of which 288 were identically zero. DEVELOPMENT DATA (rule 2.4) - the
#    intervention was designed against its measured failure, so this is mechanism
#    evidence, not generalization evidence.
run_one cr9114_h1 25 25 "16 sites binary - DEVELOPMENT DATA, fair-encoding re-run"

# 6. HER2: already pool_restricted, so candidate_pool is supplied and
#    generate_candidate_variants is never called. NULL CONTROL: nothing here can reach it.
run_one her2 25 22 "NULL CONTROL - pool_restricted, never invokes the generator"

echo ""
echo "################ sweep complete $(date +%H:%M:%S) ################"

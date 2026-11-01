<#
Baseline panel: how strong is `ridge` relative to its neighbours, and how many of the
engine's results survive against the strongest baseline rather than against ridge alone?

Every comparative claim in this repository is against exactly ONE baseline - a one-hot
ridge regression - so "beats ridge" is being read as "beats the baselines". Two facts make
that unsafe: on CR9114 `ridge` was measured as NOT SEPARABLE FROM RANDOM, so the project's
one win is a win over a comparator that may be no better than chance; and on TrpB `ridge`
found the global optimum the engine never did, so ridge is not uniformly weak either.

Pre-registration, written and committed before any arm here ran:
docs/prereg-baseline-panel.md.

THIS IS NOT A SECOND HARNESS. The harness is benchmark_panoptes_isolated.py. This is
generality_sweep.ps1's invocation - the same hard guard, the same derived-table
acknowledgements, the same preflight-before-scored-run ordering, the same Invoke-Captured
stream plumbing - with the panel's arms and the panel's three landscapes. It is a sibling
rather than a flag on generality_sweep.ps1 because that script is the record of what the
2026-11-01 cross-protein sweep actually ran, and editing it would rewrite that record.

NOTHING HERE PROMOTES A BASELINE. `ridge` remains the registered reference baseline;
docs/HANDOFF.md 5.4 makes promoting another model to reference status an amendment to a
frozen hypothesis, which needs an explicit dated amendment record and not a script.
DEFAULT_METHODS is unchanged and every added arm is opt-in.

HARD GUARD, inherited verbatim in intent: refuses to start if a benchmark is already
running. An earlier night had two sweeps execute concurrently for three hours after a
task-stop ended the task but left its shell alive; that voided ten run directories through
CPU contention, a mixed code state and a mid-run source change
(docs/day1-reachability.md section 11). Waiting is the wrong behaviour - it lets a second
sweep queue up behind a stale one nobody knew was alive.

Rule 3 of that section: read progress from the run you started, BY MANIFEST, never from
"the newest directory". Each landscape's own run directory is recorded into the sweep log
as it goes.
#>
[CmdletBinding()]
param(
    # Landscape keys to run, in order. Default is the full pre-registered panel.
    [string[]] $Only,
    # Calibration: override the seed count (the panel is registered at 48).
    [int] $Seeds = 48,
    [switch] $WhatIfOnly
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

# --- interpreter and results root, resolved for a WORKTREE as well as the real checkout ---
# generality_sweep.ps1 assumes $repo/.venv exists. In a git worktree it does not: .venv is
# gitignored and lives only in the real repository, which a worktree CANNOT SEE. That has
# cost three full agent runs (docs/HANDOFF.md 2.8), so resolution is explicit here and every
# failure NAMES EVERY PATH IT CHECKED.
#
# The source tree stays the worktree's - that is the code under test, and the harness hashes
# it - while the interpreter and the results root come from the real repository root.
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
Set-Location $repo

function Resolve-RealRepositoryRoot {
    param([string] $Checkout)
    $checked = @()
    $dotGit = Join-Path $Checkout '.git'
    $checked += $dotGit
    if (Test-Path $dotGit -PathType Container) { return @{ Root = $Checkout; Checked = $checked } }
    if (Test-Path $dotGit -PathType Leaf) {
        # Worktree: `.git` is a FILE reading "gitdir: <real>/.git/worktrees/<name>".
        $line = (Get-Content $dotGit -TotalCount 1)
        if ($line -match '^gitdir:\s*(.+)$') {
            $gitDir = $Matches[1].Trim()
            $checked += $gitDir
            # <real>/.git/worktrees/<name>  ->  <real>
            $candidate = Split-Path (Split-Path (Split-Path $gitDir -Parent) -Parent) -Parent
            $checked += $candidate
            if ($candidate -and (Test-Path (Join-Path $candidate '.git'))) {
                return @{ Root = $candidate; Checked = $checked }
            }
        }
    }
    return @{ Root = $null; Checked = $checked }
}

$resolution = Resolve-RealRepositoryRoot -Checkout $repo
$realRoot = $resolution.Root
if (-not $realRoot) {
    throw ("could not resolve the real repository root from the checkout at $repo. " +
           "Checked: " + ($resolution.Checked -join '; '))
}

$pyCandidates = @(
    (Join-Path $repo '.venv\Scripts\python.exe'),
    (Join-Path $realRoot '.venv\Scripts\python.exe')
)
$py = $pyCandidates | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $py) { throw ("no interpreter found. Checked: " + ($pyCandidates -join '; ')) }

# Artifacts go to the REAL repository's results root, which is where every script in
# scripts/analysis reads them from - a hard-coded absolute path in several of them. Writing
# them into the worktree would leave the read-out tools pointing at an empty directory while
# the run looked successful. The path is passed EXPLICITLY, as a constraint on the run.
$logDir = Join-Path $realRoot 'benchmark_comparison_results'
if (-not (Test-Path $logDir)) { throw ("results root absent. Checked: $logDir") }

# The pinned CombinGym corpus the preflight gate's corpus.lock check verifies against. It
# lives in the gitignored downloads/ at the REAL root, so from a worktree the bare gate
# reports "FAILED: corpus.lock - 29 of 29 file(s) absent" and all three landscapes look
# un-quotable when in fact the corpus verifies 29/29. That is rule 2.8 in its purest form:
# a bare "absent" is what let our strongest baseline get silently dropped from a 48-seed
# run, so the path is passed explicitly and every candidate is named on failure.
$corpusCandidates = @(
    (Join-Path $repo 'downloads\combingym'),
    (Join-Path $realRoot 'downloads\combingym')
)
$corpusDir = $corpusCandidates | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $corpusDir) { throw ("pinned CombinGym corpus absent. Checked: " + ($corpusCandidates -join '; ')) }

Write-Host "source tree : $repo"
Write-Host "interpreter : $py"
Write-Host "results root: $logDir"
Write-Host "corpus      : $corpusDir"

# --- the panel ---------------------------------------------------------------------------
# BASELINE_PANEL_METHODS in benchmark_panoptes_isolated.py, in panel order. `random` is NOT
# optional: a panel without it cannot calibrate the others, because nothing then says whether
# a landscape's top-1% count can discriminate methods at all. That is the HER2 informativeness
# lesson (docs/HANDOFF.md section 1) applied before the fact rather than after it.
$methods = 'random,ridge,random_forest,gbm_onehot,knn_hamming,gp_onehot,greedy_measured,panoptes_default'

# --- hard guard ---------------------------------------------------------------------------
$running = @(Get-CimInstance Win32_Process -Filter "Name='python.exe'" |
             Where-Object { $_.CommandLine -and $_.CommandLine -match 'benchmark_panoptes' })
if ($running.Count -ne 0) {
    Write-Host "REFUSING TO START: $($running.Count) benchmark process(es) already running." -ForegroundColor Red
    $running | ForEach-Object { Write-Host "  PID $($_.ProcessId): $($_.CommandLine)" }
    Write-Host "Tree-kill them first:  taskkill /F /T /PID <pid>" -ForegroundColor Yellow
    exit 1
}
Write-Host "no benchmark process running; starting panel at $(Get-Date -Format HH:mm:ss) local"

# --- derived-table acknowledgements -------------------------------------------------------
# Carried over unchanged from generality_sweep.ps1 so the provenance posture of the panel is
# identical to the sweep it is compared against. corpus.derived_table flags seven landscapes
# reading a CombinGym "clean" workbook while a raw table exists at the pin; reviewed by
# measuring observed-vs-expected genotype coverage PER DEPTH on all nine landscapes
# (scripts/analysis/depth_coverage.py).
$why = 'reviewed 2026-10-31 by measuring observed-vs-expected coverage per depth; no depth at zero, deep-half coverage within 4.5 points of shallow, unlike the CR9114 clean workbook which fell to zero at depth 16'
$ack = @()
foreach ($k in @('cr6261_h1','cr6261_h9','creilov','eqfp611_blue','eqfp611_combined','eqfp611_red','gb1_combingym')) {
    $ack += '--acknowledge-derived'
    $ack += "$k=$why"
}

# --- landscapes, in the pre-registered priority order -------------------------------------
# trpb FIRST and on seeds 48-95. Those seeds were never looked at during the 2^5
# configuration scan (scripts/analysis/heldout.py), and TrpB is the one landscape where the
# corrected qualification bar accepts the identical variant set, so nothing in the
# 2026-10-31 fix set changed what is being searched there.
#
# cr9114_h1 SECOND: it is the project's only win, and it is the landscape where the
# ridge-versus-random question is live. DEVELOPMENT DATA under rule 2.4 - the intervention
# was designed against its measured failure - so whatever it shows is mechanism evidence,
# never generalization evidence, and the read-out must say so.
#
# creilov THIRD: the worst loss (-50.29 on 48 of 48 seeds). It is the test of whether
# `ridge` is strong where the engine lost, which decides whether that loss is real or an
# artefact of an unfairly strong comparator.
$plan = @(
    [pscustomobject]@{ Key='trpb';      Initial=25; Batch=25; BaseSeed=48; Why='HELD-OUT seeds 48-95; the one landscape untouched by the 2026-10-31 fix set' }
    [pscustomobject]@{ Key='cr9114_h1'; Initial=25; Batch=25; BaseSeed=0;  Why='the one win, and where ridge was NOT separable from random - DEVELOPMENT DATA' }
    [pscustomobject]@{ Key='creilov';   Initial=25; Batch=25; BaseSeed=0;  Why='the worst loss, -50.29 on 48/48 seeds - is ridge strong here?' }
)
if ($Only) { $plan = $plan | Where-Object { $Only -contains $_.Key } }
if (-not $plan) { throw "no landscapes selected; -Only must name keys from the plan" }

# Every file this script writes carries the same start stamp. Without it the per-landscape
# gate and stdout files have fixed names, and a second invocation - or an abandoned first
# one whose handles have not closed - collides with them: the panel aborted on
# "Access to the path panel_gate_trpb.txt is denied" after an interrupted attempt left the
# file open, which looks like a permissions problem and is not one. Several sessions share
# this results root, so a fixed name is a shared mutable name.
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$sweepLog = Join-Path $logDir ("panel_{0}.log" -f $stamp)
"baseline panel started $(Get-Date -Format 'o')  seeds=$Seeds  methods=$methods" | Out-File $sweepLog -Encoding utf8
"source tree $repo" | Add-Content $sweepLog
Write-Host "panel log: $sweepLog"

# Landscapes whose gate failed. Collected rather than aborting: a failed gate means that
# landscape must not be quoted (rule 2.7), not that the rest of the panel is worthless.
$failed = @()
#: Landscapes whose POST-RUN gate failed, reported separately. The pre-run gate and the
#: post-run gate answer different questions - the pre-run one cannot see the run at all -
#: so collapsing them into one list would hide which checks actually ran.
$postFailed = @()

function Invoke-Captured {
    <#
      Run a native command, returning @{ Lines; ExitCode } with stdout and stderr merged.

      Deliberately NOT `& $py ... 2>&1` under 'Stop'. Windows PowerShell 5.1 wraps each
      stderr line from a native executable in a NativeCommandError ErrorRecord, which under
      $ErrorActionPreference='Stop' throws on the first line the tool logs. Both the
      benchmark and the preflight gate log at INFO to stderr, so a whole sweep once aborted
      seven seconds in, after a gate that had passed.

      The call operator with splatting, NOT Start-Process -ArgumentList: Start-Process joins
      the argument array into one string and lets the child re-parse it, which word-splits
      any argument containing spaces. The --acknowledge-derived KEY=WHY reasons are full
      sentences, so every one of the seven arrived as a dozen stray positional arguments and
      the gate exited 2 with "unrecognized arguments".
    #>
    param([string] $Exe, [string[]] $Arguments, [string] $Tag)

    $prev = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try {
        $raw = & $Exe @Arguments 2>&1
        $code = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $prev
    }
    $lines = @($raw | ForEach-Object { [string] $_ })
    return @{ Lines = $lines; ExitCode = $code }
}

foreach ($p in $plan) {
    $budget = $p.Initial + 4 * $p.Batch
    $banner = ("## $($p.Key)   initial=$($p.Initial) batch=$($p.Batch) budget=$budget " +
               "base-seed=$($p.BaseSeed) seeds=$Seeds   start $(Get-Date -Format HH:mm:ss)")
    Write-Host ''
    Write-Host '################################################################'
    Write-Host $banner
    Write-Host "## $($p.Why)"
    Write-Host '################################################################'
    Add-Content $sweepLog "`n$banner`n## $($p.Why)"

    if ($WhatIfOnly) {
        Write-Host "   (WhatIfOnly) would run preflight + benchmark for $($p.Key)"
        continue
    }

    # Rule 2.7: the gate runs before any scored run, at the budget actually about to run.
    $gateArgs = @('-X','utf8','-m','veritas.preflight','--library',$p.Key,
                  '--budget',$budget,'--initial-size',$p.Initial,
                  '--corpus-dir',$corpusDir) + $ack
    $gateRes = Invoke-Captured -Exe $py -Arguments $gateArgs -Tag "gate_$($p.Key)"
    $gate = $gateRes.ExitCode
    $gateRes.Lines | Out-File (Join-Path $logDir "panel_gate_$($p.Key)_$stamp.txt") -Encoding utf8
    ($gateRes.Lines | Select-Object -Last 8) | ForEach-Object { Write-Host "   $_"; Add-Content $sweepLog "   $_" }
    if ($gate -ne 0) {
        $msg = "!! PREFLIGHT FAILED for $($p.Key) (exit $gate) - DO NOT QUOTE THIS LANDSCAPE"
        Write-Host $msg -ForegroundColor Red
        Add-Content $sweepLog $msg
        $failed += $p.Key
    } else {
        Write-Host "   preflight OK for $($p.Key)" -ForegroundColor Green
        Add-Content $sweepLog "   preflight OK for $($p.Key)"
    }

    $runArgs = @('-X','utf8','benchmark_panoptes_isolated.py',
                 '--library',$p.Key,'--methods',$methods,
                 '--initial-size',$p.Initial,'--batch-size',$p.Batch,
                 '--seeds',$Seeds,'--base-seed',$p.BaseSeed,'--rounds','4',
                 '--provenance-mode','acknowledge-unflagged',
                 '--output-dir',$logDir)
    Add-Content $sweepLog ("   argv: " + ($runArgs -join ' '))
    $runRes = Invoke-Captured -Exe $py -Arguments $runArgs -Tag "run_$($p.Key)"
    $runOut = $runRes.Lines
    $code = $runRes.ExitCode

    # Record the run directory this landscape actually produced, so progress is never read
    # from "the newest directory" (rule 3) - during contention the newest was another
    # sweep's, and TrpB was reported as CreiLOV for an hour.
    $runOut | Out-File (Join-Path $logDir "panel_$($p.Key)_stdout_$stamp.txt") -Encoding utf8
    $runOut | Select-String -Pattern 'Run state|Artifacts|DEGENERATE RUN|^\s*(random|ridge|random_forest|gbm_onehot|knn_hamming|gp_onehot|greedy_measured|panoptes)|^\d+ round' |
        ForEach-Object { Write-Host "   $_"; Add-Content $sweepLog "   $_" }
    $hit = $runOut | Select-String -Pattern 'run_\d{8}_\d{6}_[0-9a-f]+' | Select-Object -First 1
    $dir = if ($hit) { $hit.Matches[0].Value } else { '(none reported)' }
    $done = "## $($p.Key) done $(Get-Date -Format HH:mm:ss) exit=$code run_dir=$dir"
    Write-Host $done
    Add-Content $sweepLog $done

    # SECOND GATE PASS, WITH --run-dir. generality_sweep.ps1 invokes the gate without it,
    # which makes `threshold` and `budget` report NOT-APPLICABLE - they are not static
    # checks: distinct measurements per (method, seed) live in round_metrics.csv and the
    # qualification bar lives in the manifest. Four of the gate's six checks therefore never
    # ran in the 2026-11-01 sweep, and a NOT-APPLICABLE is not a pass. Running the gate
    # again against the run it just produced is the only way those two checks are exercised
    # at all, and `budget` is precisely the check that would catch an arm of this panel
    # failing to fill its budget - the thing the paired test depends on.
    if ($dir -ne '(none reported)') {
        $postArgs = @('-X','utf8','-m','veritas.preflight','--library',$p.Key,
                      '--budget',$budget,'--initial-size',$p.Initial,
                      '--corpus-dir',$corpusDir,
                      '--run-dir',(Join-Path $logDir $dir),
                      '--arms',$methods) + $ack
        $postRes = Invoke-Captured -Exe $py -Arguments $postArgs -Tag "gate_post_$($p.Key)"
        $postRes.Lines | Out-File (Join-Path $logDir "panel_gate_postrun_$($p.Key)_$stamp.txt") -Encoding utf8
        $head = "   POST-RUN GATE (with --run-dir) exit=$($postRes.ExitCode)"
        Write-Host $head -ForegroundColor Cyan
        Add-Content $sweepLog $head
        ($postRes.Lines | Select-String -Pattern '^(PASS|FAIL|SKIP|NOT-APPLICABLE)|passed,' ) |
            ForEach-Object { Write-Host "   $_"; Add-Content $sweepLog "   $_" }
        if ($postRes.ExitCode -ne 0) {
            $msg = ("!! POST-RUN GATE FAILED for $($p.Key) (exit $($postRes.ExitCode)). An " +
                    "un-passed check must NOT be reported as a pass; see " +
                    "panel_gate_postrun_$($p.Key)_$stamp.txt")
            Write-Host $msg -ForegroundColor Yellow
            Add-Content $sweepLog $msg
            $postFailed += $p.Key
        }
    } else {
        $msg = "!! no run directory reported for $($p.Key); the post-run gate could not run"
        Write-Host $msg -ForegroundColor Yellow
        Add-Content $sweepLog $msg
    }
}

Write-Host ''
Write-Host "################ panel complete $(Get-Date -Format HH:mm:ss) ################"
Add-Content $sweepLog "panel complete $(Get-Date -Format 'o')"
Write-Host "log: $sweepLog"

# Explicit, because without it the script's exit status is whatever $LASTEXITCODE happened
# to hold - which was 255 on a clean dry run. A sweep's exit code gets read by whatever
# launched it, so it must mean something.
if ($postFailed.Count -gt 0) {
    Write-Host ("landscapes whose POST-RUN gate FAILED: $($postFailed -join ', ') - report the " +
                "failing check alongside the table; do not call it a pass") -ForegroundColor Yellow
    Add-Content $sweepLog "FAILED POST-RUN GATE: $($postFailed -join ', ')"
}
if ($failed.Count -gt 0) {
    Write-Host "landscapes whose preflight FAILED: $($failed -join ', ')" -ForegroundColor Red
    Add-Content $sweepLog "FAILED PREFLIGHT: $($failed -join ', ')"
    exit 1
}
if ($postFailed.Count -gt 0) { exit 3 }
exit 0

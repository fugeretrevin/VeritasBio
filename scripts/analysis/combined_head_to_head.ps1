<#
Combined head-to-head: the engine WITH the one-hot GP surrogate against the seven-arm
baseline panel, all nine arms in ONE run folder per landscape, on ONE seed block, from one
shared starting sample per seed.

WHY. Two 48-seed TrpB results exist and they are on DIFFERENT seed blocks, so they cannot be
combined: the baseline panel (branch worktree-agent-a8325f2a02fbcdc90, 493826b) ran seeds
48-95, and the one-hot surrogate arm (worktree-agent-a14c69a80e0afa44a, a779422) ran
--base-seed 0. Stating "engine+one-hot 57.35 beats best-baseline 38.54" would splice them.
An exactly analogous splice produced the WITHDRAWN "6.3x the best baseline" claim, which
divided one arm's score by a comparator fixed at a different budget and whose quoted range
used two DIFFERENT engine arms. The fix is not a caveat; it is to measure the head-to-head.

Pre-registration, written and committed before any arm here ran:
docs/prereg-combined-head-to-head.md.

THIS IS NOT A SECOND HARNESS. The harness is benchmark_panoptes_isolated.py. This is
baseline_panel.ps1's invocation - the same hard guard, the same derived-table
acknowledgements, the same preflight-before-scored-run ordering, the same Invoke-Captured
stream plumbing, the same stamped output names - with the ninth arm and two landscapes. It
is a sibling rather than a flag on baseline_panel.ps1 because that script is the record of
what the 2026-11-01 panel actually ran, and editing it would rewrite that record.

NOTHING HERE PROMOTES ANYTHING. DEFAULT_METHODS is unchanged (pinned by
tests/test_baseline_panel_arms.py) and `ridge` remains the registered reference baseline;
docs/HANDOFF.md 5.4 makes promoting another model to reference status an amendment to a
frozen hypothesis, which needs a dated amendment record and not a script.

HARD GUARD: refuses to start if a benchmark is already running. Other sessions run scored
campaigns in this same shared results root. An earlier night had two sweeps execute
concurrently for three hours after a task-stop ended a task but left its shell alive; that
voided ten run directories through CPU contention and a mid-run source change. WAITING is
the caller's job, not this script's - queuing behind a stale run nobody knew was alive is
how that happened. A task-stop is NOT a process kill.

Rule: read progress from the run you started, BY MANIFEST, never from "the newest
directory". Each landscape's run directory is recorded into the sweep log as it goes, and a
run directory with no manifest.json is a PARTIAL LAUNCH, not data.
#>
[CmdletBinding()]
param(
    # Landscape keys to run, in the plan's order. Default is both.
    [string[]] $Only,
    # Calibration only: the campaign is registered at 48 paired seeds.
    [int] $Seeds = 48,
    [switch] $WhatIfOnly
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

# --- interpreter and results root, resolved for a WORKTREE as well as the real checkout ---
# In a git worktree $repo/.venv does not exist: it is gitignored and lives only in the real
# repository, which a worktree CANNOT SEE (rule 2.8, three lost agent runs). The source tree
# stays the worktree's - that is the code under test, and the harness hashes it - while the
# interpreter and the results root come from the real repository root. Every failure NAMES
# EVERY PATH IT CHECKED.
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

# Artifacts go to the REAL repository's results root, which is where scripts/analysis reads
# them from. Writing them into the worktree leaves the read-out tools pointing at an empty
# directory while the run looks successful. The path is passed EXPLICITLY, as a constraint.
$logDir = Join-Path $realRoot 'benchmark_comparison_results'
if (-not (Test-Path $logDir)) { throw ("results root absent. Checked: $logDir") }

# The pinned CombinGym corpus the gate's corpus.lock check verifies against. It lives in the
# gitignored downloads/ at the REAL root, so from a worktree the bare gate reports
# "FAILED: corpus.lock - 29 of 29 file(s) absent" and every landscape looks un-quotable when
# in fact the corpus verifies 29/29 at c5dccd514a9a. That is an ENVIRONMENT ARTIFACT, not
# rule-2.5 drift: scripts/freeze_cross_protein_corpus.py resolves downloads/combingym against
# its own __file__. So the path is passed explicitly and every candidate is named on failure.
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

# --- the nine arms ------------------------------------------------------------------------
# COMBINED_HEAD_TO_HEAD_METHODS in benchmark_panoptes_isolated.py. `random` is NOT optional:
# without it nothing says whether a landscape's top-1% count can discriminate methods at all
# (the HER2 informativeness lesson applied before the fact rather than after it). `ridge` is
# not optional either: it is the registered reference and the only arm that makes this run
# comparable to the figures already published here.
$methods = 'random,ridge,random_forest,gbm_onehot,knn_hamming,gp_onehot,greedy_measured,panoptes_default,panoptes_onehot_surrogate'

# --- hard guard ---------------------------------------------------------------------------
$running = @(Get-CimInstance Win32_Process -Filter "Name='python.exe'" |
             Where-Object { $_.CommandLine -and $_.CommandLine -match 'benchmark_panoptes' })
if ($running.Count -ne 0) {
    Write-Host "REFUSING TO START: $($running.Count) benchmark process(es) already running." -ForegroundColor Red
    $running | ForEach-Object { Write-Host "  PID $($_.ProcessId): $($_.CommandLine)" }
    Write-Host "Another session may own these. WAIT for them; do not kill anything you did not start." -ForegroundColor Yellow
    exit 1
}
Write-Host "no benchmark process running; starting campaign at $(Get-Date -Format HH:mm:ss) local"

# --- derived-table acknowledgements -------------------------------------------------------
# Carried over unchanged from baseline_panel.ps1 / generality_sweep.ps1 so the provenance
# posture is identical to the runs this campaign is compared against.
$why = 'reviewed 2026-10-31 by measuring observed-vs-expected coverage per depth; no depth at zero, deep-half coverage within 4.5 points of shallow, unlike the CR9114 clean workbook which fell to zero at depth 16'
$ack = @()
foreach ($k in @('cr6261_h1','cr6261_h9','creilov','eqfp611_blue','eqfp611_combined','eqfp611_red','gb1_combingym')) {
    $ack += '--acknowledge-derived'
    $ack += "$k=$why"
}

# --- landscapes, in the pre-registered priority order -------------------------------------
# trpb FIRST and on seeds 48-95. TrpB is the ONE landscape with a designated held-out block:
# those seeds were never looked at during the 2^5 configuration scan
# (scripts/analysis/heldout.py), and TrpB is the one landscape where the corrected
# qualification bar accepts the identical variant set, so nothing in the 2026-10-31 fix set
# changed what is being searched there.
#
# gb1 SECOND, --base-seed 0. GB1 has NO designated held-out block, and base seed 0 is what
# every GB1 figure in docs/HANDOFF.md was measured on - the engine's -12.90 loss to ridge and
# the one-hot arm's +4.31 recovery. Matching it is what makes this run settle that question
# rather than open a third incomparable one. Its seed block is therefore NOT held out and the
# read-out must say so.
#
# BOTH are DEVELOPMENT DATA under rule 2.4. Nothing produced here is generalization evidence.
$plan = @(
    [pscustomobject]@{ Key='trpb'; Initial=25; Batch=25; BaseSeed=48; Why='HELD-OUT seeds 48-95; the one landscape with a designated held-out block' }
    [pscustomobject]@{ Key='gb1';  Initial=25; Batch=25; BaseSeed=0;  Why='base seed 0, matching every published GB1 figure; NOT a held-out block' }
)
if ($Only) { $plan = $plan | Where-Object { $Only -contains $_.Key } }
if (-not $plan) { throw "no landscapes selected; -Only must name keys from the plan" }

# Every file this script writes carries the same start stamp. Without it the per-landscape
# gate and stdout files have fixed names, and SEVERAL SESSIONS SHARE THIS RESULTS ROOT, so a
# fixed name is a shared mutable name: the panel once aborted on "Access to the path
# panel_gate_trpb.txt is denied" after an interrupted attempt left the file open, which looks
# like a permissions problem and is not one. A collision can also truncate an artifact.
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$sweepLog = Join-Path $logDir ("combined_h2h_{0}.log" -f $stamp)
"combined head-to-head started $(Get-Date -Format 'o')  seeds=$Seeds  methods=$methods" | Out-File $sweepLog -Encoding utf8
"source tree $repo" | Add-Content $sweepLog
Write-Host "campaign log: $sweepLog"

$failed = @()
$postFailed = @()

function Invoke-Captured {
    <#
      Run a native command, returning @{ Lines; ExitCode } with stdout and stderr merged.

      Deliberately NOT `& $py ... 2>&1` under 'Stop'. Windows PowerShell 5.1 wraps each
      stderr line from a native executable in a NativeCommandError ErrorRecord, which under
      $ErrorActionPreference='Stop' throws on the first line the tool logs. Both the
      benchmark and the gate log at INFO to stderr, so a whole sweep once aborted seven
      seconds in, after a gate that had passed.

      The call operator with splatting, NOT Start-Process -ArgumentList: Start-Process joins
      the argument array into one string and lets the child re-parse it, which word-splits
      any argument containing spaces. The --acknowledge-derived KEY=WHY reasons are full
      sentences, so every one arrived as a dozen stray positional arguments.
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
    $gateRes.Lines | Out-File (Join-Path $logDir "combined_h2h_gate_$($p.Key)_$stamp.txt") -Encoding utf8
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

    $runOut | Out-File (Join-Path $logDir "combined_h2h_$($p.Key)_stdout_$stamp.txt") -Encoding utf8
    $runOut | Select-String -Pattern 'Run state|Artifacts|DEGENERATE RUN|^\s*(random|ridge|random_forest|gbm_onehot|knn_hamming|gp_onehot|greedy_measured|panoptes)|^\d+ round' |
        ForEach-Object { Write-Host "   $_"; Add-Content $sweepLog "   $_" }
    $hit = $runOut | Select-String -Pattern 'run_\d{8}_\d{6}_[0-9a-f]+' | Select-Object -First 1
    $dir = if ($hit) { $hit.Matches[0].Value } else { '(none reported)' }
    $done = "## $($p.Key) done $(Get-Date -Format HH:mm:ss) exit=$code run_dir=$dir"
    Write-Host $done
    Add-Content $sweepLog $done

    # SECOND GATE PASS, WITH --run-dir. Without it `threshold` and `budget` report
    # NOT-APPLICABLE - they are not static checks: the distinct measurements per
    # (method, seed) live in round_metrics.csv and the qualification bar lives in the
    # manifest. Four of the gate's six checks therefore never ran in the 2026-11-01 sweep,
    # and a NOT-APPLICABLE IS NOT A PASS. `budget` is precisely the check that catches an arm
    # failing to fill its budget - the thing every paired test here depends on.
    if ($dir -ne '(none reported)') {
        $postArgs = @('-X','utf8','-m','veritas.preflight','--library',$p.Key,
                      '--budget',$budget,'--initial-size',$p.Initial,
                      '--corpus-dir',$corpusDir,
                      '--run-dir',(Join-Path $logDir $dir),
                      '--arms',$methods) + $ack
        $postRes = Invoke-Captured -Exe $py -Arguments $postArgs -Tag "gate_post_$($p.Key)"
        $postRes.Lines | Out-File (Join-Path $logDir "combined_h2h_gate_postrun_$($p.Key)_$stamp.txt") -Encoding utf8
        $head = "   POST-RUN GATE (with --run-dir) exit=$($postRes.ExitCode)"
        Write-Host $head -ForegroundColor Cyan
        Add-Content $sweepLog $head
        ($postRes.Lines | Select-String -Pattern '^(PASS|FAIL|SKIP|NOT-APPLICABLE)|passed,' ) |
            ForEach-Object { Write-Host "   $_"; Add-Content $sweepLog "   $_" }
        if ($postRes.ExitCode -ne 0) {
            $msg = ("!! POST-RUN GATE FAILED for $($p.Key) (exit $($postRes.ExitCode)). An " +
                    "un-passed check must NOT be reported as a pass; see " +
                    "combined_h2h_gate_postrun_$($p.Key)_$stamp.txt")
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
Write-Host "################ campaign complete $(Get-Date -Format HH:mm:ss) ################"
Add-Content $sweepLog "campaign complete $(Get-Date -Format 'o')"
Write-Host "log: $sweepLog"

# Explicit, because without it the script's exit status is whatever $LASTEXITCODE happened to
# hold - which was 255 on a clean dry run once. A sweep's exit code gets read.
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

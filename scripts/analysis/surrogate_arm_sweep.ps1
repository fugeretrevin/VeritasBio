<#
Runs the pre-registered one-hot surrogate arm sweep
(docs/prereg-onehot-featurisation-arm.md).

This is generality_sweep.ps1's invocation with exactly two things changed:

  * `-Methods` is a parameter, because the registered arm `panoptes_onehot_surrogate`
    has to run in the SAME run folder as `panoptes_default` - both arms then see one
    shared starting sample per seed and every comparison is paired WITHIN seed instead
    of across run folders. generality_sweep.ps1 hardcodes three methods.
  * the plan carries creilov, gb1 and trpb, which is what the pre-registration names.

Everything else is deliberately copied rather than rewritten: the same budgets, the same
seed block, the same provenance mode, the same seven --acknowledge-derived reasons, the
same pre-flight gate at the budget about to run, the same hard guard against a concurrent
benchmark, and the same Invoke-Captured helper. The .sh sweep cannot run on this machine
(no Git bash, no WSL), which is why the .ps1 exists at all.

Three notes carried over from generality_sweep.ps1 because each cost a real run:

  * HARD GUARD. Two sweeps once executed concurrently for three hours after a task-stop
    ended the task but left its shell alive, voiding ten run directories through CPU
    contention and a mid-run source change. A task-stop is not a process kill.
  * Read progress from the run you started, BY MANIFEST, never from "the newest
    directory". During that contention TrpB was reported as CreiLOV for an hour. Each
    landscape's own run directory is therefore recorded into the sweep log as it goes.
  * NOT `& $py ... 2>&1` through PowerShell's own plumbing: 5.1 wraps each stderr line
    from a native exe in a NativeCommandError, which under $ErrorActionPreference='Stop'
    throws on the first INFO line the tool logs.
#>
[CmdletBinding()]
param(
    # Landscape keys to run, in order. Default is the full pre-registered set.
    [string[]] $Only,
    # The arms. The default is the pre-registered four.
    [string] $Methods = 'random,ridge,panoptes_default,panoptes_onehot_surrogate',
    # Calibration: override the seed count (the study is registered at 48).
    [int] $Seeds = 48,
    # A label written into the sweep log, so a calibration sweep cannot be mistaken for
    # the registered one when both are on disk.
    [string] $Tag = 'prereg',
    [switch] $WhatIfOnly
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repo = Resolve-Path (Join-Path $PSScriptRoot '..\..')
Set-Location $repo
# The interpreter lives in the REAL repository root, not in a worktree (rule 2.8). Both
# paths are named when neither resolves, because a bare "not found" is what let our
# strongest baseline get silently dropped from a 48-seed run.
$checked = @(Join-Path $repo '.venv\Scripts\python.exe')
$common = & git -C $repo rev-parse --path-format=absolute --git-common-dir 2>$null
if ($common) { $checked += (Join-Path (Split-Path -Parent $common) '.venv\Scripts\python.exe') }
$py = $checked | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $py) { throw "no interpreter found. Checked: $($checked -join '; ')" }
Write-Host "interpreter: $py"

# --- hard guard ------------------------------------------------------------------------
$running = @(Get-CimInstance Win32_Process -Filter "Name='python.exe'" |
             Where-Object { $_.CommandLine -and $_.CommandLine -match 'benchmark_panoptes' })
if ($running.Count -ne 0) {
    Write-Host "REFUSING TO START: $($running.Count) benchmark process(es) already running." -ForegroundColor Red
    $running | ForEach-Object { Write-Host "  PID $($_.ProcessId): $($_.CommandLine)" }
    Write-Host "Tree-kill them first:  taskkill /F /T /PID <pid>" -ForegroundColor Yellow
    exit 1
}
Write-Host "no benchmark process running; starting sweep at $(Get-Date -Format HH:mm:ss) local"

# --- derived-table acknowledgements ----------------------------------------------------
# Identical to generality_sweep.ps1's, including the reason text, so the gate decisions
# are the same ones that sweep recorded.
$why = 'reviewed 2026-10-31 by measuring observed-vs-expected coverage per depth; no depth at zero, deep-half coverage within 4.5 points of shallow, unlike the CR9114 clean workbook which fell to zero at depth 16'
$ack = @()
foreach ($k in @('cr6261_h1','cr6261_h9','creilov','eqfp611_blue','eqfp611_combined','eqfp611_red','gb1_combingym')) {
    $ack += '--acknowledge-derived'
    $ack += "$k=$why"
}

# --- the plan --------------------------------------------------------------------------
# The two landscapes the engine loses worst on, plus TrpB as the negative control. Budgets
# are generality_sweep.ps1's, unchanged. CreiLOV is both the most informative and the
# slowest, so it is not last.
$plan = @(
    [pscustomobject]@{ Key='creilov'; Initial=25; Batch=25; Why='15 sites mixed alphabet - engine loses by -50.29 on 0/48 seeds; the headline test' }
    [pscustomobject]@{ Key='gb1';     Initial=25; Batch=25; Why='saturated 4-site - engine loses by -12.90 on 6/48; biophysical features win here under a LINEAR readout' }
    [pscustomobject]@{ Key='trpb';    Initial=25; Batch=25; Why='NEGATIVE CONTROL - the one-hot arm is a NEW ARM here, never a correction to the published TrpB figure' }
)
if ($Only) { $plan = $plan | Where-Object { $Only -contains $_.Key } }
if (-not $plan) { throw "no landscapes selected; -Only must name keys from the plan" }

# The results directory resolves against the REAL repository root, not the worktree this
# may be standing in: it is gitignored, so a worktree cannot see it and a relative path
# would silently create a second, empty one beside the real artifacts (rule 2.8). Asked of
# veritas.paths rather than assembled here.
$logDir = (& $py -X utf8 -c "import veritas.paths as p; print(p.results_dir('benchmark_comparison_results'))")
if (-not $logDir -or -not (Test-Path $logDir)) {
    throw ("could not resolve the results root. veritas.paths.results_dir" +
           "('benchmark_comparison_results') was asked from $repo and returned '$logDir'. " +
           "It resolves against the real repository root, never the worktree.")
}
$logDir = (Resolve-Path $logDir).Path
Write-Host "results root: $logDir"
$sweepLog = Join-Path $logDir ("surrogate_arm_sweep_{0}_{1}.log" -f $Tag, (Get-Date -Format 'yyyyMMdd_HHmmss'))
"sweep started $(Get-Date -Format 'o')  tag=$Tag  seeds=$Seeds  methods=$Methods" | Out-File $sweepLog -Encoding utf8
"head: $(& git -C $repo rev-parse HEAD)" | Add-Content $sweepLog
"status: $((& git -C $repo status --porcelain) -join ' | ')" | Add-Content $sweepLog
Write-Host "sweep log: $sweepLog"

$failed = @()

function Invoke-Captured {
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
    $banner = "## $($p.Key)   initial=$($p.Initial) batch=$($p.Batch) budget=$budget   start $(Get-Date -Format HH:mm:ss)"
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
                  '--budget',$budget,'--initial-size',$p.Initial) + $ack
    $gateRes = Invoke-Captured -Exe $py -Arguments $gateArgs -Tag "gate_$($p.Key)"
    ($gateRes.Lines | Select-Object -Last 4) | ForEach-Object { Write-Host "   $_"; Add-Content $sweepLog "   $_" }
    if ($gateRes.ExitCode -ne 0) {
        $msg = "!! PREFLIGHT FAILED for $($p.Key) (exit $($gateRes.ExitCode)) - DO NOT QUOTE THIS LANDSCAPE"
        Write-Host $msg -ForegroundColor Red
        Add-Content $sweepLog $msg
        $failed += $p.Key
    } else {
        Write-Host "   preflight OK for $($p.Key)" -ForegroundColor Green
        Add-Content $sweepLog "   preflight OK for $($p.Key)"
    }

    $runArgs = @('-X','utf8','benchmark_panoptes_isolated.py',
                 '--library',$p.Key,'--methods',$Methods,
                 '--initial-size',$p.Initial,'--batch-size',$p.Batch,
                 '--seeds',$Seeds,'--base-seed','0','--rounds','4',
                 '--provenance-mode','acknowledge-unflagged')
    $runRes = Invoke-Captured -Exe $py -Arguments $runArgs -Tag "run_$($p.Key)"
    $runOut = $runRes.Lines

    $runOut | Out-File (Join-Path $logDir "surrogate_arm_$($Tag)_$($p.Key)_stdout.txt") -Encoding utf8
    $runOut | Select-String -Pattern 'Run state|Artifacts|DEGENERATE RUN|^\s*(random|ridge|panoptes)|^\d+ round' |
        ForEach-Object { Write-Host "   $_"; Add-Content $sweepLog "   $_" }
    $hit = $runOut | Select-String -Pattern 'run_\d{8}_\d{6}_[0-9a-f]+' | Select-Object -First 1
    $dir = if ($hit) { $hit.Matches[0].Value } else { '(none reported)' }
    $done = "## $($p.Key) done $(Get-Date -Format HH:mm:ss) exit=$($runRes.ExitCode) run_dir=$dir"
    Write-Host $done
    Add-Content $sweepLog $done
}

Write-Host ''
Write-Host "################ sweep complete $(Get-Date -Format HH:mm:ss) ################"
Add-Content $sweepLog "sweep complete $(Get-Date -Format 'o')"
Write-Host "log: $sweepLog"

if ($failed.Count -gt 0) {
    Write-Host "landscapes whose preflight FAILED: $($failed -join ', ')" -ForegroundColor Red
    Add-Content $sweepLog "FAILED PREFLIGHT: $($failed -join ', ')"
    exit 1
}
exit 0

<#
Spatial viability veto ablation: is the veto the CAUSE of the engine's three measured
losses against one-hot ridge, or merely correlated with them?

Pre-registered in docs/prereg-spatial-veto-ablation.md, committed before this ran.

A direct port of generality_sweep.ps1's invocation - same interpreter resolution, same
hard guard, same derived-table acknowledgements, same preflight-before-scored-run
ordering, same Invoke-Captured stream plumbing, same per-landscape run-directory
recording. The ONE difference is that each landscape runs TWICE:

    run A: --methods ridge,panoptes_default                                (unablated)
    run B: --methods ridge,panoptes_default --ablate-components spatial_filter

because --ablate-components is a single global flag applied to every panoptes arm in a
run, so the two engine arms cannot share a run directory. The pairing survives it:
round 0 is `random.Random(seed).sample over the eligible set, identical for every method
within a seed`, so the same seed block gives the same initial sample in both runs.
`ridge` is deliberately carried in BOTH runs as the pairing check - if ridge is not
bit-identical across a landscape's two runs, the runs are not paired and the landscape
must not be quoted.

Interpreter: resolved against the REAL repository root, never the worktree. A linked
worktree has no .venv (HANDOFF rule 2.8; it has cost three agent runs).
#>
[CmdletBinding()]
param(
    [string[]] $Only,
    [int] $Seeds = 48,
    [switch] $WhatIfOnly
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$checkout = Resolve-Path (Join-Path $PSScriptRoot '..\..')
Set-Location $checkout

# HANDOFF rule 2.8 / docs/agents/paths.md: the interpreter and the gitignored results
# root live at the REAL root. Every failed resolution names the paths it checked.
$real = (& git -C $checkout rev-parse --path-format=absolute --git-common-dir) -replace '[\\/]\.git[\\/]?$', ''
if (-not $real) { throw "could not resolve the real repository root from $checkout" }
$py = Join-Path $real '.venv\Scripts\python.exe'
if (-not (Test-Path $py)) {
    throw ("no interpreter found. Checked: {0} (checkout {1}); {2} (real root {3}). " -f
           (Join-Path $checkout '.venv\Scripts\python.exe'), $checkout, $py, $real)
}
Write-Host "checkout : $checkout"
Write-Host "real root: $real"
Write-Host "python   : $py"

# The worktree's own sources must be what runs, not the real root's.
$env:PYTHONPATH = $checkout

# --- hard guard (verbatim from generality_sweep.ps1) -----------------------------------
$running = @(Get-CimInstance Win32_Process -Filter "Name='python.exe'" |
             Where-Object { $_.CommandLine -and $_.CommandLine -match 'benchmark_panoptes' })
if ($running.Count -ne 0) {
    Write-Host "REFUSING TO START: $($running.Count) benchmark process(es) already running." -ForegroundColor Red
    $running | ForEach-Object { Write-Host "  PID $($_.ProcessId): $($_.CommandLine)" }
    Write-Host "Tree-kill them first:  taskkill /F /T /PID <pid>" -ForegroundColor Yellow
    exit 1
}
Write-Host "no benchmark process running; starting at $(Get-Date -Format HH:mm:ss) local"

# --- derived-table acknowledgements (verbatim from generality_sweep.ps1) ---------------
$why = 'reviewed 2026-10-31 by measuring observed-vs-expected coverage per depth; no depth at zero, deep-half coverage within 4.5 points of shallow, unlike the CR9114 clean workbook which fell to zero at depth 16'
$ack = @()
foreach ($k in @('cr6261_h1','cr6261_h9','creilov','eqfp611_blue','eqfp611_combined','eqfp611_red','gb1_combingym')) {
    $ack += '--acknowledge-derived'
    $ack += "$k=$why"
}

# Budgets are generality_sweep.ps1's registered plan, unchanged. Priority order is the
# pre-registration's: red first (its ceiling is exactly 0, so it is the decisive case),
# then CreiLOV (the largest loss), then GB1.
$plan = @(
    [pscustomobject]@{ Key='eqfp611_red'; Initial=17; Batch=16; Why='veto ceiling is 0 of 82 hits - the unablated arm CANNOT score' }
    [pscustomobject]@{ Key='creilov';     Initial=25; Batch=25; Why='largest measured loss (-50.29); veto ceiling 250 of 1,655' }
    [pscustomobject]@{ Key='gb1';         Initial=25; Batch=25; Why='-12.90; cohesion rule never fires here - core/steric only' }
)
if ($Only) { $plan = $plan | Where-Object { $Only -contains $_.Key } }
if (-not $plan) { throw "no landscapes selected; -Only must name keys from the plan" }

# The benchmark writes its run folders beside the sources it ran, i.e. into the CHECKOUT's
# gitignored benchmark_comparison_results (verified by a one-seed probe). The campaign log
# goes to the same place so a landscape's log and its run directories cannot end up in
# two different trees. The interpreter above still comes from the real root; those are
# different questions and rule 2.8's point is to keep them apart.
$logDir = Join-Path $checkout 'benchmark_comparison_results'
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$log = Join-Path $logDir ("vetoablation_{0}.log" -f (Get-Date -Format 'yyyyMMdd_HHmmss'))
"veto-ablation campaign started $(Get-Date -Format 'o')  seeds=$Seeds" | Out-File $log -Encoding utf8
Write-Host "log: $log"

$failed = @()

function Invoke-Captured {
    # Deliberately NOT `& $py ... 2>&1` under $ErrorActionPreference='Stop': 5.1 wraps each
    # stderr line from a native exe in a NativeCommandError, and both the benchmark and the
    # gate log at INFO to stderr. See generality_sweep.ps1 for the full story.
    param([string] $Exe, [string[]] $Arguments)
    $prev = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try {
        $raw = & $Exe @Arguments 2>&1
        $code = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $prev
    }
    return @{ Lines = @($raw | ForEach-Object { [string] $_ }); ExitCode = $code }
}

foreach ($p in $plan) {
    $budget = $p.Initial + 4 * $p.Batch
    $banner = "## $($p.Key)   initial=$($p.Initial) batch=$($p.Batch) budget=$budget   start $(Get-Date -Format HH:mm:ss)"
    Write-Host ''
    Write-Host '################################################################'
    Write-Host $banner
    Write-Host "## $($p.Why)"
    Write-Host '################################################################'
    Add-Content $log "`n$banner`n## $($p.Why)"

    if ($WhatIfOnly) { Write-Host "   (WhatIfOnly) would gate + run two arms for $($p.Key)"; continue }

    # Rule 2.7: the gate runs before any scored run, at the budget actually about to run.
    $gateArgs = @('-X','utf8','-m','veritas.preflight','--library',$p.Key,
                  '--budget',$budget,'--initial-size',$p.Initial) + $ack
    $gateRes = Invoke-Captured -Exe $py -Arguments $gateArgs
    ($gateRes.Lines | Select-Object -Last 4) | ForEach-Object { Write-Host "   $_"; Add-Content $log "   $_" }
    if ($gateRes.ExitCode -ne 0) {
        $msg = "!! PREFLIGHT FAILED for $($p.Key) (exit $($gateRes.ExitCode)) - DO NOT QUOTE THIS LANDSCAPE"
        Write-Host $msg -ForegroundColor Red; Add-Content $log $msg
        $failed += $p.Key
    } else {
        Write-Host "   preflight OK for $($p.Key)" -ForegroundColor Green
        Add-Content $log "   preflight OK for $($p.Key)"
    }

    foreach ($arm in @(
        [pscustomobject]@{ Tag='unablated'; Extra=@() }
        [pscustomobject]@{ Tag='ablated';   Extra=@('--ablate-components','spatial_filter') }
    )) {
        Write-Host "   --- arm: $($arm.Tag)" -ForegroundColor Cyan
        $runArgs = @('-X','utf8','benchmark_panoptes_isolated.py',
                     '--library',$p.Key,'--methods','ridge,panoptes_default',
                     '--initial-size',$p.Initial,'--batch-size',$p.Batch,
                     '--seeds',$Seeds,'--base-seed','0','--rounds','4',
                     '--provenance-mode','acknowledge-unflagged') + $arm.Extra
        $runRes = Invoke-Captured -Exe $py -Arguments $runArgs
        $out = $runRes.Lines
        $out | Out-File (Join-Path $logDir "vetoablation_$($p.Key)_$($arm.Tag)_stdout.txt") -Encoding utf8
        $out | Select-String -Pattern 'Run state|Artifacts|DEGENERATE RUN|^\s*(ridge|panoptes)|^\d+ round' |
            ForEach-Object { Write-Host "      $_"; Add-Content $log "      $_" }
        $hit = $out | Select-String -Pattern 'run_\d{8}_\d{6}_[0-9a-f]+' | Select-Object -First 1
        $dir = if ($hit) { $hit.Matches[0].Value } else { '(none reported)' }
        $done = "## $($p.Key) [$($arm.Tag)] done $(Get-Date -Format HH:mm:ss) exit=$($runRes.ExitCode) run_dir=$dir"
        Write-Host "   $done"; Add-Content $log $done
    }
}

Write-Host ''
Write-Host "################ campaign complete $(Get-Date -Format HH:mm:ss) ################"
Add-Content $log "campaign complete $(Get-Date -Format 'o')"
Write-Host "log: $log"
if ($failed.Count -gt 0) {
    Write-Host "landscapes whose preflight FAILED: $($failed -join ', ')" -ForegroundColor Red
    Add-Content $log "FAILED PREFLIGHT: $($failed -join ', ')"
    exit 1
}
exit 0

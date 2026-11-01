<#
Generality sweep: does library-aware candidate generation help where reachability was the
bottleneck, and leave the saturated libraries alone?

Pre-registered prediction (docs/day1-reachability.md section 7a): improvement on the
reachability-limited landscapes, and NO CHANGE on GB1, which is 100% reachable at
depth <= 4 and so has nothing for these fixes to add.

A PowerShell port of generality_sweep.sh, which cannot run on the primary Windows
development machine: it has neither Git bash nor a WSL distribution. That is not a
cosmetic difference. The sweep the .sh script describes produced run directories that are
absent from this checkout, so its results cannot be re-verified here, and a sweep that
can only run elsewhere is a sweep whose evidence this machine cannot audit. Both scripts
are kept because the .sh is the record of what was run before.

Behaviour is intended to match the .sh exactly: same landscape order, budgets, methods,
seed block, provenance mode and acknowledgements.

HARD GUARD: refuses to start if a benchmark is already running. An earlier attempt had
two sweeps execute concurrently for three hours after a task-stop ended the task but left
its shell alive; that voided ten run directories through CPU contention, a mixed code
state and a mid-run source change (docs/day1-reachability.md section 11). Waiting was the
wrong behaviour - it let a second sweep queue up behind a stale one nobody knew was alive.

Rule 3 of that section: read progress from the run you started, BY MANIFEST, never from
"the newest directory". During contention the newest was the other sweep's, and TrpB was
reported as CreiLOV for an hour. This script therefore records each landscape's own run
directory into the sweep log as it goes.
#>
[CmdletBinding()]
param(
    # Landscape keys to run, in order. Default is the full pre-registered sweep.
    [string[]] $Only,
    # Calibration: override the seed count (the sweep is registered at 48).
    [int] $Seeds = 48,
    [switch] $WhatIfOnly
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repo = Resolve-Path (Join-Path $PSScriptRoot '..\..')
Set-Location $repo
$py = Join-Path $repo '.venv\Scripts\python.exe'
if (-not (Test-Path $py)) { throw "no interpreter at $py" }

$methods = 'random,ridge,panoptes_default'

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
# corpus.derived_table flags seven landscapes reading a CombinGym "clean" workbook while a
# raw table exists at the pin. Reviewed by measuring observed-vs-expected genotype coverage
# PER DEPTH on all nine landscapes (scripts/analysis/depth_coverage.py): no depth is at
# zero anywhere, eqFP611 and TrpB are complete, CR6261-H9's deep half is slightly higher
# than its shallow half, and the largest gradient is CreiLOV at 4.5 points - none shows the
# CR9114 pattern of attrition reaching 100% at maximum depth. See section 10.
$why = 'reviewed 2026-10-31 by measuring observed-vs-expected coverage per depth; no depth at zero, deep-half coverage within 4.5 points of shallow, unlike the CR9114 clean workbook which fell to zero at depth 16'
$ack = @()
foreach ($k in @('cr6261_h1','cr6261_h9','creilov','eqfp611_blue','eqfp611_combined','eqfp611_red','gb1_combingym')) {
    $ack += '--acknowledge-derived'
    $ack += "$k=$why"
}

# --- the sweep -------------------------------------------------------------------------
# Scope trimmed to what earns its runtime. Dropped: TrpB (already measured bit-identical on
# code carrying every generator fix, and it declares no per-site alphabet so the encoding
# change cannot touch it); CR6261 H1/H9 (19-variant top-1% tier, low power by construction,
# and one protein against two antigens so not independent of each other).
#
# Fast landscapes first, so signal arrives early, with CreiLOV third rather than last - it
# is both the most informative and the slowest, and should not be the thing lost if the
# sweep is interrupted.
$plan = @(
    [pscustomobject]@{ Key='eqfp611_red';  Initial=17; Batch=16; Why='13 sites binary, NON-functional parent, prior run found 0 hits' }
    [pscustomobject]@{ Key='eqfp611_blue'; Initial=17; Batch=16; Why='same geometry as red, functional parent - paired control' }
    [pscustomobject]@{ Key='creilov';      Initial=25; Batch=25; Why='15 sites, mixed alphabet - highest-value generality test' }
    [pscustomobject]@{ Key='gb1';          Initial=25; Batch=25; Why='NEGATIVE CONTROL - saturated 4-site, must not change' }
    [pscustomobject]@{ Key='cr9114_h1';    Initial=25; Batch=25; Why='16 sites binary - DEVELOPMENT DATA, fair-encoding re-run' }
    [pscustomobject]@{ Key='her2';         Initial=25; Batch=22; Why='NULL CONTROL - pool_restricted, never invokes the generator' }
)
if ($Only) { $plan = $plan | Where-Object { $Only -contains $_.Key } }
if (-not $plan) { throw "no landscapes selected; -Only must name keys from the plan" }

$logDir = Join-Path $repo 'benchmark_comparison_results'
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$sweepLog = Join-Path $logDir ("sweep_{0}.log" -f (Get-Date -Format 'yyyyMMdd_HHmmss'))
"sweep started $(Get-Date -Format 'o')  seeds=$Seeds  methods=$methods" | Out-File $sweepLog -Encoding utf8
Write-Host "sweep log: $sweepLog"

# Landscapes whose gate failed. Collected rather than aborting: a failed gate means that
# landscape must not be quoted (rule 2.7), not that the rest of the sweep is worthless.
$failed = @()

function Invoke-Captured {
    <#
      Run a native command, returning @{ Lines; ExitCode } with stdout and stderr merged.

      Deliberately NOT `& $py ... 2>&1`. Windows PowerShell 5.1 wraps each stderr line
      from a native executable in a NativeCommandError ErrorRecord, which under
      $ErrorActionPreference='Stop' throws on the first line the tool logs. Both the
      benchmark and the preflight gate log at INFO to stderr, so the whole sweep aborted
      seven seconds in, after a gate that had passed. Redirecting the child process's own
      handles to files sidesteps PowerShell's stream plumbing entirely.
    #>
    param([string] $Exe, [string[]] $Arguments, [string] $Tag)

    # The call operator with splatting, NOT Start-Process -ArgumentList. Start-Process
    # joins the argument array into one string and lets the child re-parse it, which
    # word-splits any argument containing spaces: the --acknowledge-derived KEY=WHY
    # reasons are full sentences, so every one of the seven arrived as a dozen stray
    # positional arguments and the gate exited 2 with "unrecognized arguments". The .sh
    # version is safe because a bash array quotes each element; -ArgumentList does not.
    #
    # 'Continue' for the duration: 5.1 wraps each stderr line from a native exe in a
    # NativeCommandError, which under 'Stop' throws on the first INFO line logged. Exit
    # status is read from $LASTEXITCODE, which stays correct, rather than from $?.
    $prev = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try {
        $raw = & $Exe @Arguments 2>&1
        $code = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $prev
    }
    # ErrorRecords stringify to their message, so this flattens both streams to text.
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
    $gate = $gateRes.ExitCode
    ($gateRes.Lines | Select-Object -Last 4) | ForEach-Object { Write-Host "   $_"; Add-Content $sweepLog "   $_" }
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
                 '--seeds',$Seeds,'--base-seed','0','--rounds','4',
                 '--provenance-mode','acknowledge-unflagged')
    $runRes = Invoke-Captured -Exe $py -Arguments $runArgs -Tag "run_$($p.Key)"
    $runOut = $runRes.Lines
    $code = $runRes.ExitCode

    # Record the run directory this landscape actually produced, so progress is never read
    # from "the newest directory" (rule 3).
    $runOut | Out-File (Join-Path $logDir "sweep_$($p.Key)_stdout.txt") -Encoding utf8
    $runOut | Select-String -Pattern 'Run state|Artifacts|DEGENERATE RUN|^\s*(random|ridge|panoptes)|^\d+ round' |
        ForEach-Object { Write-Host "   $_"; Add-Content $sweepLog "   $_" }
    $hit = $runOut | Select-String -Pattern 'run_\d{8}_\d{6}_[0-9a-f]+' | Select-Object -First 1
    $dir = if ($hit) { $hit.Matches[0].Value } else { '(none reported)' }
    $done = "## $($p.Key) done $(Get-Date -Format HH:mm:ss) exit=$code run_dir=$dir"
    Write-Host $done
    Add-Content $sweepLog $done
}

Write-Host ''
Write-Host "################ sweep complete $(Get-Date -Format HH:mm:ss) ################"
Add-Content $sweepLog "sweep complete $(Get-Date -Format 'o')"
Write-Host "log: $sweepLog"

# Explicit, because without it the script's exit status is whatever $LASTEXITCODE happened
# to hold - which was 255 on a clean dry run. A sweep's exit code gets read by whatever
# launched it, so it must mean something.
if ($failed.Count -gt 0) {
    Write-Host "landscapes whose preflight FAILED: $($failed -join ', ')" -ForegroundColor Red
    Add-Content $sweepLog "FAILED PREFLIGHT: $($failed -join ', ')"
    exit 1
}
exit 0

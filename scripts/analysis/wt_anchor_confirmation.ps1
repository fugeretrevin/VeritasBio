<#
Confirm the fabricated WT training row at 48 paired seeds.

Pre-registered in docs/prereg-wt-anchor-confirmation.md, committed BEFORE any arm here
ran. `select_next_micro_batch` injects one synthetic row per round for the wild type with
the literal label y = 0.0 whatever the assay's scale, at variance 1e-4 against a real
measurement's 1e-2. Measured once at 12 seeds on CreiLOV: omitting it is worth +29.75
top-1% hits (t = 12.17, 12 of 12 seeds). This settles it at 48 and measures how far it
generalises.

A port of generality_sweep.ps1's invocation for a two-arm-per-landscape design, following
spatial_veto_ablation.ps1, because --ablate-components is global to a run: each landscape
gets one unablated run and one ablated run, and `ridge` rides along in both so the pairing
can be PROVED rather than assumed.

Three differences from generality_sweep.ps1, each deliberate:

 1. The gate runs WITH --run-dir, after the run rather than before it, as well as before.
    generality_sweep.ps1 invokes it without, which makes both `threshold` and `budget`
    report NOT-APPLICABLE - and that is why nobody noticed that threshold.<lib> FAILS on
    the three runs HANDOFF section 1's published figures come from.

 2. The gate gets --corpus-dir pointed at the REAL root's downloads/combingym.
    freeze_cross_protein_corpus.py:55 resolves the corpus against its own __file__ and
    downloads/ is gitignored, so from a worktree all 29 pinned files report `absent` and
    corpus.lock FAILS as an ENVIRONMENT ARTIFACT, not as rule-2.5 corpus drift. With the
    flag: CORPUS VERIFIED: 29/29 files match CombinGym @ c5dccd514a9a.

 3. Landscape order is the pre-registration's priority order, not the sweep's: creilov
    (the confirmation) first, then trpb and gb1 (the built-in negative controls, and the
    two landscapes every quotable engine figure rests on), then cr9114_h1 (the
    counter-example - it violates the label scale by 14.14 SDs and the engine WINS there).
#>
[CmdletBinding()]
param(
    [string[]] $Only,
    [int] $Seeds = 48,
    [string] $Ablate = 'wt_anchor',
    [string] $Tag = 'wtanchor',
    [switch] $WhatIfOnly
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$checkout = Resolve-Path (Join-Path $PSScriptRoot '..\..')
Set-Location $checkout

# HANDOFF rule 2.8 / docs/agents/paths.md: the interpreter, the gitignored results root
# and the corpus all resolve against the REAL repository root, not this worktree. Every
# failed resolution names the paths it checked.
$real = (& git -C $checkout rev-parse --path-format=absolute --git-common-dir) -replace '[\\/]\.git[\\/]?$', ''
if (-not $real) { throw "could not resolve the real repository root from $checkout" }
$py = Join-Path $real '.venv\Scripts\python.exe'
if (-not (Test-Path $py)) {
    throw ("no interpreter found. Checked: {0} (checkout {1}); {2} (real root {3})." -f
           (Join-Path $checkout '.venv\Scripts\python.exe'), $checkout, $py, $real)
}
$corpus = Join-Path $real 'downloads\combingym'
if (-not (Test-Path $corpus)) {
    throw ("no corpus found. Checked: {0} (real root {1}); {2} (checkout, gitignored and so absent from a worktree)." -f
           $corpus, $real, (Join-Path $checkout 'downloads\combingym'))
}
Write-Host "checkout : $checkout"
Write-Host "real root: $real"
Write-Host "python   : $py"
Write-Host "corpus   : $corpus"

# The worktree's own sources must be what runs, not the real root's.
$env:PYTHONPATH = $checkout

# --- hard guard (verbatim from generality_sweep.ps1) -----------------------------------
# An earlier attempt had two sweeps execute concurrently for three hours after a task-stop
# ended the task but left its shell alive; that voided ten run directories.
$running = @(Get-CimInstance Win32_Process -Filter "Name='python.exe'" |
             Where-Object { $_.CommandLine -and $_.CommandLine -match 'benchmark_panoptes' })
if ($running.Count -ne 0) {
    Write-Host "REFUSING TO START: $($running.Count) benchmark process(es) already running." -ForegroundColor Red
    $running | ForEach-Object { Write-Host "  PID $($_.ProcessId): $($_.CommandLine)" }
    Write-Host "Tree-kill them first:  taskkill /F /T /PID <pid>" -ForegroundColor Yellow
    exit 1
}

# Rule: no repository edits while a scored run executes - the harness hashes the source
# tree and voids the run. Refuse on a dirty tree rather than silently scoring a mixed one.
$dirty = @(& git -C $checkout status --porcelain)
if ($dirty.Count -ne 0) {
    Write-Host "REFUSING TO START: the tree is dirty, and the harness hashes the source tree." -ForegroundColor Red
    $dirty | ForEach-Object { Write-Host "  $_" }
    exit 1
}
Write-Host "no benchmark process running, tree clean; starting at $(Get-Date -Format HH:mm:ss) local"

# --- derived-table acknowledgements (verbatim from generality_sweep.ps1) ---------------
$why = 'reviewed 2026-10-31 by measuring observed-vs-expected coverage per depth; no depth at zero, deep-half coverage within 4.5 points of shallow, unlike the CR9114 clean workbook which fell to zero at depth 16'
$ack = @()
foreach ($k in @('cr6261_h1','cr6261_h9','creilov','eqfp611_blue','eqfp611_combined','eqfp611_red','gb1_combingym')) {
    $ack += '--acknowledge-derived'
    $ack += "$k=$why"
}

# Budgets are generality_sweep.ps1's registered plan, unchanged. TrpB and GB1 run at
# 25 + 4x25 = 125, which is the budget every quotable TrpB figure was measured at.
$plan = @(
    [pscustomobject]@{ Key='creilov';   Initial=25; Batch=25; Why='THE CONFIRMATION - 15.68 label SDs, -50.29 vs ridge, +29.75 at 12 seeds' }
    [pscustomobject]@{ Key='trpb';      Initial=25; Batch=25; Why='NEGATIVE CONTROL - 0.35 label SDs; every quotable engine figure rests here' }
    [pscustomobject]@{ Key='gb1';       Initial=25; Batch=25; Why='NEGATIVE CONTROL - 0.20 label SDs, the least violated landscape' }
    [pscustomobject]@{ Key='cr9114_h1'; Initial=25; Batch=25; Why='THE COUNTER-EXAMPLE - 14.14 SDs violated and the engine WINS (+3.17)' }
    [pscustomobject]@{ Key='eqfp611_red';  Initial=17; Batch=16; Why='secondary - veto ceiling is 0 of 82 hits, so 0.0000 is arithmetic' }
    [pscustomobject]@{ Key='eqfp611_blue'; Initial=17; Batch=16; Why='secondary - 1.19 label SDs, paired geometry control for red' }
)
if ($Only) { $plan = $plan | Where-Object { $Only -contains $_.Key } }
if (-not $plan) { throw "no landscapes selected; -Only must name keys from the plan" }

# The benchmark writes its run folders beside the sources it ran, i.e. into the CHECKOUT's
# gitignored benchmark_comparison_results. The campaign log goes to the same place so a
# landscape's log and its run directories cannot end up in two different trees.
$logDir = Join-Path $checkout 'benchmark_comparison_results'
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$log = Join-Path $logDir ("{0}_{1}.log" -f $Tag, (Get-Date -Format 'yyyyMMdd_HHmmss'))
"wt_anchor confirmation started $(Get-Date -Format 'o')  seeds=$Seeds  ablate=$Ablate  head=$(& git -C $checkout rev-parse HEAD)" |
    Out-File $log -Encoding utf8
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

function Invoke-Gate {
    <# Rule 2.7, WITH --run-dir so `threshold` and `budget` are actually exercised. #>
    param([string] $Key, [int] $Budget, [int] $Initial, [string] $RunDir, [string] $Label)
    $a = @('-X','utf8','-m','veritas.preflight','--library',$Key,
           '--budget',$Budget,'--initial-size',$Initial,'--corpus-dir',$corpus) + $ack
    if ($RunDir) { $a += @('--run-dir', $RunDir) }
    $res = Invoke-Captured -Exe $py -Arguments $a
    Write-Host "   gate [$Label] exit=$($res.ExitCode)" -ForegroundColor Cyan
    Add-Content $log "   gate [$Label] exit=$($res.ExitCode)"
    $res.Lines | Select-String -Pattern '^(PASS|FAIL|NOT-APPLICABLE|SKIP|\d+ passed|FAILED:)' |
        ForEach-Object { Write-Host "      $_"; Add-Content $log "      $_" }
    return $res.ExitCode
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

    # Pre-run gate. Static checks only; `threshold` and `budget` need a campaign, so the
    # authoritative verdict is the post-run gate below.
    if ((Invoke-Gate -Key $p.Key -Budget $budget -Initial $p.Initial -RunDir '' -Label 'pre-run, static') -ne 0) {
        $failed += "$($p.Key):pre-run"
    }

    foreach ($arm in @(
        [pscustomobject]@{ ArmTag='default'; Extra=@() }
        [pscustomobject]@{ ArmTag='ablated'; Extra=@('--ablate-components',$Ablate) }
    )) {
        Write-Host "   --- arm: $($arm.ArmTag)" -ForegroundColor Cyan
        Add-Content $log "   --- arm: $($arm.ArmTag)"
        $runArgs = @('-X','utf8','benchmark_panoptes_isolated.py',
                     '--library',$p.Key,'--methods','ridge,panoptes_default',
                     '--initial-size',$p.Initial,'--batch-size',$p.Batch,
                     '--seeds',$Seeds,'--base-seed','0','--rounds','4',
                     '--provenance-mode','acknowledge-unflagged') + $arm.Extra
        $runRes = Invoke-Captured -Exe $py -Arguments $runArgs
        $out = $runRes.Lines
        $out | Out-File (Join-Path $logDir "$($Tag)_$($p.Key)_$($arm.ArmTag)_stdout.txt") -Encoding utf8
        $out | Select-String -Pattern 'Run state|DEGENERATE RUN|^\s*(ridge|panoptes)\s' |
            ForEach-Object { Write-Host "      $_"; Add-Content $log "      $_" }

        # Rule 3 of day1 section 11: read the run directory from the run you started, by
        # name, never from "the newest directory".
        $hit = $out | Select-String -Pattern 'run_\d{8}_\d{6}_[0-9a-f]+' | Select-Object -First 1
        $dir = if ($hit) { $hit.Matches[0].Value } else { '(none reported)' }
        $done = "## $($p.Key) [$($arm.ArmTag)] done $(Get-Date -Format HH:mm:ss) exit=$($runRes.ExitCode) run_dir=$dir"
        Write-Host "   $done"; Add-Content $log $done

        if ($runRes.ExitCode -ne 0) { $failed += "$($p.Key):$($arm.ArmTag):run" }

        # Post-run gate WITH --run-dir. This is the one that actually exercises
        # `threshold` and `budget`.
        if ($dir -ne '(none reported)') {
            $rd = Join-Path $logDir $dir
            if ((Invoke-Gate -Key $p.Key -Budget $budget -Initial $p.Initial -RunDir $rd -Label "post-run $($arm.ArmTag)") -ne 0) {
                $failed += "$($p.Key):$($arm.ArmTag):gate"
            }
        }
    }
}

Write-Host ''
Write-Host "################ campaign complete $(Get-Date -Format HH:mm:ss) ################"
Add-Content $log "campaign complete $(Get-Date -Format 'o')"
Write-Host "log: $log"
if ($failed.Count -gt 0) {
    # Collected rather than aborting: a failed gate means that landscape must not be
    # quoted (rule 2.7), not that the rest of the campaign is worthless.
    Write-Host "NON-ZERO exits: $($failed -join ', ')" -ForegroundColor Red
    Add-Content $log "NON-ZERO: $($failed -join ', ')"
    exit 1
}
exit 0

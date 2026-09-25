param([int]$RegenerationPid = 0)
$ErrorActionPreference = 'Stop'
$rscript = 'C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe'
if ($RegenerationPid -gt 0) {
    Wait-Process -Id $RegenerationPid -ErrorAction SilentlyContinue
}
$pipelineLog = Join-Path $env:TEMP 'issue181_regenerate.out.log'
if (!(Select-String -LiteralPath $pipelineLog -Pattern 'ISSUE181 COMPLETE 12_scenario_EVPPIs' -Quiet)) {
    throw 'Regeneration did not complete; inspect issue181_regenerate logs before continuing.'
}
function Run-Step([string]$Name, [string]$Executable, [string[]]$Arguments) {
    Write-Output ("START {0} {1}" -f $Name, (Get-Date -Format o))
    $logPath = Join-Path $env:TEMP ("issue181-{0}.log" -f $Name)
    $ErrorActionPreference = 'Continue' # R/Quarto write progress and warnings to stderr.
    & $Executable @Arguments *> $logPath
    $code = $LASTEXITCODE
    $ErrorActionPreference = 'Stop'
    if ($code -ne 0) {
        Get-Content -LiteralPath $logPath -Tail 40
        throw ("{0} failed with exit code {1}" -f $Name, $code)
    }
    Write-Output ("COMPLETE {0} {1}" -f $Name, (Get-Date -Format o))
}
Run-Step 'fixed-snapshot' $rscript @('analysis/13_save_snapshot.R', '181', 'fixed')
Run-Step 'snapshot-comparison' $rscript @('tests/compare_snapshots.R', '181')
foreach ($report in @('pfs_os_violations', 'psa_extrapolation_plausibility')) {
    foreach ($format in @('pdf', 'gfm')) {
        Run-Step "$report-$format" 'quarto' @('render', "reports/technical/$report.qmd", '--to', $format)
    }
}
Run-Step 'all-reports' $rscript @('publish/publish_reports.R')
Run-Step 'all-vignettes' 'powershell' @('-NoProfile', '-File', 'validation/issue168_2026-09-22/render_artifacts.ps1')
Run-Step 'report-contracts' $rscript @('tests/test_report_contracts.R')
Write-Output 'All issue181 regeneration, snapshot, render and report-contract steps completed.'

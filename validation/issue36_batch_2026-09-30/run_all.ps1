# Issue #36/#164/#56 batch: regenerate caches, render every report and
# publication vignette, check report contracts, then save the fixed snapshot,
# compare it with the baseline and re-render the bug-fix impact report.
# Run from the repository root after committing the source changes.
$ErrorActionPreference = 'Stop'
$rscript = 'C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe'
$logDir = Join-Path $env:TEMP 'issue36_batch'
New-Item -ItemType Directory -Force $logDir | Out-Null
function Run-Step([string]$Name, [string]$Executable, [string[]]$Arguments) {
    Write-Output ("START {0} {1}" -f $Name, (Get-Date -Format o))
    $logPath = Join-Path $logDir ("{0}.log" -f $Name)
    $ErrorActionPreference = 'Continue' # R/Quarto write progress and warnings to stderr.
    & $Executable @Arguments *> $logPath
    $code = $LASTEXITCODE
    $ErrorActionPreference = 'Stop'
    if ($code -ne 0) {
        Get-Content -LiteralPath $logPath -Tail 40
        throw ("FAILED {0} with exit code {1}; log {2}" -f $Name, $code, $logPath)
    }
    Write-Output ("COMPLETE {0} {1}" -f $Name, (Get-Date -Format o))
}
Run-Step 'regenerate' $rscript @('validation/issue36_batch_2026-09-30/regenerate.R')
Run-Step 'all-reports' $rscript @('publish/publish_reports.R')
Run-Step 'all-vignettes' 'powershell' @('-NoProfile', '-File', 'validation/issue168_2026-09-22/render_artifacts.ps1')
Run-Step 'report-contracts' $rscript @('tests/test_report_contracts.R')
Run-Step 'fixed-snapshot' $rscript @('analysis/13_save_snapshot.R', '36', 'fixed')
Run-Step 'snapshot-comparison' $rscript @('tests/compare_snapshots.R', '36')
foreach ($format in @('pdf', 'gfm')) {
    Run-Step "bug_fix_impact-$format" 'quarto' @('render', 'reports/technical/bug_fix_impact.qmd', '--to', $format)
}
Write-Output 'ALL DONE issue36 batch.'

param([int]$ReportsPid = 0)
$ErrorActionPreference = 'Stop'
if ($ReportsPid -gt 0) { Wait-Process -Id $ReportsPid -ErrorAction SilentlyContinue }
$reportsLog = Join-Path $env:TEMP 'issue181_reports_retry.out.log'
if (!(Select-String -LiteralPath $reportsLog -Pattern '_site/ assembled:' -Quiet)) {
    throw 'Report renders did not complete; inspect issue181_reports_retry logs.'
}
$ErrorActionPreference = 'Continue'
& powershell -NoProfile -File validation/issue168_2026-09-22/render_artifacts.ps1 *> "$env:TEMP\issue181-vignettes-retry.log"
if ($LASTEXITCODE -ne 0) { throw 'Publication vignette rendering failed.' }
& 'C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe' tests/test_report_contracts.R *> "$env:TEMP\issue181-report-contracts.log"
if ($LASTEXITCODE -ne 0) { throw 'Report contracts failed.' }
Write-Output 'All reports, publication vignettes and report contracts completed.'

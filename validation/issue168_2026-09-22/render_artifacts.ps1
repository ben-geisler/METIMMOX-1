# Run from the repository root. One writer at a time owns outputs/manifest.csv.
$ErrorActionPreference = 'Stop'
$vignettes = Get-ChildItem -LiteralPath 'outputs/vignettes' -Filter '*.qmd' | Sort-Object Name
foreach ($vignette in $vignettes) {
    $logPath = Join-Path $env:TEMP ('issue168-render-' + $vignette.BaseName + '.log')
    Write-Output ('Rendering ' + $vignette.Name)
    $ErrorActionPreference = 'Continue' # Quarto progress uses stderr.
    & quarto render $vignette.FullName *> $logPath
    $ErrorActionPreference = 'Stop'
    if ($LASTEXITCODE -ne 0) {
        Get-Content -LiteralPath $logPath -Tail 35
        throw ('Render failed: ' + $vignette.Name)
    }
}
Write-Output 'All publication vignettes rendered successfully.'

param(
    [string]$Python = "python",
    [string]$OutputDir = "dist"
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..\")).Path
Set-Location $repo

& $Python -m pip show pyinstaller *> $null
if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller is required. Install it with: $Python -m pip install pyinstaller"
}

& $Python -m PyInstaller --noconfirm --clean --windowed `
    --name AlFatigue `
    --distpath $OutputDir `
    --workpath (Join-Path $OutputDir "build") `
    --specpath (Join-Path $OutputDir "spec") `
    simulations\fem_tension_app.py

if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller failed."
}

Write-Host "Built $OutputDir\AlFatigue\AlFatigue.exe"

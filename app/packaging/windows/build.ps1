param(
    [string]$Python = "python",
    [string]$OutputDir = "dist"
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..\")).Path
Set-Location $repo

$solver = Join-Path $repo "fem1d\bin\fem1d_solver.exe"
if (-not (Test-Path $solver)) {
    throw "Build fem1d\bin\fem1d_solver.exe before packaging the desktop app."
}

& $Python -m pip show pyinstaller *> $null
if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller is required. Install it with: $Python -m pip install pyinstaller"
}

& $Python -m PyInstaller --noconfirm --clean --windowed `
    --name AlFatigue `
    --distpath $OutputDir `
    --workpath (Join-Path $OutputDir "build") `
    --specpath (Join-Path $OutputDir "spec") `
    --add-binary "$solver;fem1d\bin" `
    simulations\fem_tension_app.py

if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller failed."
}

Write-Host "Built $OutputDir\AlFatigue\AlFatigue.exe"

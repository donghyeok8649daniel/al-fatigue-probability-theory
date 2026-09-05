param(
    [string]$Python = "python",
    [string]$OutputDir = "dist"
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..\")).Path
Set-Location $repo

$solver = Join-Path $repo "fem1d\bin\fem1d_solver.exe"
$cylinder = Join-Path $repo "examples\meshes\default_tensile_cylinder.stl"
$solverArgs = @()
if (Test-Path $solver) {
    $solverArgs = @("--add-binary", "$solver;fem1d\bin")
} else {
    Write-Host "FEM binary not found; packaging the bundled Python FEM fallback."
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
    --add-data "$cylinder;examples\meshes" `
    $solverArgs `
    app\desktop_ui.py

if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller failed."
}

Write-Host "Built $OutputDir\AlFatigue\AlFatigue.exe"

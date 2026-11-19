@echo off
setlocal
set "OPENBLAS_NUM_THREADS=1"
set "OMP_NUM_THREADS=1"
set "MKL_NUM_THREADS=1"
if defined CELL2D_PYTHON goto runtime_ready
set "CELL2D_PYTHON=%LOCALAPPDATA%\Programs\Python\Python313\python.exe"
if exist "%CELL2D_PYTHON%" goto runtime_ready
set "CELL2D_PYTHON=%USERPROFILE%\AppData\Local\Programs\Python\Python313\python.exe"
if exist "%CELL2D_PYTHON%" goto runtime_ready
set "CELL2D_PYTHON=python"
:runtime_ready
pushd "%~dp0"
"%CELL2D_PYTHON%" -X utf8 independent_cells/run_family.py --out outputs/fvm-cell2d/new_run --base-s 1 --refine 3 --workers 4 %*
set "CELL2D_EXIT=%ERRORLEVEL%"
popd
exit /b %CELL2D_EXIT%

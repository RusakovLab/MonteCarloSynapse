@echo off
setlocal
REM ============================================================
REM  Launch the glutamate diffusion / receptor-kinetics GUI.
REM  Unzip the repository, keep this file next to gui.py and
REM  double-click it. No machine-specific paths: the Python
REM  interpreter is discovered and missing packages are offered
REM  for installation.
REM ============================================================

cd /d "%~dp0"
echo Starting the Monte Carlo synapse GUI...
echo Folder: %cd%
echo.

set "PYEXE="

REM 1) An activated environment, or a local .venv next to this file.
if defined VIRTUAL_ENV if exist "%VIRTUAL_ENV%\Scripts\python.exe" set "PYEXE=%VIRTUAL_ENV%\Scripts\python.exe"
if not defined PYEXE if exist "%~dp0.venv\Scripts\python.exe" set "PYEXE=%~dp0.venv\Scripts\python.exe"

REM 2) Whatever is on PATH (python, then the py launcher).
if not defined PYEXE for /f "delims=" %%I in ('where python 2^>nul') do if not defined PYEXE set "PYEXE=%%I"
if not defined PYEXE where py >nul 2>nul && set "PYEXE=py"

REM 3) Common per-user Anaconda / Miniconda / python.org installs.
if not defined PYEXE for %%D in (
    "%USERPROFILE%\anaconda3\python.exe"
    "%USERPROFILE%\miniconda3\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python313\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python311\python.exe"
) do if not defined PYEXE if exist "%%~D" set "PYEXE=%%~D"

if not defined PYEXE (
    echo *** No Python interpreter found. ***
    echo Install Python 3.9 or newer from https://www.python.org/downloads/
    echo and tick "Add python.exe to PATH" during setup, then run this file again.
    echo.
    pause
    exit /b 1
)

echo Using Python: %PYEXE%

REM --- check the three required packages, offer to install them ---
"%PYEXE%" -c "import numpy, scipy, matplotlib" >nul 2>nul
if %errorlevel% equ 0 goto run

echo.
echo Some required packages are missing (numpy / scipy / matplotlib).
set "INSTALL=Y"
set /p INSTALL="Install them now with pip? [Y/n] "
if /i "%INSTALL%"=="n" goto abort
"%PYEXE%" -m pip install -r requirements.txt
if %errorlevel% neq 0 goto pipfail

:run
echo.
"%PYEXE%" gui.py
if %errorlevel% neq 0 (
    echo.
    echo *** The program exited with an error ^(code %errorlevel%^). ***
    echo.
    pause
)
endlocal
exit /b 0

:abort
echo Cannot start without them. Exiting.
pause
endlocal
exit /b 1

:pipfail
echo.
echo *** Installation failed. Try:
echo     "%PYEXE%" -m pip install --user -r requirements.txt
echo.
pause
endlocal
exit /b 1

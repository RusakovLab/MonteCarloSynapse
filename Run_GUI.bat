@echo off
setlocal enabledelayedexpansion
REM ============================================================
REM  Launch the glutamate diffusion / receptor-kinetics GUI.
REM  Unzip the repository, keep this file next to gui.py and
REM  double-click it. Nothing here is specific to one machine:
REM  the script looks through the Python interpreters on this
REM  computer and picks one that can run the program.
REM ============================================================

cd /d "%~dp0"
echo Starting the Monte Carlo synapse GUI...
echo Folder: %cd%
echo.

set "PYEXE="
set "FALLBACK="

REM --- Collect candidates, best first, into a temporary list. ---
set "CANDLIST=%TEMP%\mcs_pythons_%RANDOM%.txt"
break > "%CANDLIST%"

if defined VIRTUAL_ENV if exist "%VIRTUAL_ENV%\Scripts\python.exe" echo %VIRTUAL_ENV%\Scripts\python.exe>> "%CANDLIST%"
if exist "%~dp0.venv\Scripts\python.exe" echo %~dp0.venv\Scripts\python.exe>> "%CANDLIST%"

REM Conda and python.org installs, which usually carry the packages.
for %%D in (
    "%USERPROFILE%\anaconda3\python.exe"
    "%USERPROFILE%\miniconda3\python.exe"
    "%USERPROFILE%\Anaconda3\python.exe"
    "%USERPROFILE%\miniforge3\python.exe"
    "%PROGRAMDATA%\anaconda3\python.exe"
    "%PROGRAMDATA%\miniconda3\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python313\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python311\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python310\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python39\python.exe"
) do if exist "%%~D" echo %%~D>> "%CANDLIST%"

REM Whatever is on PATH, and the py launcher. These come last: on Windows
REM the first python.exe on PATH is often the Microsoft Store stub, which
REM has no scientific packages installed.
for /f "delims=" %%I in ('where python 2^>nul') do echo %%I>> "%CANDLIST%"
for /f "delims=" %%I in ('where py 2^>nul') do echo %%I>> "%CANDLIST%"

REM --- Pick the first candidate that already has the packages. ---
for /f "usebackq delims=" %%P in ("%CANDLIST%") do (
    if not defined PYEXE (
        "%%P" -c "import sys; sys.exit(0)" >nul 2>nul
        if !errorlevel! equ 0 (
            if not defined FALLBACK set "FALLBACK=%%P"
            "%%P" -c "import numpy, scipy, matplotlib" >nul 2>nul
            if !errorlevel! equ 0 set "PYEXE=%%P"
        )
    )
)
del "%CANDLIST%" >nul 2>nul

if defined PYEXE goto run

REM --- Nothing ready to go: offer to install into the best interpreter. ---
if not defined FALLBACK (
    echo *** No Python interpreter found. ***
    echo Install Python 3.9 or newer from https://www.python.org/downloads/
    echo and tick "Add python.exe to PATH" during setup, then run this file again.
    echo.
    pause
    exit /b 1
)

set "PYEXE=%FALLBACK%"
echo Using Python: %PYEXE%
echo.
echo Some required packages are missing (numpy / scipy / matplotlib).
set "INSTALL=Y"
set /p INSTALL="Install them now with pip? [Y/n] "
if /i "%INSTALL%"=="n" goto abort
"%PYEXE%" -m pip install -r requirements.txt
if %errorlevel% neq 0 goto pipfail
goto start

:run
echo Using Python: %PYEXE%

:start
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

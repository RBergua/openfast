@echo off
:: TEMPORARY bisection harness for OLAF_nGridOut (PR #3414) Intel Fortran
:: console-output investigation. NOT meant to be merged upstream - remove
:: this script and the corresponding CI step once the root cause is fixed.
::
:: Runs AeroDyn_Driver.exe against a repro case whose OLAF input deliberately
:: sets DTfvw < DTaero, which triggers the following fatal, unconditional
:: validation error in FVW_IO.f90:
::   "DTfvw must be >= DTaero from AD15."
:: The driver is expected to abort in both "good" and "bad" states, so we do
:: NOT gate on the process exit code - only on whether the message was
:: actually printed to the console.

setlocal enabledelayedexpansion

set "CASE_DIR=reg_tests\manualTestCases\OLAF_bug_repro"
set "LOG_FILE=%CASE_DIR%\olaf_bisect_harness.log"
set "MARKER=DTfvw must be >= DTaero"

echo on
build\bin\AeroDyn_Driver.exe "%CASE_DIR%\OLAF_LC11.dvr" > "%LOG_FILE%" 2>&1
echo off

echo.
echo ----- Captured console output -----
type "%LOG_FILE%"
echo ------------------------------------

findstr /C:"%MARKER%" "%LOG_FILE%" >nul
if %ERRORLEVEL% EQU 0 (
    echo PASS: expected message found in console output - "%MARKER%"
    exit /b 0
) else (
    echo FAIL: expected message NOT found in console output - "%MARKER%"
    echo This reproduces the Intel Fortran missing-console-output bug.
    exit /b 1
)

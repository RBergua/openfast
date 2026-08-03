@echo off
:: Runs the OLAF test case and saves the console output to a log file.

set "CASE_DIR=reg_tests\manualTestCase"
set "LOG_FILE=%CASE_DIR%\console.log"

build\bin\AeroDyn_Driver.exe "%CASE_DIR%\OLAF_LC11.dvr" > "%LOG_FILE%" 2>&1

echo.
echo ----- Captured console output -----
type "%LOG_FILE%"
echo ------------------------------------
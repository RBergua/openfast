@call "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" %VS_VER%

for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

@REM Make the script that generates the git version description ignore dirty
@REM since building the Visual Studio projects modifies files
powershell -command "(Get-Content -Path '.\vs-build\CreateGitVersion.bat') -replace '--dirty', '' | Set-Content -Path '.\vs-build\CreateGitVersion.bat'"

echo on


setlocal enabledelayedexpansion

:: Initialize a variable to store failed solutions
set "FailedSolutions="
set "OverallErrorLevel=0"


@REM Build all solutions (release 64)
devenv vs-build/OpenFAST.sln /Build "Release|x64"
if %ERRORLEVEL% NEQ 0 (
    set "FailedSolutions=!FailedSolutions!OpenFAST Release "
    set "OverallErrorLevel=1"
    echo Build of OpenFAST.sln Release failed!
)


@REM Build all OpenMP solutions (release 64 OpenMP)
devenv vs-build/OpenFAST.sln /Build "Release_OpenMP|x64"
if %ERRORLEVEL% NEQ 0 (
    set "FailedSolutions=!FailedSolutions!OpenFAST Release_OpenMP "
    set "OverallErrorLevel=1"
    echo Build of OpenFAST.sln Release_OpenMP failed!
)


@REM Build MATLAB solution last
devenv vs-build/OpenFAST.sln /Build "Release_Matlab|x64"
if %ERRORLEVEL% NEQ 0 (
    set "FailedSolutions=!FailedSolutions!OpenFAST Release_Matlab "
    set "OverallErrorLevel=1"
    echo Build of OpenFAST.sln Release_Matlab failed!
)



@REM Copy controllers to bin directory
@REM xcopy .\reg_tests\r-test\glue-codes\openfast\5MW_Baseline\ServoData\*.dll .\build\bin\ /y

echo.
echo Build Summary:
if defined FailedSolutions (
    echo The following solutions failed to build:
    echo %FailedSolutions%
) else (
    echo All solutions built successfully.
)

:: Set the final error level based on the overall build status
exit /b %OverallErrorLevel%

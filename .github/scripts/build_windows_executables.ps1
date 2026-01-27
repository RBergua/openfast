# Stop on unhandled errors
$ErrorActionPreference = 'Stop'

# --- Set up Intel oneAPI + Visual Studio environment ---

# VS version must already be set in the environment, same as %VS_VER%
& "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" $env:VS_VER

# Find the latest Intel compiler version (excluding "latest")
$compilerRoot = "C:\Program Files (x86)\Intel\oneAPI\compiler"

$LatestVersion =
    Get-ChildItem -Path $compilerRoot -Directory |
    Where-Object { $_.Name -ne 'latest' } |
    Sort-Object Name |
    Select-Object -Last 1 -ExpandProperty Name

& "$compilerRoot\$LatestVersion\env\vars.bat"


# --- Modify CreateGitVersion.bat to ignore --dirty ---

# Visual Studio builds modify files, so remove "--dirty"
$gitVersionBat = ".\vs-build\CreateGitVersion.bat"

(Get-Content $gitVersionBat) `
    -replace '--dirty', '' |
    Set-Content $gitVersionBat


# --- Build tracking variables ---

$FailedSolutions = @()
$OverallErrorLevel = 0


# Helper function to build a solution/config
function Invoke-Build {
    param (
        [string]$ConfigName,
        [string]$BuildSpec
    )

    Write-Host "Build $ConfigName"
    & devenv vs-build/OpenFAST.sln /Build $BuildSpec

    if ($LASTEXITCODE -ne 0) {
        $script:FailedSolutions += $ConfigName
        $script:OverallErrorLevel = 1
        Write-Host "Build of OpenFAST.sln $ConfigName failed!" -ForegroundColor Red
    }
}


# --- Build steps ---

Invoke-Build "Release"         "Release|x64"
Invoke-Build "OpenMP_Release"  "OpenMP_Release|x64"
Invoke-Build "Matlab_Release"  "Matlab_Release|x64"


# --- Build summary ---

Write-Host "Build Summary:"
if ($FailedSolutions.Count -gt 0) {
    Write-Host "The following solutions failed to build:"
    $FailedSolutions -join '  ' | Write-Host
} else {
    Write-Host "All solutions built successfully."
}


# --- Rename output files ---

Write-Host "Remove '_Release' and '_Matlab' from file names"

Push-Location .\build\bin
Get-ChildItem -File -Filter '*_Release*' |
    Rename-Item -NewName { $_.Name -replace '_Release', '' }

Get-ChildItem -File -Filter '*_Matlab*' |
    Rename-Item -NewName { $_.Name -replace '_Matlab', '' }
Pop-Location


# --- List executables ---

Write-Host "List executables in build\bin"
Get-ChildItem build\bin


# --- Exit with overall status ---

exit $OverallErrorLevel

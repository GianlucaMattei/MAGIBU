# Builds the Magibu Windows installer. Run it on a Windows PC, from the App folder:
#   powershell -ExecutionPolicy Bypass -File build\build_windows.ps1
# Needs: internet connection, Node.js (https://nodejs.org). R is downloaded here, no need to have it installed.
$ErrorActionPreference = "Stop"
$RVersion = "4.4.2"
$root  = Split-Path -Parent $PSScriptRoot
$build = $PSScriptRoot
$rdir  = Join-Path $build "R-portable"

if (-not (Test-Path (Join-Path $rdir "bin\Rscript.exe"))) {
  $inst = Join-Path $build "R-$RVersion-win.exe"
  if (-not (Test-Path $inst)) {
    Write-Host "Downloading R $RVersion..."
    Invoke-WebRequest "https://cran.r-project.org/bin/windows/base/old/$RVersion/R-$RVersion-win.exe" -OutFile $inst
  }
  Write-Host "Installing R in $rdir (no administrator rights needed)..."
  Start-Process $inst -ArgumentList "/VERYSILENT","/CURRENTUSER","/NOICONS","/DIR=`"$rdir`"","/MERGETASKS=!desktopicon,!recordversion,!associate" -Wait
}
$rscript = Join-Path $rdir "bin\Rscript.exe"

Write-Host "Installing the R packages (it can take 20-40 minutes)..."
& $rscript --vanilla (Join-Path $build "install_packages.R")
if ($LASTEXITCODE -ne 0) { throw "Package installation failed" }

# slimmer: drop package documentation and tests
Get-ChildItem (Join-Path $rdir "library") -Directory | ForEach-Object {
  foreach ($d in "doc","help","html","tests","demo") { $p = Join-Path $_.FullName $d; if (Test-Path $p) { Remove-Item $p -Recurse -Force } }
}
foreach ($d in "doc","tests","Tcl") { $p = Join-Path $rdir $d; if (Test-Path $p) { Remove-Item $p -Recurse -Force } }

& $rscript --vanilla (Join-Path $build "make_licenses.R") (Join-Path $rdir "library") (Join-Path $build "THIRD_PARTY_LICENSES.txt")

Write-Host "Building the installer..."
Set-Location (Join-Path $root "electron")
npm install --no-audit --no-fund
npm run dist
Write-Host "Done. Installer in: $(Join-Path $root 'dist')"

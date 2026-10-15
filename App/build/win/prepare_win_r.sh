#!/usr/bin/env bash
# Prepares, on Linux, the R for Windows bundled in the zip (build/win/R-extract/app), used by
# electron/builder-win-cross.json. Needs: curl, R with the data packages listed in zips/data_pkgs.txt installed
# (Bioconductor annotation/manifest packages, platform independent, copied as they are).
# Then: cd ../../electron && npm install && npx electron-builder --win zip --x64 --config builder-win-cross.json
set -e
cd "$(dirname "$0")"
RVER=4.4.2
[ -x innoextract-1.9-linux/innoextract ] || {
  curl -sL -o innoextract.tar.xz https://github.com/dscharrer/innoextract/releases/download/1.9/innoextract-1.9-linux.tar.xz
  tar xf innoextract.tar.xz; }
[ -f R-$RVER-win.exe ] || curl -sL -o R-$RVER-win.exe https://cran.r-project.org/bin/windows/base/old/$RVER/R-$RVER-win.exe
[ -d R-extract/app ] || ./innoextract-1.9-linux/innoextract -s -d R-extract R-$RVER-win.exe
Rscript fetch_win_pkgs.R
WL=R-extract/app/library
for p in $(cat zips/data_pkgs.txt); do
  d=$(Rscript -e "cat(find.package('$p'))")
  [ -d "$WL/$p" ] || cp -r "$d" "$WL/"
done
# slimmer: BH is only needed at compile time; drop documentation and tests
rm -rf "$WL/BH" R-extract/app/doc R-extract/app/tests R-extract/app/Tcl
for d in doc help html tests demo; do find "$WL" -mindepth 2 -maxdepth 2 -type d -name $d -exec rm -rf {} +; done
Rscript ../make_licenses.R "$PWD/$WL" ../THIRD_PARTY_LICENSES.txt
du -sh R-extract/app

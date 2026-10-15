#!/usr/bin/env bash
# Starts Magibu on a Linux server (no executable) at http://127.0.0.1:<port>.
# Usage: ./run_local.sh [port]    then from your PC: ssh -L 7654:127.0.0.1:7654 <server> and open http://127.0.0.1:7654
HERE="$(cd "$(dirname "$0")" && pwd)"
PORT=${1:-7654}
# rlib-server (optional): extra R library searched first, e.g. irlba rebuilt against the installed Matrix
export R_LIBS="$HERE/rlib-server${R_LIBS:+:$R_LIBS}"
export MAGIBU_APP_DIR="$HERE/app" MAGIBU_PORT="$PORT" MAGIBU_MODELS="$HERE/models"
exec Rscript "$HERE/app/launch.R"

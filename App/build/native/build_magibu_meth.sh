#!/bin/bash
# Rebuilds magibu-meth (BAM methylation reader) for Linux and Windows and copies it into app/bin.
# Needs the static htslib already built in htslib-linux/ (gcc) and htslib-win/ (zig, mingw target),
# plus win/{lib,include} (zlib, winpthreads, regex.h stub). Sources in src/, downloads in dl/.
set -e
cd "$(dirname "$0")"
N=$PWD
Z=$N/src/zig-linux-x86_64-0.13.0/zig
W=$N/win
mkdir -p out ../../app/bin/linux-x64 ../../app/bin/win-x64
gcc -O2 -Wall -Ihtslib-linux -o out/magibu-meth_linux prog/magibu_meth.c htslib-linux/libhts.a -lz -lm -lpthread
$Z cc -target x86_64-windows-gnu -O2 -s -D_FILE_OFFSET_BITS=64 -DWINPTHREAD_STATIC -Ihtslib-win -I$W/include \
  -o out/magibu-meth.exe prog/magibu_meth.c htslib-win/libhts.a -L$W/lib -lz -lpthread -lws2_32
cp out/magibu-meth_linux ../../app/bin/linux-x64/magibu-meth && strip ../../app/bin/linux-x64/magibu-meth
cp out/magibu-meth.exe ../../app/bin/win-x64/magibu-meth.exe
ls -la ../../app/bin/*/

#!/usr/bin/env bash
# Builds CAPD 6.1.0 (the commit used for the paper) into $1 (default: ./capd-install).
set -euo pipefail
PREFIX=$(realpath -m "${1:-capd-install}")
SRC=$(mktemp -d)
git clone --quiet https://github.com/CAPDGroup/CAPD "$SRC/CAPD"
git -C "$SRC/CAPD" checkout --quiet 2f060980d0685cc9bf88352e08c242197cd7d686
cmake -S "$SRC/CAPD" -B "$SRC/build" -DCMAKE_BUILD_TYPE=Release -DCAPD_ENABLE_MULTIPRECISION=OFF \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" > "$SRC/cmake.log"
cmake --build "$SRC/build" -j"$(nproc)" > "$SRC/build.log"
cmake --install "$SRC/build" > "$SRC/install.log"
echo "CAPD installed in $PREFIX; capd-config: $PREFIX/bin/capd-config"

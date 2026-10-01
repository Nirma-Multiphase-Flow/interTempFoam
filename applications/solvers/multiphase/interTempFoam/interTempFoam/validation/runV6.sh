#!/bin/sh
# Usage: validation/runV6.sh <tag> <R0/h> [nProcs] [tEnd]
set -e
here=$(cd "$(dirname "$0")" && pwd)
dst=${FOAM_RUN:?}/itf/$1/V6_$2
rm -rf "$dst"; mkdir -p "$(dirname "$dst")"
cp -r "$here/V6_d2law_axi" "$dst"
cd "$dst" && ./Allrun "$2" "${3:-2}" "${4:-2e-3}"
cat "$dst/V6_metrics.txt"

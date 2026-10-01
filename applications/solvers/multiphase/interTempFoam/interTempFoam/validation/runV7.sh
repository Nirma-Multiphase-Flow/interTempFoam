#!/bin/sh
# Usage: validation/runV7.sh <tag> <h_um> [nProcs]
set -e
here=$(cd "$(dirname "$0")" && pwd)
dst=${FOAM_RUN:?}/itf/$1/V7_$2
rm -rf "$dst"; mkdir -p "$(dirname "$dst")"
cp -r "$here/V7_bubble_rise_axi" "$dst"
cd "$dst" && ./Allrun "$2" "${3:-4}"
cat "$dst/V7_metrics.txt"

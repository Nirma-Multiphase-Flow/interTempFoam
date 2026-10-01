#!/bin/sh
# Usage: validation/runV5.sh <tag> <R0/h> [nProcs]
set -e
here=$(cd "$(dirname "$0")" && pwd)
dst=${FOAM_RUN:?}/itf/$1/V5_$2
rm -rf "$dst"; mkdir -p "$(dirname "$dst")"
cp -r "$here/V5_scriven_axi" "$dst"
cd "$dst" && ./Allrun "$2" "${3:-4}"
cat "$dst/V5_metrics.txt"

#!/bin/sh
# Usage: validation/runV8.sh <tag> <nx> [nProcs]
set -e
here=$(cd "$(dirname "$0")" && pwd)
dst=${FOAM_RUN:?}/itf/$1/V8_$2
rm -rf "$dst"; mkdir -p "$(dirname "$dst")"
cp -r "$here/V8_film_boiling" "$dst"
cd "$dst" && ./Allrun "$2" "${3:-4}"
cat "$dst/V8_metrics.txt"

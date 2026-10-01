#!/bin/sh
# Usage: validation/runV34.sh <tag> <stefan|sucking> <N>
set -e
here=$(cd "$(dirname "$0")" && pwd)
dst=${FOAM_RUN:?}/itf/$1/V34_$2_$3
rm -rf "$dst"; mkdir -p "$(dirname "$dst")"
cp -r "$here/V34_stefan1d" "$dst"
cd "$dst" && ./Allrun "$2" "$3"
cat "$dst/V34_metrics.txt"

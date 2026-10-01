#!/bin/sh
# Usage: validation/runTranslate.sh <tag> <R/h>
set -e
here=$(cd "$(dirname "$0")" && pwd)
dst=${FOAM_RUN:?}/itf/$1/translate_$2
rm -rf "$dst"; mkdir -p "$(dirname "$dst")"
cp -r "$here/V1_static_droplet" "$dst"
cd "$dst" && ./Allrun.translate "$2"
cat "$dst/translate_metrics.txt"

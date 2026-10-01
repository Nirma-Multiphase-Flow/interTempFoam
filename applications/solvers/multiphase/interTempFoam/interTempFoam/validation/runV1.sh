#!/bin/sh
# Usage: validation/runV1.sh <tag> <planar|wedge> <R/h> [nSteps] [curvatureSettingsFile]
# Runs one V1 instance in $FOAM_RUN/itf/<tag>/V1_<geom>_<R/h>.
set -e
here=$(cd "$(dirname "$0")" && pwd)
dst=${FOAM_RUN:?}/itf/$1/V1_$2_$3
rm -rf "$dst"; mkdir -p "$(dirname "$dst")"
cp -r "$here/V1_static_droplet" "$dst"
[ -n "$5" ] && cp "$5" "$dst/constant/curvatureSettings"
cd "$dst" && ./Allrun "$2" "$3" "$4"
cat "$dst/V1_metrics.txt"

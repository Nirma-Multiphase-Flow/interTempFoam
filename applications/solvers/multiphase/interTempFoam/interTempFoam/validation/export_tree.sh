#!/bin/sh
# Export interTempFoam (phase-change models are part of the solver directory)
# and the Sun et al. (2014) film-boiling case into an OpenFOAM user tree laid
# out as
#   <root>/applications/solvers/interTempFoam  (incl. thermalPhaseChangeModels/)
#   <root>/tutorials/filmBoiling_Sun2014
# Usage: export_tree.sh <root>
set -e
root=${1:?root}
here=$(cd "$(dirname "$0")/.." && pwd)                 # solver dir

mkdir -p "$root/applications/solvers" "$root/tutorials"

# solver: sources and Make/{files,options} only
rm -rf "$root/applications/solvers/interTempFoam"
mkdir -p "$root/applications/solvers/interTempFoam/Make"
cp "$here"/*.C "$here"/*.H "$root/applications/solvers/interTempFoam/"
mkdir -p "$root/applications/solvers/interTempFoam/thermalPhaseChangeModels"
cp "$here"/thermalPhaseChangeModels/*.C "$here"/thermalPhaseChangeModels/*.H "$root/applications/solvers/interTempFoam/thermalPhaseChangeModels/"
cp "$here/Make/files" "$here/Make/options" "$root/applications/solvers/interTempFoam/Make/"
cp "$here/README.md" "$root/applications/solvers/interTempFoam/" 2>/dev/null || true

# case
rm -rf "$root/tutorials/filmBoiling_Sun2014"
python3 "$here/validation/V9_filmBoiling_Sun2014/make_case.py" "$root/tutorials/filmBoiling_Sun2014" 100 3.0 50 8
echo "exported to $root"

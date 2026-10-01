#!/bin/sh
# Usage: [SIGMA=<value>] [RUNNER=<wrapper>] validation/runCase.sh <caseDir> <tag>
# Copies the case template to $FOAM_RUN/itf/<tag>/<case> and runs its Allrun.
# SIGMA overrides the surface tension coefficient (diagnostic runs only);
# RUNNER prefixes Allrun, e.g. the baseline-binary wrapper.
# Templates stay clean in git; results live outside the repo tree.
set -e
src=$(cd "$(dirname "$0")/$1" && pwd)
dst=${FOAM_RUN:?}/itf/$2/$(basename "$1")
rm -rf "$dst"
mkdir -p "$(dirname "$dst")"
cp -r "$src" "$dst"
cd "$dst"
if [ -n "$SIGMA" ]; then
    sed -i "s/^sigma .*/sigma $SIGMA;   \/\/ overridden by runCase.sh SIGMA/" constant/transportProperties
    grep "^sigma" constant/transportProperties
fi
$RUNNER ./Allrun
echo "done: $dst"

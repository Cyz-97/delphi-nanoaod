#!/bin/bash

# Print the arguments
echo "Arguments passed to this script:"
echo "  input file  : $1"
echo "  output file : $2"
echo "  move to     : $3"
echo "  data or mc  : $4"

# Create a temporary working directory
WORKDIR=$(mktemp -d /tmp/delphi_job_XXXXXX)
echo "Working directory: $WORKDIR"
trap 'rm -rf "$WORKDIR"' EXIT

# Copy everything into the temp directory
BASEDIR="/afs/cern.ch/user/z/zhangj/private/DELPHI/delphi-nanoaod"
cp "$BASEDIR/build/delphi-nanoaod/delphi-nanoaod" "$WORKDIR/"
cp "$BASEDIR/config/delphi-nanoaod.yaml" "$WORKDIR/"
cp "$BASEDIR/scripts/treefy.C" "$WORKDIR/"
cd "$WORKDIR" || exit 1

# Set up environments
source /cvmfs/delphi.cern.ch/setup.sh
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.34.04/x86_64-almalinux9.5-gcc115-opt/bin/thisroot.sh

# Write unique PDLINPUT file
echo "FILE=$1" > "${2}_dummy"

# Run nanoAOD producer
if [ "$4" = "MC" ]; then
    delphi-nanoaod -P "${2}_dummy" --mc --config delphi-nanoaod.yaml --output "$2"
else
    delphi-nanoaod -P "${2}_dummy" --config delphi-nanoaod.yaml --output "$2"
fi

# Run treefy step
#root -q -b -l "treefy.C+(\"$2\")"

# Move the output ROOT file to desired location
tpc="$2"
#nanotree="${tpc/.root/_ttree.root}"

mv "$tpc" "$3"
#mv "$nanotree" "$3"



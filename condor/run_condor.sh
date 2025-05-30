#!/bin/bash

# Move to the working directory
cd /afs/cern.ch/user/z/zhangj/private/DELPHI/delphi-nanoaod || exit 1

# Set up DELPHI and ROOT environments
source /cvmfs/delphi.cern.ch/setup.sh
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.34.04/x86_64-almalinux9.5-gcc115-opt/bin/thisroot.sh

# Print the arguments
echo "Arguments passed to this script:"
echo "  input file  : $1"
echo "  output file : $2"
echo "  move to     : $3"

# Write PDLINPUT file for file-based input
echo "FILE=$1" > dummy

# Run the nanoAOD producer
build/delphi-nanoaod/delphi-nanoaod -P dummy --mc --config config/delphi-nanoaod.yaml --output "$2"

# Run treefy step
root -q -b -l "scripts/treefy.C+(\"$2\")"

# Move the output ROOT file to desired location
mv nanoaod_* $3



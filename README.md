# DELPHI NanoAOD

TThe DELPHI NanoAOD is built on the DELPHI [SKELANA Analysis Framework](https://opendata-qa.cern.ch/record/80502), written in Fortran. SKELANA processes data from the DELPHI Full, Short, and XShort DST formats, storing extracted information in COMMON blocks for further analysis. Various DELPHI analysis programs are integrated, allowing users to control data processing and re-execute selected algorithms to study the impact of parameter adjustments and fine-tuning.

DELPHI data is structured around a particle-based model, where each object is represented as a particle with associated blocklets containing detector-specific information. These blocklets include tracking data, calorimeter measurements, and particle identification details.

At the core of SKELANA are COMMON blocks, such as `VECP(10, NVECP)`, where each column corresponds to a specific physical quantity, including momentum components (`px`, `py`, `pz`), mass, energy, absolute momentum (`|p|`), PDG ID, charge, and more. Users must correctly interpret column types, as floating-point arrays are often mapped to integer arrays. Similar array structures exist for various blocklets. The definitions and descriptions of these structures, along with control flags that steer program execution, can be found in [stdcdes.car](http://github.com/delphi/maxi/stdcdes.car).

Tools are provided to map some of these COMMON blocks to corresponding C++ structures while maintaining a syntax similar to Fortran. For example, adjustments are made to accommodate Fortran’s 1-based indexing, ensuring that `VECP(1,1)` in Fortran corresponds to `skelana::VECP(1,1)` in C++.

In the next step, this information is stored in the still experimental ROOT RNtuple format. Further details on the new NTuple format can be found in [reference to be added].

The recommended approach for analyzing this NTuple data is through the ROOT DataFrame framework, which integrates well with Python. While direct access to the data outside of DataFrames in Python is possible, the implementation is still evolving, and interfaces may change.  


# Installation on LXPLUS

The installation requires the DELPHI framework setup and access to a recent ROOT version. To avoid complication, it is recommanded to use the standard os compiler from the installation. 

Here the setuo on a EL9 installation, as on lxplus, where DELPHI software and ROOT
are provided on CVMFS.

DELPHI software
```bash
# DELPHI Software
source /cvmfs/delphi.cern.ch/setup.sh
# Recent ROOT with the system compiler
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.34.04/x86_64-almalinux9.5-gcc115-opt/bin/thisroot.sh
```

Checkout the repository and compile

```bash
git clone https://gitlab.cern.ch/delphi/delphi-nanoaod.git
cd delphi-nanoaod

# Configure make files
cmake -B build

# Build the software
cmake --build build

# Run test
python3 python/run.py

# Batch production
Edit condor/sample_list.yaml
Edit condor/submit_nano.py
python3 condor/submit_nano.py
```

# on SVMIT03

```bash
# DELPHI Software environment
source setup.sh
# produce a small sample and convert to TTree
python3 python/run.py
```

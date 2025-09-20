# DELPHI NanoAOD

TThe DELPHI NanoAOD is built on the DELPHI [SKELANA Analysis Framework](https://opendata-qa.cern.ch/record/80502), written in Fortran. SKELANA processes data from the DELPHI Full, Short, and XShort DST formats, storing extracted information in COMMON blocks for further analysis. Various DELPHI analysis programs are integrated, allowing users to control data processing and re-execute selected algorithms to study the impact of parameter adjustments and fine-tuning.

DELPHI data is structured around a particle-based model, where each object is represented as a particle with associated blocklets containing detector-specific information. These blocklets include tracking data, calorimeter measurements, and particle identification details.

At the core of SKELANA are COMMON blocks, such as `VECP(10, NVECP)`, where each column corresponds to a specific physical quantity, including momentum components (`px`, `py`, `pz`), mass, energy, absolute momentum (`|p|`), PDG ID, charge, and more. Users must correctly interpret column types, as floating-point arrays are often mapped to integer arrays. Similar array structures exist for various blocklets. The definitions and descriptions of these structures, along with control flags that steer program execution, can be found in [stdcdes.car](http://github.com/delphi/maxi/stdcdes.car).

Tools are provided to map some of these COMMON blocks to corresponding C++ structures while maintaining a syntax similar to Fortran. For example, adjustments are made to accommodate Fortran’s 1-based indexing, ensuring that `VECP(1,1)` in Fortran corresponds to `skelana::VECP(1,1)` in C++.

In the next step, this information is stored in the still experimental ROOT RNtuple format. Further details on the new NTuple format can be found in [reference to be added].

The recommended approach for analyzing this NTuple data is through the ROOT DataFrame framework, which integrates well with Python. While direct access to the data outside of DataFrames in Python is possible, the implementation is still evolving, and interfaces may change.  


---

# Running on LXPLUS

The installation requires the DELPHI framework setup and access to a recent ROOT version.  
On **EL9 lxplus**, both are available on CVMFS.

## Environment Setup

```bash
# DELPHI Software
source /cvmfs/delphi.cern.ch/setup.sh

# ROOT (recent version with system compiler)
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.34.04/x86_64-almalinux9.5-gcc115-opt/bin/thisroot.sh
```
## Build the Software
bash
```
# Clone repository
git clone https://gitlab.cern.ch/delphi/delphi-nanoaod.git
cd delphi-nanoaod

# Configure build
cmake -B build

# Build
cmake --build build

# Run test
python3 python/run.py
```

## 🚀 Running
The main executable is:
```
./build/delphi-nanoaod/delphi-nanoaod
```
#### Arguments:

| Option | Description |
| :--- | :--- |
| `-h`, `--help` | Show help message and exit |
| `-v`, `--version` | Print version information and exit |
| `-N`, `--nickname NICKNAME` | FATMEN nickname of the sample |
| `-P`, `--pdlinput FILE` | Path to the PDL input file |
| `-C`, `--config FILE` | Path to the YAML configuration file [required] |
| `-O`, `--output FILE` | Output file name [required] |
| `-m`, `--max-events N` | Maximum number of events to process |
| `--mc` | Write simulation information |
| `--data` | Do not write simulation information (default) |

> **⚠️ Important:** At least one of `--nickname` or `--pdlinput` must be specified.

**Input Types:**
- **Nickname:** Refers to DELPHI FATMEN sample names
- **PDL input:** Should follow the format `./dummy`
  
## 🏭 Batch Production

#### 📋 Sample Lists:
Edit condor/sample_list.yaml to look up and add samples.

These lists map each DELPHI FATMEN sample nickname to its /eos path, so that a job is submitted per file.
As of Sep 20, 2025, only some data and MC samples (1994, 1995, 1999) are mapped.
**You will need to extend this mapping to include more samples.**

> 🔗 **Complete sample lists are available:**

[DELPHI CASTOR catalog](https://delphi-www.web.cern.ch/delphi-www//offline/data/castor/html/index.html)

[DELPHI Open Data Portal](https://opendata.cern.ch/search?q=DELPHI&l=list&order=asc&p=1&s=10&sort=bestmatch)

📝 Note: Each year’s samples have multiple simulation campaigns. Only the recommended ones are needed for data analysis (e.g., 1994 → v94c2 is recommended).

#### 📤 Submitting Jobs:
Edit condor/submit_nano.py to specify the FATMEN nickname you want to process. 

```
python3 condor/submit_nano.py
```
📝 Note: You will need to update the relevant path in condor/run_condor.sh and condor/submit_nano.py to match your input and output path. 

#### ⚡ Quick start example: 
The default setting converts the DELPHI official tautau sample to ROOT. 
```
python3 condor/submit_nano.py
```
📝 Note: You will need to update the relevant path in condor/run_condor.sh and condor/submit_nano.py to match your input and output path. 

# on SVMIT03

```bash
# DELPHI Software environment
source setup.sh
# produce a small sample and convert to TTree
python3 python/run.py
```

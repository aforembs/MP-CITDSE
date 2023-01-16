# ParTDSE

A set of ab-initio parallel programs for the solution of the Time-Dependent Schrödinger Equation (TDSE) in the case of one- and/or two-electron atom-laser interactions.

The included programs are, in order of operation: \
`h1e`      - Solves the 1-electron Time-Independent Schrödinger Equation (TISE) on a B-splines basis\
`w1e`      - Returns the wave function values and integration parameters for subsequent programs\
`dmx1e`    - Computes the 1-electron dipole matrices\
`tdse1e`   - Solves the 1-electron TDSE for a given laser pulse\
`pes-td1`  - Computes the photoelectron energy distribution from the 1-electron time dependent coefficients\
`gen2eidx` - Indexes the configurations provided for each total (2-electron) angular momentum\
`r12`      - Calculates the inter-electronic correlation for the 2-electron states\
`dmx2e`    - Computes the 2-electron dipole matrices using the 1-electron dipoles\
`cibasis`  - Forms the Configuration Interaction (CI) 2-electron basis\
`tdse2e`   - Solves the 2-electron TDSE for a given laser pulse on the CI basis\
`pes-td2`  - Computes the photoelectron energy distribution from the 2-electron time dependent coefficients

## Requirements

A C++17 or greater compatible compiler \
tbb                       - library for C++17 multithreading (apt-get, pacman, etc...)  \
hdf5                      - data storage library (apt-get, pacman, etc...) \
yaml-cpp                  - library for human readable (yaml) format configuration files \
any parallel BLAS library - e.g. intel MKL, OpenBLAS \
wigxjpf                   - library for wigner symbols; cmake downloads and links it automatically, if using make install from github:

```
[[ -d lib ]] || mkdir lib && cd lib
git clone https://github.com/nd-nuclear-theory/wigxjpf.git
cd wigxjpf && make
```

## Installation

cmake compilation:

```
mkdir build && cd build
cmake ..
make
```

make compilation:

```
cd src
make
```

## Running a 1-electron example

## Running a 2-electron example

To run go back to the main directory and run

```
bin/run_structure.sh
bin/tdse2e -f inp/general.yaml
```

All input parameters can be found and modified in inp/general.yaml
The 2e configurations are located in the inp/cfg-\<L\>.inp files

To produce the energy spectra run (for now this only works with bin/tdse2e not bin/td*diag)
`bin/pes2e -f inp/general -i dat/he_ct*<time>.dat`

All data and output files can be found in
`dat/`

`he_ct_<time>.dat` - contain the coefficients

`he_pop.dat` - contains the ground state population

`he_pes<L>.dat` - contains the energy spectrum for a given angular momentum L

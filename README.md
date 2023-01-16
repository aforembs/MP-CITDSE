# ParTDSE

## Decsription

Parallel programs for the time propagation of 1- and 2-electron systems

## Requirements

A C++17 or greater compatible compiler
tbb - library for C++17 multithreading (apt-get, pacman, etc...)
hdf5 - data storage library (apt-get, pacman, etc...)
yaml-cpp - library for human readable (yaml) format configuration files
any parallel BLAS library - e.g. intel MKL, OpenBLAS
wigxjpf - library for wigner symbols; cmake downloads and links it automatically, if using make install from github:
`[[ -d lib ]] || mkdir lib && cd lib`
`git clone https://github.com/nd-nuclear-theory/wigxjpf.git`
`cd wigxjpf && make`

## Installation

cmake compilation:
`mkdir build && cd build`
`cmake ..`
`make`

make compilation:
`cd src`
`make`

## Running a 1-electron example

## Running a 2-electron example

To run go back to the main directory and run
`bin/run_structure.sh`
`bin/tdse2e -f inp/general.yaml`

All input parameters can be found and modified in inp/general.yaml
The 2e configurations are located in the inp/cfg-\<L\>.inp files

To produce the energy spectra run (for now this only works with bin/tdse2e not bin/td*diag)
`bin/pes2e -f inp/general -i dat/he_ct*<time>.dat`

All data and output files can be found in
`dat/`
`he_ct_<time>.dat` - contain the coefficients
`he_pop.dat` - contains the ground state population
`he_pes<L>.dat` - contains the energy spectrum for a given angular momentum L

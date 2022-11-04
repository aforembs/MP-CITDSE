# 2e-DMX
2 electron dipole matrices from 1 electron partial waves

## Requirements
A C++17 or greater compatible compiler
tbb  - library for C++17 multithreading (apt-get, pacman, etc...)
hdf5 - data storage library (apt-get, pacman, etc...)


## Installation 
To compile do 
`cd src`
`make`

To run go back to the main directory and run
`bin/run_structure.sh`
`bin/tdse2e -f inp/general.yaml`

To produce the energy spectra run (for now this only works with bin/tdse2e not bin/td\_diag)
`bin/pes2e -f inp/general -i dat/he_ct_<time>.dat`

All data and output files can be found in 
`dat/`
`he_ct_<time>.dat` contain the coefficients
`he_pop.dat` contains the ground state population

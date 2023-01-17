# ParTDSE

A set of parallel programs for the ab-initio solution of the Time-Dependent Schrödinger Equation (TDSE) in the case of 1- and/or 2-electron atom-laser interactions. Initially the 1-electron Time-Independent Schrödinger Equation (TISE) is solved on a basis of B-splines, as described in [1](#appl_B-spl). The 1-electron dipole transition elements are then obtained, using the 1-electron eigenstates of the TISE. From here, solution of the 1-electron TDSE is possible. \
Moving onto the 2-electron problem, the inter-electronic interactions and 2-electron dipole transition elements, are calculated before forming the Configuration Interaction (CI) basis. Once the 2-electron time dependent Hamiltonian is obtained in the CI basis the 2-electron TDSE ma be solved.\
As a form of result analysis this suite provides two programs for the calculation of photoelectron spectra from time dependent coefficients, one for the 1-electron, and another for the 2-electron problem. These programs also provide the ground state populations and ionisation yields.

#### The included programs are, in order of operation:
`h1e`      - Solves the 1-electron TISE on a B-splines basis\
`w1e`      - Returns the wave function values and integration parameters for subsequent programs\
`dmx1e`    - Computes the 1-electron dipole matrices\
`tdse1e`   - Solves the 1-electron TDSE for a given laser pulse\
`pes-td1`  - Computes the photoelectron energy distribution from the 1-electron time dependent coefficients\
`gen2eidx` - Indexes the configurations provided for each total (2-electron) angular momentum\
`r12`      - Calculates the inter-electronic interactions for the 2-electron states\
`dmx2e`    - Computes the 2-electron dipole matrices using the 1-electron dipoles\
`cibasis`  - Forms the CI 2-electron basis\
`tdse2e`   - Solves the 2-electron TDSE for a given laser pulse on the CI basis\
`pes-td2`  - Computes the photoelectron energy distribution from the 2-electron time dependent coefficients

## Requirements

A C++17 or greater compatible compiler (eg. g++, clang++, ... etc.) \
tbb                       - library for C++17 multithreading (apt-get, pacman, etc...)  \
hdf5                      - data storage library (apt-get, pacman, etc...) \
yaml-cpp                  - library for human readable (yaml) format configuration files \
any parallel BLAS library - e.g. intel MKL, OpenBLAS \
[wigxjpf](https://github.com/nd-nuclear-theory/wigxjpf)                   - library for wigner symbols; the cmake compilation downloads and links it automatically, if using make you need to download it from github and compile manually:

```
[[ -d lib ]] || mkdir lib && cd lib
git clone https://github.com/nd-nuclear-theory/wigxjpf.git
cd wigxjpf && make
```

## Compilation
#### Cmake
Requires Cmake version 3.20 or greater.\
To compile the code, including the wigxjpf library:
```
mkdir build && cd build
cmake -DBLA_VENDOR=<chosen blas library (default=OpenBLAS)> ..
make
```

If your BLAS library is installed in a non-standard location do:
```
cmake -DBLA_VENDOR=<blas library> -DBLA_HINT=<path to blas library> -DBLA_INC=<path to blas/lapack headers> ..
```
instead of the above cmake call.

#### Make
First set the variables BLAS_LIB, BLAS_INC, WIG_LIB, WIG_INC to the correct paths in src/Makefile.\
Then simpy run:
```
cd src
make
```

## Running a 1-electron example
```
[[ -d dat ]] || mkdir dat
bin/h1e -f inp/H_test.yaml
bin/w1e -f inp/H_test.yaml
bin/dmx1e -f inp/H_test.yaml
mkdir dat/H_tdat
bin/tdse1e -f inp/H_test.yaml -o dat/H_tdat
```
Text files contatining the time dependent coefficients will then be available in dat/H_tdat.\
To obtain the photoelectron energy distribution run:
```
bin/pes-td1 -f inp/H_test.yaml -i dat/H_tdat/h_ct_<time>.dat
```
This will produce a photoelectron distribution in the file dat/H_tdat/h_pes.dat

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

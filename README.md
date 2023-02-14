# ParTDSE

A set of parallel programs for the ab-initio solution of the Time-Dependent Schrödinger Equation (TDSE) in the case of 1- and/or 2-electron atom-laser interactions. In many ways this package is a spiritual successor to \[[1](<https://doi.org/10.1016/S0010-4655(02)00684-7>)\].

Initially the 1-electron Time-Independent Schrödinger Equation (TISE) is solved on a basis of B-splines, as described in \[[2](https://iopscience.iop.org/article/10.1088/0034-4885/64/12/205)\]. The 1-electron dipole transition elements are then obtained, using the 1-electron eigenstates of the TISE. From here, solution of the 1-electron TDSE is possible. \
Moving onto the 2-electron problem, the inter-electronic interactions and 2-electron dipole transition elements, are calculated before forming the Configuration Interaction (CI) basis. Once the 2-electron time dependent Hamiltonian is obtained in the CI basis the 2-electron TDSE ma be solved.\
As a form of result analysis this suite provides two programs for the calculation of photoelectron spectra from time dependent coefficients, one for the 1-electron, and another for the 2-electron problem. These programs also provide the ground state populations and ionisation yields.

#### The included programs are, in order of operation:

`h1e` - Solves the 1-electron TISE on a B-splines basis\
`w1e` - Returns the wave function values and integration parameters for subsequent programs\
`dmx1e` - Computes the 1-electron dipole matrices\
`gen2eidx` - Indexes the configurations provided for each total (2-electron) angular momentum\
`r12` - Calculates the inter-electronic interactions for the 2-electron states\
`dmx2e` - Computes the 2-electron dipole matrices using the 1-electron dipoles\
`cibasis` - Forms the CI 2-electron basis\
`tdse` - Solves either the '-e 1' (1-electron) or '-e 2' (2-electron) TDSE for a given laser pulse\
`pes` - Computes the photoelectron energy distribution from the provided '-e 1' or '-e 2' time dependent coefficients

## Requirements

A C++17 or greater compatible compiler (eg. g++, clang++, ... etc.) \
tbb - library for C++17 multithreading (apt-get, pacman, etc...) \
hdf5 - data storage library (apt-get, pacman, etc...) \
yaml-cpp - library for human readable (yaml) format configuration files (apt-get, pacman, etc...) \
boost - specifically ::odeint for the tdse program (apt-get, pacman, etc...) \
any parallel BLAS library - e.g. intel MKL, OpenBLAS \
[wigxjpf](https://github.com/nd-nuclear-theory/wigxjpf) - library for wigner symbols \[[3](https://doi.org/10.1137/15M1021908)\]; the cmake compilation downloads and links it automatically, if using make you need to download it from github and compile manually:

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

If you are installing on Debian you may need to add the location of libhdf5.so to your PATH, you can then include the headers by using the:

```
-DH5_INC=<path/to/hdf5/include>
```

#### Make

First set the variables BLAS_LIB, BLAS_INC, WIG_LIB, WIG_INC to the correct paths in src/Makefile. If you're using some non-standard install locations, you may also need to set: HDF5_LIB, HDF5_INC, YAML_LIB, YAML_INC, LAPACKE_LIB and LAPACKE_INC.\
Then simpy run:

```
cd src
make
```

## Running a 1-electron example

#### Structure

To run the electronic structure programs you can use the provided bin/run-structure.sh script. \
For the 1-electron case do:

```
bin/run-structure.sh -e 1 -f inp/H_test.yaml
```

Alternatively you can run the programs individually

```
bin/h1e -f inp/H_test.yaml
bin/w1e -f inp/H_test.yaml
bin/dmx1e -f inp/H_test.yaml
```

#### TDSE Solution

```
bin/tdse -e 1 -f inp/H_test.yaml -o dat/H_tdat
```

If you wish to track the population of a specific state at each time step then use the optional parameters `-n <primary quantum number> -l <angular quantum number>`:

```
bin/tdse -e 1 -f inp/H_test.yaml -o dat/H_tdat -n 1 -l 1
```

Text files contatining the time dependent coefficients will then be available in dat/H_tdat. The folder "dat/H_tdat" will be created at runtime, any UNIX-valid folder name can be used. \
To obtain the photoelectron energy distribution run:

```
bin/pes -e 1 -f inp/H_test.yaml -i dat/H_tdat/h_ct_<time>.dat
```

By default this creates separate energy spectra for each 'l' in "dat/H_tdat/h_pes<l>.dat" files. A a total spectrum can also be obtained by using the -s option.

```
bin/pes -e 1 -f inp/H_test.yaml -i dat/H_tdat/h_ct_<time>.dat -s
```

which creates a "h_pes.dat" file.

## Running a 2-electron example

#### Structure

To run the electronic structure programs you can use the provided bin/run-structure.sh script. \
For the 2-electron case do:

```
bin/run-structure.sh -e 2 -f inp/He_test.yaml -i inp/
```

Alternatively you can run the programs individually

```
bin/h1e -f inp/He_test.yaml
bin/w1e -f inp/He_test.yaml
bin/dmx1e -f inp/He_test.yaml
bin/gen2eidx -f inp/He_test.yaml -i inp/
bin/r12 -f inp/He_test.yaml -i inp/
bin/dmx2e -f inp/He_test.yaml -i inp/
bin/cibasis -f inp/He_test.yaml
```

#### TDSE Solution

```
bin/tdse -e 2 -f inp/general.yaml -o dat/He_tdat
```

If you wish to track the population of a specific state at each time step then use the optional parameters `-n <primary quantum number> -l <angular quantum number>`:

```
bin/tdse -e 2 -f inp/general.yaml -o dat/He_tdat -n 1 -l 1
```

Text fil.es contatining the time dependent coefficients will then be available in dat/H_tdat.\
To obtain the photoelectron energy distribution run:

```
bin/pes -e 2 -f inp/He_test.yaml -i dat/He_tdat/he_ct_<time>.dat
```

By default this creates separate energy spectra for each 'l' in "dat/He_tdat/he_pes<l>.dat" files. A a total spectrum can also be obtained by using the -s option.

```
bin/pes -e 1 -f inp/He_test.yaml -i dat/He_tdat/he_ct_<time>.dat -s
```

which creates a "he_pes.dat" file.

## I/O

All input parameters can be found in the supplied inp/H_test.yaml and inp/He_test.yaml files. \
The 2e configurations are located in the inp/cfg-\<L\>.inp files. \
To provide custom inputs youcan just copy the inp/ folder and rename/modify as needed

All structure calculation files can be found in the dat directory.\
Most of the data is stored in the HDF5 format.

## References

[1] Nikolopoulos L. A., A package for the ab-initio calculation of one-and
two-photon cross sections of two-electron atoms, using a CI B-splines method.
Computer physics communications. 2003 Feb 1;150(2):140-65.

[2] Bachau H, Cormier E, Decleva P, Hansen JE, Martín F. Applications of
B-splines in atomic and molecular physics. Reports on progress in physics.
2001 Nov 19;64(12):1815.

[3] Johansson H. T. and Forssén C., Fast and Accurate Evaluation of Wigner 3j,
6j, and 9j Symbols Using Prime Factorization and Multiword Integer Arithmetic,
SIAM J. Sci. Comput., 38(1) (2016), A376-A384.

## Cite as (Add once published)

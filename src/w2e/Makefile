ifndef CPPC
CPPC = g++
endif

EXEC = ../../bin/w2e

IDIR = -Iinclude
SH_IDIR = -I../include
ODIR = obj

SRC =cibasis_main.cpp cibasis.cpp

ifndef BLAS_LIB
BLAS_LIB = -lopenblas -L/opt/OpenBLAS/lib
endif

ifndef BLAS_INC
BLAS_INC = -I/opt/OpenBLAS/include
endif

ifndef HDF5_LIB
HDF5_LIB = -lhdf5 -lhdf5_cpp
endif

ifndef YAML_LIB
YAML_LIB = -lyaml-cpp
endif

ifndef LAPACK_LIB
LAPACK_LIB = -llapack
endif

CPPFLAGS = -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -O3 -march=native -std=c++17 \
 -fopenmp $(IDIR) $(SH_IDIR) $(BLAS_INC) $(HDF5_INC) $(YAML_INC) $(LAPACKE_INC)
LDFLAGS = $(HDF5_LIB) $(YAML_LIB) $(LAPACK_LIB) $(BLAS_LIB) -fopenmp
OBJ = $(patsubst %.cpp,$(ODIR)/%.o,$(SRC))

all: $(SRC) $(EXEC)

$(EXEC): $(OBJ)
	$(CPPC) $(OBJ) $(LDFLAGS) -o $@

$(ODIR)/%.o: %.cpp
	$(CPPC) $(CPPFLAGS) -c -o $@ $<

../$(ODIR)/%.o: %.cpp
	$(CPPC) $(CPPFLAGS) -c -o $@ $<

.PHONY : clean
clean :
	rm -f $(EXEC) $(OBJ)

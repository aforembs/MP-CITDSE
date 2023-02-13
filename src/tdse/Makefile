ifndef CPPC
CPPC = g++
endif

EXEC = ../../bin/tdse

IDIR = -Iinclude
SH_IDIR = -I../include
ODIR = obj

VPATH := ../

SRC =tdse_main.cpp td_read.cpp tdse.cpp pulse.cpp
SRC_SH=

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

CPPFLAGS = -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -O3 -march=native -std=c++17 \
-fopenmp $(IDIR) $(SH_IDIR) $(BLAS_INC) $(HDF5_INC) $(YAML_INC)
LDFLAGS = $(HDF5_LIB) $(YAML_LIB) -ltbb -fopenmp $(BLAS_LIB)
OBJ = $(patsubst %.cpp,$(ODIR)/%.o,$(SRC))
OBJ+=$(join $(addsuffix ../obj/, $(dir $(SRC_SH))), $(notdir $(SRC_SH:.cpp=.o)))

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
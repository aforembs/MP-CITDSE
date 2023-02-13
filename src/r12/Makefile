ifndef CPPC
CPPC = g++
endif

EXEC = ../../bin/r12

IDIR = -Iinclude
SH_IDIR = -I../include
ODIR = obj

VPATH := ../

SRC = r12_main.cpp r12.cpp integrator.cpp
SRC_SH = fastgl.cpp cfg_in.cpp

ifndef BLAS_LIB
BLAS_LIB = -lopenblas -L/opt/OpenBLAS/lib
endif

ifndef WIG_LIB
WIG_LIB = -lwigxjpf -L../../lib
endif

ifndef WIG_INC
WIG_INC = -I../../lib/include
endif

ifndef HDF5_LIB
HDF5_LIB = -lhdf5 -lhdf5_cpp
endif

ifndef YAML_LIB
YAML_LIB = -lyaml-cpp
endif

CPPFLAGS = -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -O3 -march=native -std=c++17 \
-fopenmp $(IDIR) $(SH_IDIR) $(WIG_INC) $(HDF5_INC) $(YAML_INC)
LDFLAGS = $(HDF5_LIB) $(YAML_LIB) -ltbb -fopenmp -L$(VPATH) $(WIG_LIB) -lgsl $(BLAS_LIB)
OBJ = $(patsubst %.cpp,$(ODIR)/%.o,$(SRC))
OBJ+=$(join $(addsuffix ../$(ODIR)/, $(dir $(SRC_SH))), $(notdir $(SRC_SH:.cpp=.o)))

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
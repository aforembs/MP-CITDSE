ifndef CPPC
CPPC = g++
endif

EXEC = ../../bin/dmx1e

IDIR = -Iinclude
SH_IDIR = -I../include
ODIR = obj

VPATH := ../

SRC =dmx1e_main.cpp dmx1e.cpp integrator.cpp
SRC_SH=fastgl.cpp

ifndef HDF5_LIB
HDF5_LIB = -lhdf5 -lhdf5_cpp
endif

ifndef YAML_LIB
YAML_LIB = -lyaml-cpp
endif

CPPFLAGS = -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -O3 -march=native -std=c++17 \
-fopenmp $(IDIR) $(SH_IDIR) $(HDF5_INC) $(YAML_INC)
LDFLAGS = $(HDF5_LIB) $(YAML_LIB) -fopenmp
OBJ = $(patsubst %.cpp,$(ODIR)/%.o,$(SRC))
#OBJ+= $(addprefix ../, $(patsubst %.cpp,$(ODIR)/%.o,$(SRC_SHARED)))
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
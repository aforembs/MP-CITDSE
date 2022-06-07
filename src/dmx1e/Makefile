CPPC = g++
EXEC = ../../bin/dmx1e

IDIR = include
SH_IDIR = ../include
ODIR = obj

VPATH := ../

SRC =dmx1e_main.cpp dmx1e.cpp 
SRC_SH=fastgl.cpp bsplines.cpp

CPPFLAGS = -Wall -pedantic -O3 -march=native -std=c++17 -fopenmp -I$(IDIR) -I$(SH_IDIR) -I/opt/OpenBLAS/include -g
LDFLAGS = -lhdf5 -lhdf5_cpp -ltbb -lyaml-cpp -llapacke -fopenmp -L../ -lslatec -lgfortran #-lopenblas -L/opt/OpenBLAS/lib
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
	rm $(EXEC) $(OBJ)
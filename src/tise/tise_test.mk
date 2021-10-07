CPPC = g++
EXEC = ../bin/h1e_test#dmx_test

IDIR = include
SH_IDIR = ../include
ODIR = obj

SRC =tise_test.cpp tise.cpp ../fastgl.cpp ../bsplines.cpp

CPPFLAGS = -Wall -pedantic -O3 -march=native -std=c++17 -fopenmp -I$(IDIR) -I$(SH_IDIR) -I/opt/OpenBLAS/include -g
LDFLAGS = -lhdf5 -lhdf5_cpp -lyaml-cpp -llapacke -fopenmp -L. -lslatec -lgfortran #-lopenblas -L/opt/OpenBLAS/lib
OBJ = $(patsubst %.cpp,$(ODIR)/%.o,$(SRC))

all: $(SRC) $(EXEC)

$(EXEC): $(OBJ)
	$(CPPC) $(OBJ) $(LDFLAGS) -o $@

$(ODIR)/%.o: %.cpp
	$(CPPC) $(CPPFLAGS) -c -o $@ $<

.PHONY : clean
clean :
	rm $(EXEC) $(OBJ)
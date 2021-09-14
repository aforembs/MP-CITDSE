CPPC = g++
EXEC = V12_test#dmx_test

IDIR = include
ODIR = obj

SRC =V12_test.cpp V_12.cpp wigner_sym.cpp fastgl.cpp bsplines.cpp
#dmx_test.cpp dmx_calc.cpp wigner6j.cpp

CPPFLAGS = -Wall -pedantic -O3 -march=native -std=c++17 -fopenmp -I$(IDIR) -g
LDFLAGS = -lhdf5 -lhdf5_cpp -ltbb -fopenmp -L. -lslatec -lgfortran
OBJ = $(patsubst %.cpp,$(ODIR)/%.o,$(SRC))

all: $(SRC) $(EXEC)

$(EXEC): $(OBJ)
	$(CPPC) $(OBJ) $(LDFLAGS) -o $@

$(ODIR)/%.o: %.cpp
	$(CPPC) $(CPPFLAGS) -c -o $@ $<

.PHONY : clean
clean :
	rm $(EXEC) $(OBJ)
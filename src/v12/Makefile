CPPC = g++
EXEC = ../../bin/V12_test

IDIR = include
SH_IDIR = ../include
ODIR = obj

VPATH := ../

SRC = V12_test.cpp V_12.cpp 
SRC_SH = wigner_sym.cpp fastgl.cpp bsplines.cpp

CPPFLAGS = -Wall -pedantic -O3 -march=native -std=c++17 -fopenmp -I$(IDIR) -I$(SH_IDIR) -g
LDFLAGS = -lhdf5 -lhdf5_cpp -ltbb -fopenmp -L../ -lslatec -lgfortran -g
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
	rm $(EXEC) $(OBJ)
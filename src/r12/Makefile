CPPC = g++
EXEC = ../../bin/r12

IDIR = include
SH_IDIR = ../include
ODIR = obj

VPATH := ../

SRC = r12_main.cpp r12.cpp integrator.cpp
SRC_SH = fastgl.cpp bsplines.cpp cfg_in.cpp

CPPFLAGS = -Wall -pedantic -O3 -march=native -std=c++17 -fopenmp -I$(IDIR) -I$(SH_IDIR) -I../../lib/wigxjpf/inc -g
LDFLAGS = -lhdf5 -lhdf5_cpp -ltbb -lyaml-cpp -fopenmp -L$(VPATH) -L../../lib/wigxjpf/lib -lwigxjpf -lslatec -lgfortran -g
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
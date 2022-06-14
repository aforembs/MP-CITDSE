CPPC = g++
EXEC = ../../bin/dmx2e

IDIR = include
SH_IDIR = ../include
ODIR = obj

VPATH := ../

SRC = dmx2e_main.cpp dmx2e.cpp
SRC_SH = wigner_sym.cpp cfg_in.cpp

CPPFLAGS = -Wall -pedantic -O3 -march=native -std=c++17 -fopenmp -I$(IDIR) -I$(SH_IDIR) -I../../lib/wigxjpf/inc
LDFLAGS = -lhdf5 -lhdf5_cpp -ltbb -lyaml-cpp -L$(VPATH) -L../../lib/wigxjpf/lib -lwigxjpf -fopenmp
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
ifndef CPPC
CPPC = g++
endif

EXEC = ../../bin/w1e

IDIR = include
SH_IDIR = ../include
ODIR = obj

VPATH := ../

SRC = w1e_main.cpp w1e.cpp
SRC_SH = fastgl.cpp

ifndef BLAS_LIB
BLAS_LIB = -lopenblas -L/opt/OpenBLAS/lib
endif

CPPFLAGS = -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -O3 -march=native -std=c++17 -I$(IDIR) -I$(SH_IDIR) -fopenmp
LDFLAGS = -lhdf5 -lhdf5_cpp -lyaml-cpp -L$(VPATH) -lgsl $(BLAS_LIB) -fopenmp
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
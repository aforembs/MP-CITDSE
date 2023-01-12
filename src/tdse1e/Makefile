CPPC = g++
EXEC = ../../bin/tdse1e

IDIR = include
SH_IDIR = ../include
ODIR = obj

VPATH := ../

SRC = tdse_main.cpp tdse.cpp td_read.cpp pulse.cpp
SRC_SH = 

ifndef BLAS_LIB
BLAS_LIB = -lopenblas -L/opt/OpenBLAS/lib
endif

ifndef BLAS_INC
BLAS_INC = /opt/OpenBLAS/include
endif

CPPFLAGS = -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -O3 -march=native -std=c++17 -I$(IDIR) -I$(SH_IDIR) -I$(BLAS_INC) -fopenmp
LDFLAGS = -llapack -lyaml-cpp -lhdf5 -lhdf5_cpp -fopenmp $(BLAS_LIB)

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
	rm $(EXEC) $(OBJ)

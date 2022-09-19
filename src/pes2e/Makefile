CPPC = g++
EXEC = ../../bin/tdse2e

IDIR = include
SH_IDIR = ../include
ODIR = obj

VPATH := ../

SRC =pes2e_main.cpp pes2e.cpp
SRC_SH=

CPPFLAGS = -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -O3 -march=native -std=c++17 -fopenmp -I$(IDIR) -I$(SH_IDIR) -g
LDFLAGS = -lhdf5 -lhdf5_cpp -ltbb -lyaml-cpp -fopenmp
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
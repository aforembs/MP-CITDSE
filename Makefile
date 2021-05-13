CPPC = g++
EXEC = dmx_test

IDIR = include
ODIR = obj

SRC = dmx_test.cpp dmx_calc.cpp wigner6j.cpp

CPPFLAGS = -Wall -pedantic -O3 -march=native -std=c++11  -I$(IDIR) 
LDFLAGS = -lhdf5 -lhdf5_cpp
OBJ = $(patsubst %.cpp,$(ODIR)/%.o,$(SRC))

all: $(SRC) $(EXEC)

$(EXEC): $(OBJ)
	$(CPPC) $(LDFLAGS) $(OBJ) -o $@

$(ODIR)/%.o: %.cpp
	$(CPPC) $(CPPFLAGS) -c -o $@ $<

.PHONY : clean
clean :
	rm $(EXEC) $(OBJ)
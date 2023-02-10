BIN_DIR = ../bin
DAT_DIR = ../dat
OBJ_DIR = obj

CPPC = g++

BLAS_LIB = -lopenblas -L/opt/OpenBLAS/lib -Wno-c99-extensions
BLAS_INC = /opt/OpenBLAS/include -Wno-c99-extensions

WIG_LIB = -lwigxjpf -L../../lib
WIG_INC = ../../lib/include

HDF5_LIB = -lhdf5 -lhdf5-cpp
HDF5_INC = 

YAML_LIB = -lyaml-cpp
YAML_INC = 

LAPACKE_LIB = -llapack
LAPACKE_INC = 

export

all:
	[[ -d $(BIN_DIR) ]] || mkdir $(BIN_DIR)
	[[ -d $(DAT_DIR) ]] || mkdir $(DAT_DIR)
	[[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR)
	cd h1e/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd w1e/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd dmx1e/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd gen2eidx/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd r12/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd dmx2e/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd cibasis/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd tdse/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd pes/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)

clean:
	$(MAKE) -C h1e/ clean
	$(MAKE) -C w1e/ clean
	$(MAKE) -C dmx1e/ clean
	$(MAKE) -C gen2eidx/ clean
	$(MAKE) -C r12/ clean
	$(MAKE) -C dmx2e/ clean
	$(MAKE) -C cibasis/ clean
	$(MAKE) -C tdse/ clean
	$(MAKE) -C pes/ clean
	rm -f $(OBJ_DIR)/*
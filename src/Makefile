BIN_DIR = ../bin
DAT_DIR = ../dat
OBJ_DIR = obj

LIB_DIR = ../../lib

CPPC = g++ #for Intel MKL headers add -D__INTEL_MKL__

BLAS_LIB = -lopenblas -L/opt/OpenBLAS/lib -Wno-c99-extensions
BLAS_INC = -I/opt/OpenBLAS/include -Wno-c99-extensions

WIG_LIB = -lwigxjpf -L$(LIB_DIR)/wigxjpf/lib
WIG_INC = -I$(LIB_DIR)/wigxjpf/include/

HDF5_LIB = -lhdf5 -lhdf5_cpp -L/usr/lib/x86_64-linux-gnu/hdf5/serial
HDF5_INC = -I/usr/include/hdf5/serial

YAML_LIB = -lyaml-cpp
YAML_INC = 

LAPACK_LIB = -llapacke
LAPACK_INC = 

export

all:
	[[ -d $(BIN_DIR) ]] || mkdir -p $(BIN_DIR)
	[[ -d $(DAT_DIR) ]] || mkdir -p $(DAT_DIR)
	[[ -d $(OBJ_DIR) ]] || mkdir -p $(OBJ_DIR)
	cd h1e/ && [[ -d $(OBJ_DIR) ]] || mkdir -p $(OBJ_DIR) && $(MAKE)
	cd w1e/ && [[ -d $(OBJ_DIR) ]] || mkdir -p $(OBJ_DIR) && $(MAKE)
	cd d1e/ && [[ -d $(OBJ_DIR) ]] || mkdir -p $(OBJ_DIR) && $(MAKE)
	cd id2ec/ && [[ -d $(OBJ_DIR) ]] || mkdir -p $(OBJ_DIR) && $(MAKE)
	cd r12/ && [[ -d $(OBJ_DIR) ]] || mkdir -p $(OBJ_DIR) && $(MAKE)
	cd d2e/ && [[ -d $(OBJ_DIR) ]] || mkdir -p $(OBJ_DIR) && $(MAKE)
	cd w2e/ && [[ -d $(OBJ_DIR) ]] || mkdir -p $(OBJ_DIR) && $(MAKE)
	cd tdse/ && [[ -d $(OBJ_DIR) ]] || mkdir -p $(OBJ_DIR) && $(MAKE)
	cd pes/ && [[ -d $(OBJ_DIR) ]] || mkdir -p $(OBJ_DIR) && $(MAKE)

clean:
	$(MAKE) -C h1e/ clean
	$(MAKE) -C w1e/ clean
	$(MAKE) -C d1e/ clean
	$(MAKE) -C id2ec/ clean
	$(MAKE) -C r12/ clean
	$(MAKE) -C d2e/ clean
	$(MAKE) -C w2e/ clean
	$(MAKE) -C tdse/ clean
	$(MAKE) -C pes/ clean
	rm -f $(OBJ_DIR)/*
	cd h1e/ && rm -rf $(OBJ_DIR)
	cd w1e/ && rm -rf $(OBJ_DIR)
	cd d1e/ && rm -rf $(OBJ_DIR)
	cd id2ec/ && rm -rf $(OBJ_DIR)
	cd r12/ && rm -rf $(OBJ_DIR)
	cd d2e/ && rm -rf $(OBJ_DIR)
	cd w2e/ && rm -rf $(OBJ_DIR)
	cd tdse/ && rm -rf $(OBJ_DIR)
	cd pes/ && rm -rf $(OBJ_DIR)

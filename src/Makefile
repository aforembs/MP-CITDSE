BIN_DIR = ../bin
DAT_DIR = ../dat
OBJ_DIR = obj

all:
	[[ -d $(BIN_DIR) ]] || mkdir $(BIN_DIR)
	[[ -d $(DAT_DIR) ]] || mkdir $(DAT_DIR)
	[[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR)
	cd h1e/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd w1e/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd dmx1e/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd sort2eEn/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd r12/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd dmx2e/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd gr2e/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd tdse2e/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)
	cd tdse_diag/ && [[ -d $(OBJ_DIR) ]] || mkdir $(OBJ_DIR) && $(MAKE)

clean:
	$(MAKE) -C h1e/ clean
	$(MAKE) -C w1e/ clean
	$(MAKE) -C dmx1e/ clean
	$(MAKE) -C sort2eEn/ clean
	$(MAKE) -C r12/ clean
	$(MAKE) -C dmx2e/ clean
	$(MAKE) -C gr2e/ clean
	$(MAKE) -C tdse2e/ clean
	$(MAKE) -C tdse_diag/ clean
	rm -f $(OBJ_DIR)/*
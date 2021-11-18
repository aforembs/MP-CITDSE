#ifndef V_12_H_
#define V_12_H_

#include <vector>
#include <algorithm>
#include <iomanip>
#include "fastgl.h"
#include "H5Cpp.h"
#include "dmx_typ.h"
#include "wigner_sym.h"
#include "bsplines.h"
#include "cfg_in.h"

int V12(std::string cpot, int L_max, std::vector<uint> &N_sz);

int V12(std::string cpot, int L_max, std::string dir);

#endif // V_12_H_
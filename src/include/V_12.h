#ifndef V_12_H_
#define V_12_H_

#include <vector>
#include <algorithm>
#include "fastgl.h"
#include "H5Cpp.h"
#include "dmx_typ.h"
#include "wigner_sym.h"
#include "bsplines.h"

int V12(std::string cpot, uint L_max, std::vector<uint> &N_sz);

#endif // V_12_H_
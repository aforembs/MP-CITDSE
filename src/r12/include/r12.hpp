#ifndef V_12_H_
#define V_12_H_

#include "cfg_in.hpp"
#include "dmx_typ.hpp"
#include "fastgl.hpp"
#include "integrator.hpp"
#include "wigxjpf.h"
#include <H5Cpp.h>
#include <algorithm>
#include <iomanip>
#include <omp.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace r_12 {

int readConfig(std::string file, int &qsz, std::string &pot, int &L_max,
               std::string &k_limit, bool &lim_flag);

int r12Glob(std::string cpot, int L_max, int qsz, std::string dir,
            bool lim_flag, std::string k_limit);

} // namespace r_12

#endif // V_12_H_
#ifndef V_12_H_
#define V_12_H_

#include "bsp_gsl.hpp"
#include "cfg_in.hpp"
#include "dmx_typ.hpp"
#include "fastgl.hpp"
#include "integrator_new.hpp"
#include "wigxjpf.h"
#include <H5Cpp.h>
#include <algorithm>
#include <iomanip>
#include <omp.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace r_12 {

int readConfig(std::string file, int &glq_pt, std::string &pot, int &L_max,
               char &gauge, std::string &integrator);

int r12Glob(std::string cpot, int L_max, int glq_pt, std::string dir);

} // namespace r_12

#endif // V_12_H_
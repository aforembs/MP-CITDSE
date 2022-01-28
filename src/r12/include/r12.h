#ifndef V_12_H_
#define V_12_H_

#include <vector>
#include <algorithm>
#include <iomanip>
#include <yaml-cpp/yaml.h>
#include <omp.h>
#include "fastgl.h"
#include "H5Cpp.h"
#include "dmx_typ.h"
#include "bsplines.h"
#include "cfg_in.h"
#include "integrator.h"
#include "wigxjpf.h"

namespace r_12{

int ReadConfig(std::string file, int &glq_pt,
              std::string &pot, int &L_max,
              char &gauge, std::string &integrator);

int R12Glob3(std::string cpot, int L_max, int glq_pt, std::string dir);

int R12Trap(std::string cpot, int L_max, int glq_pt, std::string dir);

}

#endif // V_12_H_
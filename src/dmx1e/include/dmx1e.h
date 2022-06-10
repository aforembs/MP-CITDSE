#include <iostream>
#include <cmath>
#include <vector>
#include <yaml-cpp/yaml.h>
#include "fastgl.h"
#include "H5Cpp.h"
#include "dmx_typ.h"
#include "integrator.h"

namespace dmx1e{
  int ReadConfig(std::string file, int &glq_pt,
                std::string &pot, int &l_max,
                char &gauge);

  int GenDipole(std::string cpot, int l_max, 
                int glq_pt, char gauge);

}
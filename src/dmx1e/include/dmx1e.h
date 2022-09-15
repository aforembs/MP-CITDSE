#ifndef DMX1E_H_
#define DMX1E_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <yaml-cpp/yaml.h>
#include "fastgl.h"
#include "H5Cpp.h"
#include "dmx_typ.h"
#include "integrator.h"

namespace dmx1e{
  int ReadConfig(std::string file, std::string &pot, 
                int &glq_pt, char &gauge, int &l_max);

  int GenDipole(std::string cpot, int glq_pt, char gauge, int l_max);
}

#endif // DMX2E_H_
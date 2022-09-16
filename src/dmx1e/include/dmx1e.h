#ifndef DMX1E_H_
#define DMX1E_H_

#include "dmx_typ.h"
#include "fastgl.h"
#include "integrator.h"
#include <H5Cpp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace dmx1e {
int ReadConfig(std::string file, std::string &pot, int &glq_pt, char &gauge,
               int &l_max);

int GenDipole(std::string cpot, int glq_pt, char gauge, int l_max);
} // namespace dmx1e

#endif // DMX2E_H_
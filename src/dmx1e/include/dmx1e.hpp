#ifndef DMX1E_H_
#define DMX1E_H_

#include "dmx_typ.hpp"
#include "fastgl.hpp"
#include "integrator.hpp"
#include <H5Cpp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace dmx1e {
int ReadConfig(std::string file, std::string &pot, int &qsz, char &gauge,
               int &l_max);

int GenDipole(std::string cpot, int qsz, char gauge, int l_max);
} // namespace dmx1e

#endif // DMX2E_H_
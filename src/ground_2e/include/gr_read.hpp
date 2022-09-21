#ifndef GR_READ_HPP_
#define GR_READ_HPP_

#include "dmx_typ.hpp"
#include <H5Cpp.h>
#include <iostream>
#include <memory>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace grrd {
int readConfig(std::string file, std::string &pot, char &gauge, int &L0_sz,
               double &dt);

int readStructure(std::string pot, char gauge, int L0_sz,
                  std::vector<double> &ens, std::vector<double> &block);
} // namespace grrd

#endif // GR_READ_HPP_
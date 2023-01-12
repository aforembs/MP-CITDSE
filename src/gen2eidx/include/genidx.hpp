#ifndef GENIDX_HPP_
#define GENIDX_HPP_

#include "cfg_in.hpp"
#include "dmx_typ.hpp"
#include <H5Cpp.h>
#include <iostream>
#include <memory>
#include <yaml-cpp/yaml.h>

namespace genidx {
int readConfig(std::string file, std::string &pot, int &L_max);

int sortEn(std::string pot, int L_max, std::string dir);
} // namespace genidx

#endif // GENIDX_HPP_
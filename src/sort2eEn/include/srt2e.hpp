#ifndef SRT_2E_HPP_
#define SRT_2E_HPP_

#include <memory>
#include <execution>
#include <algorithm>
#include <iostream>
#include <H5Cpp.h>
#include <yaml-cpp/yaml.h>
#include "dmx_typ.hpp"
#include "cfg_in.hpp"

namespace srt2e {
  int readConfig(std::string file, std::string &pot, int &L_max);

  int sortEn(std::string pot, int L_max, std::string dir);
} // namespace srt2e

#endif // SRT_2E_HPP_
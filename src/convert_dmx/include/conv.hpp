#ifndef CONV_HPP_
#define CONV_HPP_

#include "dmx_typ.hpp"
#include <vector>
#include <iostream>
#include <cassert>
#include <fstream>
#include <H5Cpp.h>
#include <memory>
#include <yaml-cpp/yaml.h>
#include <lapacke.h>
extern "C" {
#include <cblas.h>
}


using stvupt = std::vector<std::unique_ptr<std::vector<double>>>;

namespace conv {
  int readConfig(std::string file, std::string &pot, char &gauge,
                 int &L_max, std::vector<int> &state_sz);

  int calcEvecs(std::string pot, char gauge, int L_max, 
                std::vector<int> &state_sz, stvupt &vecs);

  int readDipoles(std::string pot, char gauge, int L_max,
                      std::vector<int> &state_sz, stvupt &dipoles);

  int transDip(std::vector<int> &state_sz, stvupt &vecs, stvupt &dipoles);
} // namespace conv

#endif // CONV_HPP_
#ifndef TD_READ_HPP_
#define TD_READ_HPP_

#include "dmx_typ.h"
#include <H5Cpp.h>
#include <cassert>
#include <complex>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>
#include <yaml-cpp/yaml.h>

using stvupt = std::vector<std::unique_ptr<std::vector<double>>>;

namespace tdrd {
int readConfig(std::string file, std::string &pot, char &gauge, int &l_max,
               std::vector<int> &state_sz, int &Lanc_iter, int &num_eval,
               double &timestep, double &w, double &Io, double &cepd,
               int &cycles);

int readStructure(std::string pot, char gauge, int L_max, int &ct_sz,
                  std::vector<int> &state_sz, std::vector<int> &offs,
                  stvupt &blocks);

int readDipoles(std::string pot, char gauge, int L_max,
                std::vector<int> &state_sz, stvupt &dipoles);
} // namespace tdrd

#endif // TD_READ_HPP_
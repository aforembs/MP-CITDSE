#ifndef TD_READ_HPP_
#define TD_READ_HPP_

#include "dmx_typ.hpp"
#include <H5Cpp.h>
#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>
#include <yaml-cpp/yaml.h>
extern "C" {
#include <cblas.h>
}

using stvupt = std::vector<std::unique_ptr<std::vector<double>>>;

namespace tdrd {
int readConfig(std::string file, std::string &pot, char &gauge, int &l_max,
               std::vector<int> &state_sz, double &timestep, double &w,
               double &Io, double &cepd, int &cycles);

int readStructure(std::string pot, int L_max, int &ct_sz,
                  std::vector<int> &state_sz, std::vector<int> &offs,
                  stvupt &ens);

int readDipoles(std::string pot, char gauge, int L_max,
                std::vector<int> &state_sz, stvupt &dipoles);
} // namespace tdrd

#endif // TD_READ_HPP_
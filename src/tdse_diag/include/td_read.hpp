#ifndef TD_READ_HPP_
#define TD_READ_HPP_

#include "dmx_typ.hpp"
#include <H5Cpp.h>
#include <cassert>
#include <complex>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <lapacke.h>
extern "C" {
  #include <cblas.h>
}

using stvupt = std::vector<std::unique_ptr<std::vector<double>>>;

namespace tdrd {
int readConfig(std::string file, std::string &pot, char &gauge, int &l_max,
               std::vector<int> &state_sz, int &Lanc_iter, int &num_eval,
               double &timestep, double &w, double &Io, double &cepd,
               int &cycles);

int readStructure(std::string pot, int L_max, int &ct_sz,
                  std::vector<int> &state_sz, std::vector<int> &offs,
                  stvupt &ens, stvupt &blocks);

int readDipoles(std::string pot, char gauge, int L_max,
                std::vector<int> &state_sz, stvupt &blocks, stvupt &dipoles);

int readGrCt(std::string pot, std::vector<int> &state_sz,
             std::vector<std::complex<double>> &ct);
} // namespace tdrd

#endif // TD_READ_HPP_
#ifndef PES2E_HPP_
#define PES2E_HPP_

#include "cfg_in.hpp"
#include "dmx_typ.hpp"
#include <H5Cpp.h>
#include <cassert>
#include <complex>
#include <cstring>
#include <iostream>
#include <lapacke.h>
#include <vector>
#include <yaml-cpp/yaml.h>
extern "C" {
#include <cblas.h>
}

namespace pes2e {
int readConfig(std::string file, std::string &pot, int &L_max, int &l_max,
               std::vector<int> &state_sz);

int readCt(std::string file, std::vector<std::complex<double>> &ct);

int genPES(std::string pot, std::string dir, int L_max, int l_max,
           std::vector<int> &state_sz, std::vector<std::complex<double>> &ct);

int genPES2eb(std::string pot, int L_max, std::vector<int> &state_sz,
              std::vector<std::complex<double>> &ct);
} // namespace pes2e

#endif // PES2E_HPP_
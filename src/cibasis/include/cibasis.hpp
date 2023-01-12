#ifndef CIBASIS_HPP_
#define CIBASIS_HPP_

#include "dmx_typ.hpp"
#include <H5Cpp.h>
#include <cassert>
#include <iostream>
#include <lapacke.h>
#include <numeric>
#include <vector>
#include <yaml-cpp/yaml.h>
extern "C" {
#include <cblas.h>
}

using stvupt = std::vector<std::unique_ptr<std::vector<double>>>;

namespace cib {
int readConfig(std::string file, std::string &pot, char &gauge, int &L_max);

int formCIh0(std::string pot, int L_max, stvupt &vecs);

int formCIDipoles(std::string pot, char gauge, int L_max, stvupt &vecs);
} // namespace cib

#endif // CIBASIS_HPP_
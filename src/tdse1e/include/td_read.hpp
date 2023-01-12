#ifndef TD_READ_H_
#define TD_READ_H_

#include <H5Cpp.h>
#include <cassert>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace tdrd {
int readConfig(std::string file, std::string &pot, char &gauge, int &l_max,
               std::vector<int> &state_sz, double &timestep, double &w,
               double &Io, double &cepd, int &cycles);

int readEnergies(std::string pot, int l_max, std::vector<int> &state_sz,
                 std::vector<double> &eps, std::vector<int> &offs, int &eps_sz,
                 std::vector<double *> &eps_off);

int readDipoles(std::string pot, char gauge, int l_max,
                std::vector<int> &state_sz, std::vector<double> &dip,
                std::vector<int> &doffs, std::vector<double *> &dip_off);
} // namespace tdrd

#endif // TD_READ_H_
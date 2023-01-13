#ifndef PES_H_
#define PES_H_

#include <H5Cpp.h>
#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace pes {
int readConfig(std::string file, std::string &pot, int &l_max,
               std::vector<int> &state_sz);

int readCt(std::string file, std::vector<std::complex<double>> &ct);

int genPES(std::string pot, int l_max, std::vector<int> &state_sz,
           std::vector<std::complex<double>> &ct, std::string output);
} // namespace pes

#endif // PES_H_
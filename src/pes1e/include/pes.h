#ifndef PES_H_
#define PES_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cassert>
#include <yaml-cpp/yaml.h>
#include <H5Cpp.h>

namespace pes {
  int ReadConfig(std::string file, 
                 std::string &pot, 
                 int &l_max, 
                 std::vector<int> &state_sz);

  int ReadCt(std::string file, 
             std::vector<std::complex<double>> &ct);

  int GenPES(std::string pot, int l_max,
             std::vector<int> &state_sz,
             std::vector<std::complex<double>> &ct);
}

#endif // PES_H_
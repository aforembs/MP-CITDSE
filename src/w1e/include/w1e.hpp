#ifndef W1E_H
#define W1E_H

#include "bsp_gsl.hpp"
#include "fastgl.hpp"
#include <H5Cpp.h>
#include <numeric>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace w1e {

int ReadConfig(std::string file, int &qsz, int &R_max, int &l_max,
               std::string &pot);

int GenWfn(std::string pot, int qsz, int R_max, int l_max);

} // namespace w1e

#endif // W1E_H
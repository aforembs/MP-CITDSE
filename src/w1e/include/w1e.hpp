#ifndef W1E_H
#define W1E_H

#include "bsp_gsl.hpp"
#include "fastgl.hpp"
#include <H5Cpp.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace w1e {

int ReadConfig(std::string file, int &glq_pt, int &l_max, std::string &pot,
               std::string &integrator);

int GenWfn(std::string pot, int glq_pt, int l_max, std::string integrator);

} // namespace w1e

#endif // W1E_H
#ifndef W1E_H
#define W1E_H

#include <vector>
#include <yaml-cpp/yaml.h>
#include "fastgl.h"
#include "bsplines.h"
#include "H5Cpp.h"

namespace w1e {

int ReadConfig(std::string file, int &glq_pt, int &l_max,
              std::string &pot, std::string &integrator);

int GenWfn(std::string pot, int glq_pt, int l_max, std::string integrator);

}

#endif //W1E_H
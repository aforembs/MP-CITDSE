#ifndef INPUT_R_H_
#define INPUT_R_H_

#include <iostream>
#include <yaml-cpp/yaml.h>

namespace inp{

int ReadConfig(std::string file,
              int &n, int &k,
              std::string &grid,
              int &l_max);

}

#endif // INPUT_R_H_
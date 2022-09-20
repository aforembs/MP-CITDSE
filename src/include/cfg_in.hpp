#ifndef CFG_IN_H_
#define CFG_IN_H_

#include "H5Cpp.h"
#include "dmx_typ.hpp"
#include <algorithm>
#include <execution>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace cfg {
struct line {
  int n1;
  int l1;
  int l2;
  int n2min;
  int n2max;
};

int readCfg(std::string dir, int L, int &sym, int &ncf,
            std::vector<line> &cfgs);

int genL_idx(std::string pot, char gauge, int L_max, std::string dir);
} // namespace cfg

#endif // CFG_IN_H_
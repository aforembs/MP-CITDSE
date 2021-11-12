#ifndef CFG_IN_H_
#define CFG_IN_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <filesystem>

namespace cfg {
  int ReadCfg(std::string dir, int L, int &sym, int &ncf, std::vector<int> &cfgs);
}

#endif // CFG_IN_H_
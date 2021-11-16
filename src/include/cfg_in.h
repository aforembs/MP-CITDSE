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
  struct line {
    int n1;
    int l1;
    int l2;
    int n2min;
    int n2max;
  };

  int ReadCfg(std::string dir, int L, int &sym, int &ncf, std::vector<line> &cfgs);
}

#endif // CFG_IN_H_
#include "cfg_in.h"

int cfg::ReadCfg(std::string dir, int L, int &sym, int &ncf, std::vector<int> &cfgs) {
  std::ifstream cfgfile(dir + "/cfg-" + std::to_string(L) + ".inp");
  std::string line;
  std::vector<int> vals;
  int val;

  do {
    std::getline(cfgfile, line);
    std::istringstream iss(line);
    vals.clear();
    std::copy(std::istream_iterator<int>(iss),
              std::istream_iterator<int>(),
              std::back_inserter(vals));
  } while(vals.size()==1);

  sym = 2*vals[1]+1;

  std::getline(cfgfile, line);
  std::istringstream iss(line);
  iss >> ncf;

  cfgs.reserve(ncf*5);

  for(auto i=0; i<ncf; ++i) {
    std::getline(cfgfile, line);
    std::istringstream iss(line);
    for(auto j=0; j<5; ++j) { // read n1, l1, l2, n2_min, n2_max
      if(iss >> val)
        cfgs[i*5+j] = val;
      else {
        --i; 
        break;
      }
    }
  }
  return 0;
}
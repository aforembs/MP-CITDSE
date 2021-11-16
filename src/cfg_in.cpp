#include "cfg_in.h"

int cfg::ReadCfg(std::string dir, int L, int &sym, int &ncf, std::vector<cfg::line> &cfgs) {
  std::string filename = dir + "/cfg-" + std::to_string(L) + ".inp"; 
  std::ifstream cfgfile(filename);
  std::string line;
  std::vector<int> vals;
  int val;

  if(!std::filesystem::exists(filename)) {
    std::cout << "Input file: " << filename << " does not exist!\n";
    return -1;
  }

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

  cfg::line cf_l;
  cfgs.reserve(ncf);

  for(auto i=0; i<ncf; ++i) {
    std::getline(cfgfile, line);
    std::istringstream iss(line);
    if(iss >> val) {
      cf_l.n1 = val-1; // -1 for C indexing
      iss >> val;
      cf_l.l1 = val;
      iss >> val;
      cf_l.l2 = val;
      iss >> val;
      cf_l.n2min = val-1;
      iss >> val;
      cf_l.n2max = val;
      cfgs.emplace_back(cf_l);
    } else {
      --i; // ignore empty lines
    }
  }
  return 0;
}
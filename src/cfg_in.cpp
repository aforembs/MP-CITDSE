#include "cfg_in.hpp"

int cfg::readCfg(std::string dir, int L, int &sym, int &ncf,
                 std::vector<cfg::line> &cfgs) {
  std::string filename = dir + "/cfg-" + std::to_string(L) + ".inp";
  std::ifstream cfgfile(filename);
  std::string line;
  std::vector<int> vals;
  int val;

  if (!std::filesystem::exists(filename)) {
    std::cout << "Input file: " << filename << " does not exist!\n";
    return -1;
  }

  do {
    std::getline(cfgfile, line);
    std::istringstream iss(line);
    vals.clear();
    std::copy(std::istream_iterator<int>(iss), std::istream_iterator<int>(),
              std::back_inserter(vals));
  } while (vals.size() == 1);

  sym = 2 * vals[1] + 1;

  std::getline(cfgfile, line);
  std::istringstream iss(line);
  iss >> ncf;

  cfg::line cf_l;
  cfgs.resize(ncf);

  for (auto i = 0; i < ncf; ++i) {
    std::getline(cfgfile, line);
    std::istringstream liss(line);
    if (liss >> val) {
      cf_l.n1 = val - 1; // -1 for C indexing
      liss >> val;
      cf_l.l1 = val;
      liss >> val;
      cf_l.l2 = val;
      liss >> val;
      cf_l.n2min = val - 1;
      liss >> val;
      cf_l.n2max = val;
      if (abs(cf_l.l1 - cf_l.l2) <= L && L <= (cf_l.l1 + cf_l.l2) &&
          ((L >> 0) & 1) == (((cf_l.l1 + cf_l.l2) >> 0) & 1))
        cfgs[i] = cf_l;
      else
        std::cout << "cfg-" << L << ".inp, line: (" << i + 1
                  << ") has an invalid configuration, skipping\n";
    } else {
      --i; // ignore empty lines
    }
  }
  return 0;
}
#include "gr2e.hpp"
#include "gr_read.hpp"
#include <cstdlib>
#include <fenv.h>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  // feenableexcept(FE_INVALID | FE_OVERFLOW);

  std::string opt_file;
  char gauge;
  int L0_sz;
  double dt;
  std::string pot;
  std::string file_prefix;

  for (;;) {
    switch (getopt(argc, argv, "hf:")) {
    case 'h':
      std::cout << "Program for propagating the 2e tdse\n"
                << "-f <path> yaml input file with the input settings\n";
      return -1;
    case 'f':
      opt_file = optarg;
      continue;
    }
    break;
  }

  grrd::readConfig(opt_file, pot, gauge, L0_sz, dt);

  file_prefix = "dat/" + pot;

  std::vector<double> ens, block;

  grrd::readStructure(file_prefix, L0_sz, ens, block);

  gr2e::prop(file_prefix, L0_sz, 0.0, 0.01, block);
  // gr2e::eig(file_prefix, L0_sz, block);

  return 0;
}
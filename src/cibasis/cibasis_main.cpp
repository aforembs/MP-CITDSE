#include "cibasis.hpp"
#include <cstdlib>
#include <fenv.h>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file;
  std::string pot;
  std::string file_prefix;
  int L_max;
  char gauge;

  for (;;) {
    switch (getopt(argc, argv, "hf:")) {
    case 'h':
      std::cout << "Program for forming the CI correlated 2e^- basis\n"
                << "-f <path> yaml input file with the input settings\n";
      return -1;
    case 'f':
      opt_file = optarg;
      continue;
    }
    break;
  }

  cib::readConfig(opt_file, pot, gauge, L_max);
  file_prefix = "dat/" + pot;

  stvupt vecs;
  for (auto i = 0; i <= L_max; ++i) {
    vecs.push_back(
        std::make_unique<std::vector<double>>(std::vector<double>()));
  }

  cib::formCIh0(file_prefix, L_max, vecs);

  cib::formCIDipoles(file_prefix, gauge, L_max, vecs);

  return 0;
}
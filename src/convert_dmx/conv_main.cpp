#include "conv.hpp"
#include <cstdlib>
#include <fenv.h>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file, pot, file_prefix;
  int L_max;
  char gauge;
  std::vector<int> state_sz;

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

  conv::readConfig(opt_file, pot, L_max, gauge, state_sz);

  file_prefix = "dat/" + pot;

  stvupt vecs;
  stvupt dipoles;

  vecs.push_back(std::make_unique<std::vector<double>>(std::vector<double>()));
  for (auto i = 0; i < L_max; ++i) {
    vecs.push_back(
        std::make_unique<std::vector<double>>(std::vector<double>()));
    dipoles.push_back(
        std::make_unique<std::vector<double>>(std::vector<double>()));
  }

  conv::calcEvecs(file_prefix, L_max, state_sz, vecs);

  conv::readDipoles(file_prefix, gauge, L_max, state_sz, dipoles);

  conv::transDip(file_prefix, gauge, L_max, state_sz, vecs, dipoles);

  return 0;
}
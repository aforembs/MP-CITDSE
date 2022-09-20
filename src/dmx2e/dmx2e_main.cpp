#include "dmx2e.hpp"
#include <cstdlib>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file;
  int L_max, l_max;
  std::string pot;
  std::string out_prefix;
  char gauge;

  for (;;) {
    switch (getopt(argc, argv, "hf:")) {
    case 'h':
      std::cout << "Program for calculating the inter-electronic interaction\n"
                << "coefficients <n1l1;n2l2|r_12|n'1l'1;n'2l'2>\n"
                << "-f <path> yaml input file with the input settings\n";
      return -1;
    case 'f':
      opt_file = optarg;
      continue;
    }
    break;
  }

  dmx2e::ReadConfig(opt_file, pot, L_max, l_max, gauge);

  out_prefix = "dat/" + pot;

  dmx2e::SortL(out_prefix, L_max, gauge, "inp");

  dmx2e::GenDipole(out_prefix, L_max, l_max, gauge, "inp");

  return 0;
}
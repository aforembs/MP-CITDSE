
#include "dmx1e.hpp"
#include <cstdlib>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file;
  int glq_pt, l_max;
  std::string pot, integrator;
  std::string out_prefix;
  char gauge;

  for (;;) {
    switch (getopt(argc, argv, "hf:")) {
    case 'h':
      std::cout << "Program for calculating the 1e dipole matrix\n"
                << "elements <nl|d_(v/l)|n'l'>\n"
                << "-f <path> yaml input file with the input settings\n";
      return -1;
    case 'f':
      opt_file = optarg;
      continue;
    }
    break;
  }

  dmx1e::ReadConfig(opt_file, pot, glq_pt, gauge, l_max);

  out_prefix = "dat/" + pot;

  dmx1e::GenDipole(out_prefix, glq_pt, gauge, l_max);

  return 0;
}
#include "dmx2e.hpp"
#include <cstdlib>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file;
  std::string inp_dir;
  int L_max, l_max;
  std::string pot;
  std::string out_prefix;
  char gauge;

  for (;;) {
    switch (getopt(argc, argv, "hf:i:")) {
    case 'h':
      std::cout << "Program for calculating the 2e^- dipole matrices\n"
                << "-f <path> yaml input file with the input settings\n"
                << "-i <path> directory containing cfg-<L>.inp files\n";
      return -1;
    case 'f':
      opt_file = optarg;
      continue;
    case 'i':
      inp_dir = optarg;
      continue;
    }
    break;
  }

  dmx2e::readConfig(opt_file, pot, L_max, l_max, gauge);

  out_prefix = "dat/" + pot;

  dmx2e::genDipole(out_prefix, L_max, l_max, gauge, inp_dir);

  return 0;
}
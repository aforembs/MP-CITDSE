
#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include "dmx1e.h"

int main(int argc, char *argv[]) {
  std::string opt_file;
  int glq_pt, l_max;
  std::string pot, integrator;
  std::string out_prefix;
  char gauge;

  for(;;) {
    switch(getopt(argc, argv, "hf:")) {
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

  dmx1e::ReadConfig(opt_file, glq_pt, pot, l_max, gauge);

  out_prefix = "dat/" + pot;

  dmx1e::GenDipole(out_prefix, l_max, glq_pt, gauge);

  return 0;
}
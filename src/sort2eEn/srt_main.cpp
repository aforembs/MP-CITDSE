#include "srt2e.hpp"
#include <cstdlib>
#include <filesystem>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file, pot, out_prefix;
  int L_max;

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

  srt2e::readConfig(opt_file, pot, L_max);

  out_prefix = "dat/" + pot;

  srt2e::sortEn(out_prefix, L_max, "inp");

  return 0;
}
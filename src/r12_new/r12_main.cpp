#include "r12_new.hpp"
#include <cstdlib>
#include <filesystem>
#include <unistd.h>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
  std::string opt_file;
  int qsz, L_max;
  std::string pot, integrator;
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

  r_12::readConfig(opt_file, qsz, pot, L_max, gauge, integrator);

  out_prefix = "dat/" + pot;

  r_12::r12Glob(out_prefix, L_max, qsz, "inp");

  return 0;
}
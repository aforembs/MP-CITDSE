#include "w1e.hpp"
#include <cstdlib>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file;
  int qsz, R_max, l_max;
  std::string pot;
  std::string out_prefix;

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

  w1e::ReadConfig(opt_file, qsz, R_max, l_max, pot);

  out_prefix = "dat/" + pot;

  w1e::GenWfn(out_prefix, qsz, R_max, l_max);

  return 0;
}
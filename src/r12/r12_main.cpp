#include "r12.hpp"
#include <cstdlib>
#include <filesystem>
#include <unistd.h>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
  std::string opt_file;
  std::string inp_dir;
  int qsz, L_max;
  std::string pot, integrator;
  std::string out_prefix;
  bool lim_flag = false;
  std::string k_limit;

  for (;;) {
    switch (getopt(argc, argv, "hf:i:")) {
    case 'h':
      std::cout << "Program for calculating the inter-electronic interaction\n"
                << "coefficients <n1l1;n2l2|r_12|n'1l'1;n'2l'2>\n"
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

  r_12::readConfig(opt_file, qsz, pot, L_max, k_limit, lim_flag);

  out_prefix = "dat/" + pot;

  r_12::r12Glob(out_prefix, L_max, qsz, inp_dir, lim_flag, k_limit);

  return 0;
}
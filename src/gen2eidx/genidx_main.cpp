#include "genidx.hpp"
#include <cstdlib>
#include <filesystem>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file, inp_dir, pot, out_prefix;
  int L_max;

  for (;;) {
    switch (getopt(argc, argv, "hf:i:")) {
    case 'h':
      std::cout << "Program for sorting the direct product 2e^- states in \n"
                << "terms of ascending total energy\n"
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

  genidx::readConfig(opt_file, pot, L_max);

  auto cfg_prefix = inp_dir + "/cfg-";
  for (int L = 0; L <= L_max; ++L) {
    std::filesystem::copy_file(cfg_prefix + std::to_string(L) + ".inp",
                               "dat/cfg-" + std::to_string(L) + ".inp",
                               std::filesystem::copy_options::update_existing);
  }

  out_prefix = "dat/" + pot;

  genidx::saveIdx(out_prefix, L_max, inp_dir);

  return 0;
}
#include "r12.hpp"
#include <cstdlib>
#include <filesystem>
#include <unistd.h>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
  std::string opt_file;
  int glq_pt, L_max;
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

  r_12::readConfig(opt_file, glq_pt, pot, L_max, gauge, integrator);

  out_prefix = "dat/" + pot;

  // check if he<L_max>idx.h5 exists if not create index files
  if (!fs::exists(out_prefix + std::to_string(L_max) + "idx.h5")) {
    cfg::genL_idx(out_prefix, gauge, L_max, "inp");
  }

  if (integrator.compare("mixed") == 0)
    r_12::r12MM(out_prefix, L_max, glq_pt, "inp");
  else if (integrator.compare("glob4") == 0)
    r_12::r12Glob4(out_prefix, L_max, glq_pt, "inp");
  else if (integrator.compare("glob3") == 0)
    r_12::r12Glob3(out_prefix, L_max, glq_pt, "inp");
  else if (integrator.compare("trapezoid") == 0)
    r_12::r12Trap(out_prefix, L_max, glq_pt, "inp");
  else
    std::cout << "Invalid integration scheme\n";

  return 0;
}
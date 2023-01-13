#include "pes.hpp"
#include <cstdlib>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file, in_file = "dat/ct.dat";
  int L_max, l_max;
  std::string pot;
  std::string input_prefix, output_prefix;
  std::vector<int> state_sz;

  for (;;) {
    switch (getopt(argc, argv, "hf:i:")) {
    case 'h':
      std::cout << "Program for generating the PES from coeffs\n"
                << "-f <path> yaml input file with the input settings\n"
                << "-i <path> file containing the coefficients\n";
      return -1;
    case 'f':
      opt_file = optarg;
      continue;
    case 'i':
      in_file = optarg;
      continue;
    }
    break;
  }

  pes::readConfig(opt_file, pot, L_max, l_max, state_sz);

  input_prefix = "dat/" + pot;

  std::vector<std::complex<double>> ct;

  pes::readCt(in_file, ct);

  std::filesystem::path inpath(in_file);
  output_prefix = inpath.parent_path().u8string() + "/" + pot;

  pes::genPES(input_prefix, L_max, state_sz, ct, output_prefix);

  return 0;
}
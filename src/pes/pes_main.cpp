#include "pes.hpp"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file, in_file = "dat/ct.dat";
  int L_max, e_num = 0;
  std::string pot;
  std::string base, option;
  std::string input_prefix, output_prefix;
  std::vector<int> state_sz;
  bool s_flag = 0;

  for (;;) {
    switch (getopt(argc, argv, "he:f:i:sl")) {
    case 'h':
      std::cout << "Program for generating the PES from coeffs\n"
                << "-e <1 or 2> the number of electrons in the tdse\n"
                << "-f <path> yaml input file with the input settings\n"
                << "-i <path> file containing the coefficients\n"
                << "-s generate PES summed over all 'l/L'\n";
      return -1;
    case 'e':
      e_num = std::stoi(optarg);
      if (e_num != 1 && e_num != 2) {
        std::cout << "Invalid number of electrons, use 1 or 2\n";
        return -1;
      }
      continue;
    case 'f':
      opt_file = optarg;
      continue;
    case 'i':
      in_file = optarg;
      continue;
    case 's':
      s_flag = 1;
      continue;
    }
    break;
  }

  std::filesystem::path inpath(in_file);
  if (!std::filesystem::exists(inpath) ||
      !std::filesystem::is_regular_file(in_file)) {
    std::cout << "Invalid coefficient file path, -i <path>!!\n";
    return -1;
  }

  switch (e_num) {
  case 1:
    base = "Basis_Settings";
    option = "l_max";
    break;
  case 2:
    base = "Global_Settings";
    option = "L_max";
    break;
  }

  pes::readConfig(opt_file, pot, base, option, L_max, state_sz);

  input_prefix = "dat/" + pot;

  std::vector<std::complex<double>> ct;

  pes::readCt(in_file, ct);

  output_prefix = inpath.parent_path().u8string() + "/" + pot;

  switch (e_num) {
  case 1:
    pes::genPES1e(input_prefix, s_flag, L_max, state_sz, ct, output_prefix);
    break;
  case 2:
    pes::genPES2e(input_prefix, s_flag, L_max, state_sz, ct, output_prefix);
    break;
  }

  return 0;
}
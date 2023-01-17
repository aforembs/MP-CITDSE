#include "td_read.hpp"
#include "tdse.hpp"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file;
  std::string out_dir;
  int l_max, cycles;
  std::string pot;
  std::string i_file_prefix;
  std::string o_file_prefix;
  char gauge;
  double dt, w, Io, cepd;
  std::vector<int> state_sz;

  for (;;) {
    switch (getopt(argc, argv, "hf:o:")) {
    case 'h':
      std::cout << "Program for propagating the 1e tdse\n"
                << "-f <path> yaml input file with the input settings\n"
                << "-o <path> output folder\n";
      return -1;
    case 'f':
      opt_file = optarg;
      continue;
    case 'o':
      out_dir = optarg;
      continue;
    }
    break;
  }

  const std::filesystem::path out_path{out_dir};
  if (!std::filesystem::exists(out_path)) {
    std::filesystem::create_directory(out_path);
  }

  tdrd::readConfig(opt_file, pot, gauge, l_max, state_sz, dt, w, Io, cepd,
                   cycles);

  i_file_prefix = "dat/" + pot;
  o_file_prefix = out_dir + "/" + pot;

  int e_sz = 0;
  std::vector<double> eps;
  std::vector<int> offs;
  std::vector<double *> eig;
  stvupt dipoles;

  for (auto i = 0; i < l_max; ++i) {
    dipoles.push_back(
        std::make_unique<std::vector<double>>(std::vector<double>()));
  }

  tdrd::readEnergies(i_file_prefix, l_max, state_sz, eps, offs, e_sz, eig);

  tdrd::readDipoles(i_file_prefix, gauge, l_max, state_sz, dipoles);

  std::vector<std::complex<double>> ct(e_sz);
  ct[0] = std::complex<double>(1.0, 0.0);

  double t = 0.0;
  double tau = pulse::sineT(w * conv::En_ev_au_, cycles);
  int steps = tau / dt;

  switch (gauge) {
  case 'v':
    tdse::propV(o_file_prefix, l_max, t, dt, steps, pulse::sineASetup,
                pulse::sineAA, Io, w, cepd, cycles, e_sz, offs, state_sz, eig,
                dipoles, ct);
    break;
  case 'l':
    tdse::propL(o_file_prefix, l_max, t, dt, steps, pulse::sineESetup,
                pulse::sineEE, Io, w, cepd, cycles, e_sz, offs, state_sz, eig,
                dipoles, ct);
    break;
  }

  return 0;
}
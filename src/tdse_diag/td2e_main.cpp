#include "td2e.hpp"
#include "td_read.hpp"
#include <cstdlib>
#include <fenv.h>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  // feenableexcept(FE_INVALID | FE_OVERFLOW);

  std::string opt_file;
  std::string out_dir;
  int L_max, Lanc_iter, num_eval, cycles;
  std::string pot;
  std::string i_file_prefix;
  std::string o_file_prefix;
  char gauge;
  double dt, w, Io, cepd;
  std::vector<int> state_sz;

  for (;;) {
    switch (getopt(argc, argv, "hf:o:")) {
    case 'h':
      std::cout << "Program for propagating the 2e tdse\n"
                << "-f <path> yaml input file with the input settings\n";
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

  tdrd::readConfig(opt_file, pot, gauge, L_max, state_sz, Lanc_iter, num_eval,
                   dt, w, Io, cepd, cycles);

  i_file_prefix = "dat/" + pot;
  o_file_prefix = out_dir + "/" + pot;

  int ct_sz = 0;
  std::vector<int> offs;
  stvupt eig;
  stvupt blocks;
  stvupt dipoles;

  eig.push_back(std::make_unique<std::vector<double>>(std::vector<double>()));
  blocks.push_back(
      std::make_unique<std::vector<double>>(std::vector<double>()));
  for (auto i = 0; i < L_max; ++i) {
    eig.push_back(std::make_unique<std::vector<double>>(std::vector<double>()));
    blocks.push_back(
        std::make_unique<std::vector<double>>(std::vector<double>()));
    dipoles.push_back(
        std::make_unique<std::vector<double>>(std::vector<double>()));
  }

  tdrd::readStructure(i_file_prefix, L_max, ct_sz, state_sz, offs, eig, blocks);

  tdrd::readDipoles(i_file_prefix, gauge, L_max, state_sz, blocks, dipoles);

  std::vector<std::complex<double>> ct(ct_sz);

  ct[0] = std::complex<double>(1.0, 0.0);

  double t = 0.0;
  double tau = pulse::Sine_T(w * conv::En_ev_au_, cycles);
  int steps = tau / dt;

  switch (gauge) {
  case 'v':
    td2e::propV(o_file_prefix, L_max, t, dt, steps, pulse::SineA_Setup,
                pulse::SineA_A, w, Io, cepd, cycles, ct_sz, offs, state_sz, eig,
                dipoles, ct);
    break;
  case 'l':
    td2e::propL(o_file_prefix, L_max, t, dt, steps, pulse::SineA_Setup,
                pulse::SineA_E, w, Io, cepd, cycles, ct_sz, offs, state_sz, eig,
                dipoles, ct);
    break;
  }

  return 0;
}
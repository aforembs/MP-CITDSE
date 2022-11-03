#include "td2e.hpp"
#include "td_read.hpp"
#include <cstdlib>
#include <fenv.h>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  // feenableexcept(FE_INVALID | FE_OVERFLOW);

  std::string opt_file;
  int L_max, Lanc_iter, num_eval, cycles;
  std::string pot;
  std::string file_prefix;
  char gauge;
  double dt, w, Io, cepd;
  std::vector<int> state_sz;

  for (;;) {
    switch (getopt(argc, argv, "hf:")) {
    case 'h':
      std::cout << "Program for propagating the 2e tdse\n"
                << "-f <path> yaml input file with the input settings\n";
      return -1;
    case 'f':
      opt_file = optarg;
      continue;
    }
    break;
  }

  tdrd::readConfig(opt_file, pot, gauge, L_max, state_sz, Lanc_iter, num_eval,
                   dt, w, Io, cepd, cycles);

  file_prefix = "dat/" + pot;

  int ct_sz = 0;
  std::vector<int> offs;
  stvupt blocks;
  stvupt dipoles;

  blocks.push_back(
      std::make_unique<std::vector<double>>(std::vector<double>()));
  for (auto i = 0; i < L_max; ++i) {
    blocks.push_back(
        std::make_unique<std::vector<double>>(std::vector<double>()));
    dipoles.push_back(
        std::make_unique<std::vector<double>>(std::vector<double>()));
  }

  tdrd::readStructure(file_prefix, L_max, ct_sz, state_sz, offs, blocks);

  tdrd::readDipoles(file_prefix, gauge, L_max, state_sz, dipoles);

  std::vector<std::complex<double>> ct(ct_sz);

  tdrd::readGrCt(file_prefix, state_sz, ct);

  double t = 0.0;
  double tau = pulse::Sine_T(w * conv::En_ev_au_, cycles);
  int steps = tau / dt;
  steps = steps + steps / 2;

  switch (gauge) {
  case 'v':
    td2e::propV(file_prefix, L_max, t, dt, steps, pulse::SineA_Setup,
                pulse::SineA_A, w, Io, cepd, cycles, ct_sz, offs, state_sz,
                blocks, dipoles, ct);
    break;
  case 'l':
    td2e::propL(file_prefix, L_max, t, dt, steps, pulse::SineA_Setup,
                pulse::SineA_E, w, Io, cepd, cycles, ct_sz, offs, state_sz,
                blocks, dipoles, ct);
    break;
  }

  return 0;
}
#include "td_read.hpp"
#include "tdse.hpp"
#include <cstdlib>
#include <fenv.h>
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

  tdrd::readConfig(opt_file, pot, gauge, l_max, state_sz, dt, w, Io, cepd,
                   cycles);

  i_file_prefix = "dat/" + pot;
  o_file_prefix = out_dir + "/" + pot;

  int e_sz = 0;
  std::vector<double> eps, dip;
  std::vector<int> offs, doffs;
  std::vector<double *> eps_off;
  std::vector<double *> dip_off;
  tdrd::readEnergies(i_file_prefix, l_max, state_sz, eps, offs, e_sz, eps_off);

  tdrd::readDipoles(i_file_prefix, gauge, l_max, state_sz, dip, doffs, dip_off);

  std::vector<std::complex<double>> ct(e_sz);
  ct[0] = std::complex<double>(1.0, 0.0);

  double t = 0.0;
  double tau = pulse::Sine_T(w * conv::En_ev_au_, cycles);
  int steps = tau / dt;

  switch (gauge) {
  case 'v':
    tdse::propV(o_file_prefix, l_max, t, dt, steps, pulse::SineA_Setup,
                pulse::SineA_A, Io, w, cepd, cycles, e_sz, offs, doffs,
                state_sz, eps_off, dip_off, ct);
    break;
  case 'l':
    tdse::propL(o_file_prefix, l_max, t, dt, steps, pulse::SineA_Setup,
                pulse::SineA_A, Io, w, cepd, cycles, e_sz, offs, doffs,
                state_sz, eps_off, dip_off, ct);
    break;
  }

  return 0;
}
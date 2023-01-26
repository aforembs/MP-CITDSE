#include "td_read.hpp"
#include "tdse.hpp"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file;
  std::string out_dir;
  int L_max, cycles, e_num = 0;
  std::string pot;
  std::string base, option;
  std::string en_pot, dip_pot;
  std::string en_set, dip_set;
  std::string o_file_prefix;
  char gauge;
  double dt, w, Io, cepd;
  std::vector<int> state_sz;

  for (;;) {
    switch (getopt(argc, argv, "he:f:o:")) {
    case 'h':
      std::cout << "Program for propagating either the 1- or 2-e^- tdse\n"
                << "-e <1 or 2> the number of electrons in the tdse\n"
                << "-f <path> yaml input file with the input settings\n"
                << "-o <path> output directory for the coefficient files\n";
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

  tdrd::readConfig(opt_file, pot, gauge, base, option, L_max, state_sz, dt, w,
                   Io, cepd, cycles);

  switch (e_num) {
  case 1:
    en_pot = "dat/" + pot;
    en_set = "En";
    dip_pot = "dat/" + pot;
    dip_set = "d_if";
    break;
  case 2:
    en_pot = "dat/" + pot + "CI";
    en_set = "En_CI";
    dip_pot = "dat/" + pot + "CI_";
    dip_set = "CI_dip";
    break;
  }

  o_file_prefix = out_dir + "/" + pot;

  int ct_sz = 0;
  std::vector<int> offs;
  stvupt eig;
  stvupt dipoles;

  eig.push_back(std::make_unique<std::vector<double>>(std::vector<double>()));
  for (auto i = 0; i < L_max; ++i) {
    eig.push_back(std::make_unique<std::vector<double>>(std::vector<double>()));
    dipoles.push_back(
        std::make_unique<std::vector<double>>(std::vector<double>()));
  }

  tdrd::readEnergies(en_pot, en_set, L_max, ct_sz, state_sz, offs, eig);

  tdrd::readDipoles(dip_pot, dip_set, gauge, L_max, state_sz, dipoles);

  std::vector<std::complex<double>> ct(ct_sz);
  ct[0] = std::complex<double>(1.0, 0.0);

  double t = 0.0;
  double tau = pulse::sineT(w * conv::En_ev_au_, cycles);
  int steps = tau / dt;

  switch (gauge) {
  case 'v':
    tdse::propV(o_file_prefix, L_max, t, dt, steps, pulse::sineASetup,
                pulse::sineAA, w, Io, cepd, cycles, ct_sz, offs, state_sz, eig,
                dipoles, ct);
    break;
  case 'l':
    tdse::propL(o_file_prefix, L_max, t, dt, steps, pulse::sineESetup,
                pulse::sineEE, w, Io, cepd, cycles, ct_sz, offs, state_sz, eig,
                dipoles, ct);
    break;
  }

  return 0;
}
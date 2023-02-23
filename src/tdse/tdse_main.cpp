#include "td_read.hpp"
#include "tdse.hpp"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <unistd.h>

int main(int argc, char *argv[]) {
  std::string opt_file;
  std::string out_dir;
  int L_max, cycles, e_num = 0;
  int pop_n = 0, pop_l = 0;
  std::string pot, init_ct = "ground";
  std::string base, option;
  std::string en_pot, dip_pot;
  std::string en_set, dip_set;
  std::string o_file_prefix;
  char gauge;
  std::string shape;
  double dt, w, Io, cepd;
  std::vector<int> state_sz;

  for (;;) {
    switch (getopt(argc, argv, "he:f:o:s:n:l:")) {
    case 'h':
      std::cout << "Program for propagating either the 1- or 2-e^- tdse\n"
                << "-e <1 or 2> the number of electrons in the tdse\n"
                << "-f <path> yaml input file with the input settings\n"
                << "-o <path> output directory for the coefficient files\n"
                << "-s <path> file with the start coefficients\n"
                << "-n <int> primary quantum number of state to track\n"
                << "-l <int> angular quantum number of state to track\n";
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
    case 's':
      init_ct = optarg;
      continue;
    case 'n':
      pop_n = std::stoi(optarg);
      if (pop_n < 1) {
        std::cout << "Invalid primary quantum number -n, use int > 0\n";
        return -1;
      }
      --pop_n;
      continue;
    case 'l':
      pop_l = std::stoi(optarg);
      if (pop_n < 0) {
        std::cout << "Invalid angular momentum number -l, use int > 0\n";
        return -1;
      }
      continue;
    }
    break;
  }

  const std::filesystem::path out_path{out_dir};
  if (!std::filesystem::exists(out_path)) {
    std::filesystem::create_directory(out_path);
  }

  auto set_file = out_dir + "/settings.yaml";
  const std::filesystem::path file_path{set_file};
  std::filesystem::copy_options cpopt = std::filesystem::copy_options::none;
  if (std::filesystem::exists(file_path)) {
    char ovw = 'n';
    std::cout << "File " << set_file << " exists. Overwrite? [y/n]\n";
    std::cin >> ovw;
    switch (ovw) {
    case 'y':
      cpopt = std::filesystem::copy_options::update_existing;
      break;
    case 'n':
      std::cout << out_dir
                << " contians old data, please choose a different output"
                << " directory!\n";
      return -1;
    }
  }

  std::filesystem::copy_file(opt_file, set_file, cpopt);

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

  tdrd::readConfig(opt_file, pot, base, option, L_max, gauge, state_sz, dt,
                   shape, w, Io, cepd, cycles);

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
    for (int L = 0; L <= L_max; ++L) {
      std::filesystem::copy_file("dat/cfg-" + std::to_string(L) + ".inp",
                                 out_dir + "/cfg-" + std::to_string(L) + ".inp",
                                 cpopt);
    }
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

  if (init_ct.compare("ground") != 0) {
    auto t_start = init_ct.find_first_of("0123456789");
    auto t_end = init_ct.find_last_of("123456789");
    auto t_sz = t_end - t_start + 1;

    t = std::stod(&init_ct[t_start], &t_sz);
    std::cout << "Start time (a.u.): " << t << "\n";
    tdrd::readInitCt(init_ct, ct_sz, ct);
  }

  double tau = pulse::period(w * conv::En_ev_au_, cycles);
  int steps = (tau - t) / dt;

  double auIo, auw;
  pulse::toAU(Io, w, auIo, auw);
  pulse::params pars;

  if (shape.compare("sine") == 0) {
    switch (gauge) {
    case 'v':
      pulse::sineASetup(auIo, auw, 0.0, cycles, pars);
      tdse::propV(o_file_prefix, L_max, t, dt, steps, pop_n, pop_l,
                  pulse::sineAA, pars, ct_sz, offs, state_sz, eig, dipoles, ct);
      break;
    case 'l':
      pulse::sineESetup(auIo, auw, 0.0, cycles, pars);
      tdse::propL(o_file_prefix, L_max, t, dt, steps, pop_n, pop_l,
                  pulse::sineEE, pars, ct_sz, offs, state_sz, eig, dipoles, ct);
      break;
    }
  } else if (shape.compare("gaussian") == 0) {
    switch (gauge) {
    case 'v':
      pulse::gaussASetup(auIo, auw, tau, cycles, 0.0, pars);
      tdse::propV(o_file_prefix, L_max, t, dt, steps, pop_n, pop_l,
                  pulse::gauss, pars, ct_sz, offs, state_sz, eig, dipoles, ct);
      break;
    case 'l':
      pulse::gaussESetup(auIo, auw, tau, cycles, 0.0, pars);
      tdse::propL(o_file_prefix, L_max, t, dt, steps, pop_n, pop_l,
                  pulse::gauss, pars, ct_sz, offs, state_sz, eig, dipoles, ct);
      break;
    }
  } else {
    std::cout << "Invalid envelope shape, use \"sine\" or \"gaussian\"\n";
    return -1;
  }

  return 0;
}

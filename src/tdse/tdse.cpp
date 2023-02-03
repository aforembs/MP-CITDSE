#include "tdse.hpp"

typedef std::vector<std::complex<double>> state_type;

class MatVecV {
public:
  int L_max;
  double field;
  std::vector<int> state_sz;
  std::vector<int> offs;
  std::vector<double *> eig;
  std::vector<std::complex<double> *> cdipole;

  void operator()(state_type &x, state_type &dxdt,
                  [[maybe_unused]] double t) const {
    constexpr std::complex<double> mI(0.0, -1.0);
    auto alp_fl = std::complex<double>(field, 0.0);
    auto alp_flm = std::complex<double>(-field, 0.0);
    auto beta = std::complex<double>(0.0, 0.0);
    auto bt2 = std::complex<double>(1.0, 0.0);

    for (auto i = 0; i < state_sz[0]; ++i) {
      dxdt[i] = eig[0][i] * mI * x[i];
    }

    cblas_zgemv(CblasRowMajor, CblasTrans, state_sz[1], state_sz[0],
                reinterpret_cast<double *>(&alp_flm),
                reinterpret_cast<double *>(cdipole[0]), state_sz[0],
                reinterpret_cast<double *>(&x[offs[1]]), 1,
                reinterpret_cast<double *>(&bt2),
                reinterpret_cast<double *>(&dxdt[0]), 1);

    for (auto L = 1; L < L_max; ++L) {
      cblas_zgemv(CblasRowMajor, CblasNoTrans, state_sz[L], state_sz[L - 1],
                  reinterpret_cast<double *>(&alp_fl),
                  reinterpret_cast<double *>(cdipole[L - 1]), state_sz[L - 1],
                  reinterpret_cast<double *>(&x[offs[L - 1]]), 1,
                  reinterpret_cast<double *>(&beta),
                  reinterpret_cast<double *>(&dxdt[offs[L]]), 1);

      for (auto i = 0; i < state_sz[L]; ++i) {
        dxdt[offs[L] + i] =
            eig[L][i] * mI * x[offs[L] + i] + bt2 * dxdt[offs[L] + i];
      }

      cblas_zgemv(CblasRowMajor, CblasTrans, state_sz[L + 1], state_sz[L],
                  reinterpret_cast<double *>(&alp_flm),
                  reinterpret_cast<double *>(cdipole[L]), state_sz[L],
                  reinterpret_cast<double *>(&x[offs[L + 1]]), 1,
                  reinterpret_cast<double *>(&bt2),
                  reinterpret_cast<double *>(&dxdt[offs[L]]), 1);
    }

    cblas_zgemv(CblasRowMajor, CblasNoTrans, state_sz[L_max],
                state_sz[L_max - 1], reinterpret_cast<double *>(&alp_fl),
                reinterpret_cast<double *>(cdipole[L_max - 1]),
                state_sz[L_max - 1],
                reinterpret_cast<double *>(&x[offs[L_max - 1]]), 1,
                reinterpret_cast<double *>(&beta),
                reinterpret_cast<double *>(&dxdt[offs[L_max]]), 1);

    for (auto i = 0; i < state_sz[L_max]; ++i) {
      dxdt[offs[L_max] + i] =
          eig[L_max][i] * mI * x[offs[L_max] + i] + bt2 * dxdt[offs[L_max] + i];
    }
  }
};

class MatVecL {
public:
  int L_max;
  double field;
  std::vector<int> state_sz;
  std::vector<int> offs;
  std::vector<double *> eig;
  std::vector<std::complex<double> *> cdipole;

  void operator()(state_type &x, state_type &dxdt,
                  [[maybe_unused]] double t) const {
    constexpr std::complex<double> mI(0.0, -1.0);
    auto alp_fl = std::complex<double>(0.0, -field);
    auto beta = std::complex<double>(0.0, 0.0);
    auto bt2 = std::complex<double>(1.0, 0.0);

    for (auto i = 0; i < state_sz[0]; ++i) {
      dxdt[i] = eig[0][i] * mI * x[i];
    }

    cblas_zgemv(CblasRowMajor, CblasTrans, state_sz[1], state_sz[0],
                reinterpret_cast<double *>(&alp_fl),
                reinterpret_cast<double *>(cdipole[0]), state_sz[0],
                reinterpret_cast<double *>(&x[offs[1]]), 1,
                reinterpret_cast<double *>(&bt2),
                reinterpret_cast<double *>(&dxdt[0]), 1);

    for (auto L = 1; L < L_max; ++L) {
      cblas_zgemv(CblasRowMajor, CblasNoTrans, state_sz[L], state_sz[L - 1],
                  reinterpret_cast<double *>(&alp_fl),
                  reinterpret_cast<double *>(cdipole[L - 1]), state_sz[L - 1],
                  reinterpret_cast<double *>(&x[offs[L - 1]]), 1,
                  reinterpret_cast<double *>(&beta),
                  reinterpret_cast<double *>(&dxdt[offs[L]]), 1);

      for (auto i = 0; i < state_sz[L]; ++i) {
        dxdt[offs[L] + i] =
            eig[L][i] * mI * x[offs[L] + i] + bt2 * dxdt[offs[L] + i];
      }

      cblas_zgemv(CblasRowMajor, CblasTrans, state_sz[L + 1], state_sz[L],
                  reinterpret_cast<double *>(&alp_fl),
                  reinterpret_cast<double *>(cdipole[L]), state_sz[L],
                  reinterpret_cast<double *>(&x[offs[L + 1]]), 1,
                  reinterpret_cast<double *>(&bt2),
                  reinterpret_cast<double *>(&dxdt[offs[L]]), 1);
    }

    cblas_zgemv(CblasRowMajor, CblasNoTrans, state_sz[L_max],
                state_sz[L_max - 1], reinterpret_cast<double *>(&alp_fl),
                reinterpret_cast<double *>(cdipole[L_max - 1]),
                state_sz[L_max - 1],
                reinterpret_cast<double *>(&x[offs[L_max - 1]]), 1,
                reinterpret_cast<double *>(&beta),
                reinterpret_cast<double *>(&dxdt[offs[L_max]]), 1);

    for (auto i = 0; i < state_sz[L_max]; ++i) {
      dxdt[offs[L_max] + i] =
          eig[L_max][i] * mI * x[offs[L_max] + i] + bt2 * dxdt[offs[L_max] + i];
    }
  }
};

int tdse::propV(std::string output, int L_max, double t, double dt, int steps,
                int pop_n, int pop_l, fieldInit fieldst, fieldFcn field,
                double w, double Io, double cepd, int cycles, int ct_sz,
                std::vector<int> &offs, std::vector<int> &state_sz, stvupt &eig,
                stvupt &dipoles, std::vector<std::complex<double>> &ct) {
  double IoA, wA, Ao, cepds, Wenv;
  pulse::toAU(Io, w, IoA, wA);

  int print = steps / 10;

  fieldst(IoA, wA, cepd, cycles, Ao, cepds, Wenv);

  std::cout << "Io (a.u.): " << IoA << " w (a.u.): " << wA << "\n";

  std::fstream f_out, field_fl, f_pop;
  f_out.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
  for (auto k = 0; k < ct_sz; ++k) {
    f_out << k << "  " << ct[k].real() << "  " << ct[k].imag() << " "
          << std::norm(ct[k]) << "\n";
  }
  f_out.close();

  std::vector<state_type> cdipole;
  for (auto L = 0; L < L_max; ++L) {
    auto L_sz = state_sz[L];
    auto L1_sz = state_sz[L + 1];
    cdipole.push_back(state_type(L_sz * L1_sz));

    for (auto i = 0; i < L_sz; ++i) {
      for (auto jd = 0; jd < L1_sz; ++jd) {
        cdipole[L][i * L1_sz + jd] =
            std::complex<double>(dipoles[L].get()->at(i * L1_sz + jd), 0.0);
      }
    }
  }

  MatVecV MV;
  MV.L_max = L_max;
  MV.state_sz = state_sz;
  MV.offs = offs;

  for (auto L = 0; L < L_max; ++L) {
    MV.eig.push_back(eig[L]->data());
    MV.cdipole.push_back(cdipole[L].data());
  }
  MV.eig.push_back(eig[L_max]->data());

  boost::numeric::odeint::runge_kutta_fehlberg78<state_type> rkf;

  field_fl.open(output + "_field.dat", std::ios::out);
  f_pop.open(output + "_pop" + std::to_string(pop_n + 1) +
                 std::to_string(pop_l) + ".dat",
             std::ios::out);

  f_pop << t << " " << std::norm(ct[offs[pop_l] + pop_n]) << "\n";

  for (auto st = 0; st < steps; ++st) {
    field_fl << t + dt << " " << field(Ao, wA, cepds, Wenv, t + dt) << "\n";
    MV.field = field(Ao, wA, cepds, Wenv, t + dt);

    rkf.do_step(MV, ct, t, dt);

    t += dt;

    auto ctnrm = cblas_dznrm2(ct_sz, reinterpret_cast<double *>(&ct[0]), 1);

    for (auto &n : ct) {
      n /= ctnrm;
    }

    f_pop << std::setprecision(16) << t << " "
          << std::norm(ct[offs[pop_l] + pop_n]) << " " << ctnrm << "\n";

    if (st % print == 0) {
      std::cout << field(Ao, wA, cepds, Wenv, t) << "\n";
      std::fstream f_ct;
      f_ct.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
      for (auto k = 0; k < ct_sz; ++k) {
        f_ct << k << "  " << ct[k].real() << "  " << ct[k].imag() << " "
             << std::norm(ct[k]) << "\n";
      }
      f_ct.close();
    }
  }
  field_fl.close();
  f_pop.close();
  return 0;
}

int tdse::propL(std::string output, int L_max, double t, double dt, int steps,
                int pop_n, int pop_l, fieldInit fieldst, fieldFcn field,
                double w, double Io, double cepd, int cycles, int ct_sz,
                std::vector<int> &offs, std::vector<int> &state_sz, stvupt &eig,
                stvupt &dipoles, std::vector<std::complex<double>> &ct) {
  double IoA, wA, Ao, cepds, Wenv;
  pulse::toAU(Io, w, IoA, wA);

  int print = steps / 10;

  fieldst(IoA, wA, cepd, cycles, Ao, cepds, Wenv);

  std::cout << "Io (a.u.): " << IoA << " w (a.u.): " << wA << "\n";

  std::fstream f_out, field_fl, f_pop;
  f_out.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
  for (auto k = 0; k < ct_sz; ++k) {
    f_out << k << "  " << ct[k].real() << "  " << ct[k].imag() << " "
          << std::norm(ct[k]) << "\n";
  }
  f_out.close();

  std::vector<state_type> cdipole;
  for (auto L = 0; L < L_max; ++L) {
    auto L_sz = state_sz[L];
    auto L1_sz = state_sz[L + 1];
    cdipole.push_back(state_type(L_sz * L1_sz));

    for (auto i = 0; i < L_sz; ++i) {
      for (auto jd = 0; jd < L1_sz; ++jd) {
        cdipole[L][i * L1_sz + jd] =
            std::complex<double>(dipoles[L].get()->at(i * L1_sz + jd), 0.0);
      }
    }
  }

  MatVecL MV;
  MV.L_max = L_max;
  MV.state_sz = state_sz;
  MV.offs = offs;

  for (auto L = 0; L < L_max; ++L) {
    MV.eig.push_back(eig[L]->data());
    MV.cdipole.push_back(cdipole[L].data());
  }
  MV.eig.push_back(eig[L_max]->data());

  boost::numeric::odeint::runge_kutta_fehlberg78<state_type> rkf;

  field_fl.open(output + "_field.dat", std::ios::out);
  f_pop.open(output + "_pop" + std::to_string(pop_n + 1) +
                 std::to_string(pop_l) + ".dat",
             std::ios::out);

  f_pop << t << " " << std::norm(ct[offs[pop_l] + pop_n]) << "\n";

  for (auto st = 0; st < steps; ++st) {
    field_fl << t + dt << " " << field(Ao, wA, cepds, Wenv, t + dt) << "\n";
    MV.field = field(Ao, wA, cepds, Wenv, t + dt);

    rkf.do_step(MV, ct, t, dt);

    t += dt;

    auto ctnrm = cblas_dznrm2(ct_sz, reinterpret_cast<double *>(&ct[0]), 1);

    for (auto &n : ct) {
      n /= ctnrm;
    }

    f_pop << std::setprecision(16) << t << " "
          << std::norm(ct[offs[pop_l] + pop_n]) << "\n";

    if (st % print == 0) {
      std::cout << field(Ao, wA, cepds, Wenv, t) << "\n";
      std::fstream f_ct;
      f_ct.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
      for (auto k = 0; k < ct_sz; ++k) {
        f_ct << k << "  " << ct[k].real() << "  " << ct[k].imag() << " "
             << std::norm(ct[k]) << "\n";
      }
      f_ct.close();
    }
  }
  field_fl.close();
  f_pop.close();
  return 0;
}
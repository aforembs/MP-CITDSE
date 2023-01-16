#include "tdse.hpp"

typedef std::vector<std::complex<double>> state_type;

class MatVecV {
public:
  int e_sz;
  int l_max;
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
    auto bt2 = std::complex<double>(1.0, 0.0);

    for (auto i = 0; i < e_sz; ++i) {
      dxdt[i] = eig[0][i] * mI * x[i];
    }

    cblas_zgemv(CblasRowMajor, CblasTrans, state_sz[1], state_sz[0],
                reinterpret_cast<double *>(&alp_flm),
                reinterpret_cast<double *>(cdipole[0]), state_sz[1],
                reinterpret_cast<double *>(&x[offs[1]]), 1,
                reinterpret_cast<double *>(&bt2),
                reinterpret_cast<double *>(&dxdt[0]), 1);

    for (auto l = 1; l < l_max; ++l) {
      cblas_zgemv(CblasRowMajor, CblasNoTrans, state_sz[l], state_sz[l - 1],
                  reinterpret_cast<double *>(&alp_fl),
                  reinterpret_cast<double *>(cdipole[l - 1]), state_sz[l],
                  reinterpret_cast<double *>(&x[offs[l - 1]]), 1,
                  reinterpret_cast<double *>(&bt2),
                  reinterpret_cast<double *>(&dxdt[offs[l]]), 1);

      cblas_zgemv(CblasRowMajor, CblasTrans, state_sz[l + 1], state_sz[l],
                  reinterpret_cast<double *>(&alp_flm),
                  reinterpret_cast<double *>(cdipole[l]), state_sz[l + 1],
                  reinterpret_cast<double *>(&x[offs[l + 1]]), 1,
                  reinterpret_cast<double *>(&bt2),
                  reinterpret_cast<double *>(&dxdt[offs[l]]), 1);
    }

    cblas_zgemv(CblasRowMajor, CblasNoTrans, state_sz[l_max],
                state_sz[l_max - 1], reinterpret_cast<double *>(&alp_fl),
                reinterpret_cast<double *>(cdipole[l_max - 1]), state_sz[l_max],
                reinterpret_cast<double *>(&x[offs[l_max - 1]]), 1,
                reinterpret_cast<double *>(&bt2),
                reinterpret_cast<double *>(&dxdt[offs[l_max]]), 1);
  }
};

class MatVecL {
public:
  int e_sz;
  int l_max;
  double field;
  std::vector<int> state_sz;
  std::vector<int> offs;
  std::vector<double *> eig;
  std::vector<std::complex<double> *> cdipole;

  void operator()(state_type &x, state_type &dxdt,
                  [[maybe_unused]] double t) const {
    constexpr std::complex<double> mI(0.0, -1.0);
    auto alp_fl = std::complex<double>(0.0, -field);
    auto bt2 = std::complex<double>(1.0, 0.0);

    for (auto i = 0; i < e_sz; ++i) {
      dxdt[i] = eig[0][i] * mI * x[i];
    }

    cblas_zgemv(CblasRowMajor, CblasTrans, state_sz[1], state_sz[0],
                reinterpret_cast<double *>(&alp_fl),
                reinterpret_cast<double *>(cdipole[0]), state_sz[1],
                reinterpret_cast<double *>(&x[offs[1]]), 1,
                reinterpret_cast<double *>(&bt2),
                reinterpret_cast<double *>(&dxdt[0]), 1);

    for (auto l = 1; l < l_max; ++l) {
      cblas_zgemv(CblasRowMajor, CblasNoTrans, state_sz[l], state_sz[l - 1],
                  reinterpret_cast<double *>(&alp_fl),
                  reinterpret_cast<double *>(cdipole[l - 1]), state_sz[l],
                  reinterpret_cast<double *>(&x[offs[l - 1]]), 1,
                  reinterpret_cast<double *>(&bt2),
                  reinterpret_cast<double *>(&dxdt[offs[l]]), 1);

      cblas_zgemv(CblasRowMajor, CblasTrans, state_sz[l + 1], state_sz[l],
                  reinterpret_cast<double *>(&alp_fl),
                  reinterpret_cast<double *>(cdipole[l]), state_sz[l + 1],
                  reinterpret_cast<double *>(&x[offs[l + 1]]), 1,
                  reinterpret_cast<double *>(&bt2),
                  reinterpret_cast<double *>(&dxdt[offs[l]]), 1);
    }

    cblas_zgemv(CblasRowMajor, CblasNoTrans, state_sz[l_max],
                state_sz[l_max - 1], reinterpret_cast<double *>(&alp_fl),
                reinterpret_cast<double *>(cdipole[l_max - 1]), state_sz[l_max],
                reinterpret_cast<double *>(&x[offs[l_max - 1]]), 1,
                reinterpret_cast<double *>(&bt2),
                reinterpret_cast<double *>(&dxdt[offs[l_max]]), 1);
  }
};

int tdse::propV(std::string output, int l_max, double t, double dt, int steps,
                fieldInit fieldst, fieldFcn field, double Io, double w,
                double cepd, int cycles, int e_sz, std::vector<int> &offs,
                std::vector<int> &state_sz, std::vector<double *> &eig,
                stvupt &dipoles, std::vector<std::complex<double>> &ct) {
  double IoA, wA, Ao, cepds, Wenv;
  pulse::toAU(Io, w, IoA, wA);

  int print = steps / 10;

  fieldst(IoA, wA, cepd, cycles, Ao, cepds, Wenv);

  state_type exv(e_sz), cvec(e_sz * e_sz), ctn(e_sz);

  std::fstream f_out, field_fl, f_pop;
  f_out.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
  for (auto k = 0; k < e_sz; ++k) {
    f_out << k << "  " << ct[k].real() << "  " << ct[k].imag() << "\n";
  }
  f_out.close();

  std::vector<state_type> cdipole;
  for (auto l = 0; l < l_max; ++l) {
    auto l_sz = state_sz[l];
    auto l1_sz = state_sz[l + 1];
    cdipole.push_back(state_type(l_sz * l1_sz));

    for (auto i = 0; i < l_sz; ++i) {
      for (auto jd = 0; jd < l1_sz; ++jd) {
        cdipole[l][i * l1_sz + jd] =
            std::complex<double>(dipoles[l]->at(i * l1_sz + jd), 0.0);
      }
    }
  }

  MatVecV MV;
  MV.e_sz = e_sz;
  MV.l_max = l_max;
  MV.state_sz = state_sz;
  MV.offs = offs;

  for (auto l = 0; l <= l_max; ++l) {
    MV.eig.push_back(eig[l]);
    MV.cdipole.push_back(cdipole[l].data());
  }

  boost::numeric::odeint::runge_kutta_fehlberg78<state_type> rkf;

  field_fl.open(output + "_field.dat", std::ios::out);
  f_pop.open(output + "_pop.dat", std::ios::out);
  for (auto st = 0; st < steps; ++st) {
    MV.field = field(Ao, wA, cepds, Wenv, t + dt);

    rkf.do_step(MV, ct, t, dt);

    t += dt;

    auto ctnrm = cblas_dznrm2(e_sz, reinterpret_cast<double *>(&ct[0]), 1);

    for (auto &n : ct) {
      n /= ctnrm;
    }

    f_pop << t << " " << std::norm(ct[0]) << "\n";

    if (st % print == 0) {
      std::cout << field(Ao, wA, cepds, Wenv, t) << "\n";
      std::fstream f_ct;
      f_ct.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
      for (auto k = 0; k < e_sz; ++k) {
        f_ct << k << "  " << ct[k].real() << "  " << ct[k].imag() << "\n";
      }
      f_ct.close();
    }
  }
  field_fl.close();
  f_pop.close();
  return 0;
}

int tdse::propL(std::string output, int l_max, double t, double dt, int steps,
                fieldInit fieldst, fieldFcn field, double Io, double w,
                double cepd, int cycles, int e_sz, std::vector<int> &offs,
                std::vector<int> &state_sz, std::vector<double *> &eig,
                stvupt &dipoles, std::vector<std::complex<double>> &ct) {
  double IoA, wA, Eo, cepds, Wenv;
  pulse::toAU(Io, w, IoA, wA);

  int print = steps / 10;

  fieldst(IoA, wA, cepd, cycles, Eo, cepds, Wenv);

  state_type exv(e_sz), cvec(e_sz * e_sz), ctn(e_sz);

  std::fstream f_out, field_fl, f_pop;
  f_out.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
  for (auto k = 0; k < e_sz; ++k) {
    f_out << k << "  " << ct[k].real() << "  " << ct[k].imag() << "\n";
  }
  f_out.close();

  std::vector<state_type> cdipole;
  for (auto l = 0; l < l_max; ++l) {
    auto l_sz = state_sz[l];
    auto l1_sz = state_sz[l + 1];
    cdipole.push_back(state_type(l_sz * l1_sz));

    for (auto i = 0; i < l_sz; ++i) {
      for (auto jd = 0; jd < l1_sz; ++jd) {
        cdipole[l][i * l1_sz + jd] =
            std::complex<double>(dipoles[l]->at(i * l1_sz + jd), 0.0);
      }
    }
  }

  MatVecL MV;
  MV.e_sz = e_sz;
  MV.l_max = l_max;
  MV.state_sz = state_sz;
  MV.offs = offs;

  for (auto l = 0; l <= l_max; ++l) {
    MV.eig.push_back(eig[l]);
    MV.cdipole.push_back(cdipole[l].data());
  }

  boost::numeric::odeint::runge_kutta_fehlberg78<state_type> rkf;

  field_fl.open(output + "_field.dat", std::ios::out);
  f_pop.open(output + "_pop.dat", std::ios::out);
  for (auto st = 0; st < steps; ++st) {
    MV.field = field(Eo, wA, cepds, Wenv, t + dt);

    rkf.do_step(MV, ct, t, dt);

    t += dt;

    auto ctnrm = cblas_dznrm2(e_sz, reinterpret_cast<double *>(&ct[0]), 1);

    for (auto &n : ct) {
      n /= ctnrm;
    }

    f_pop << t << " " << std::norm(ct[0]) << "\n";

    if (st % print == 0) {
      std::cout << field(Eo, wA, cepds, Wenv, t) << "\n";
      std::fstream f_ct;
      f_ct.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
      for (auto k = 0; k < e_sz; ++k) {
        f_ct << k << "  " << ct[k].real() << "  " << ct[k].imag() << "\n";
      }
      f_ct.close();
    }
  }
  field_fl.close();
  f_pop.close();
  return 0;
}
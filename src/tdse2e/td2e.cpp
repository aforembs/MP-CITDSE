#include "td2e.hpp"

typedef std::vector<std::complex<double>> state_type;

class MatVec {
public:
  int L_max;
  double field;
  std::vector<int> state_sz;
  std::vector<int> offs;
  std::vector<std::complex<double> *> cblock;
  std::vector<std::complex<double> *> cdipole;

  void operator()(state_type &x, state_type &dxdt,
                  [[maybe_unused]] double t) const {
    constexpr std::complex<double> mI(0.0, -1.0);
    auto alp_fl = std::complex<double>(field, 0.0);
    auto alpha = mI;
    auto beta = std::complex<double>(0.0, 0.0);
    auto bt2 = std::complex<double>(1.0, 0.0);

    cblas_zhemv(CblasRowMajor, CblasUpper, state_sz[0],
                reinterpret_cast<double *>(&alpha),
                reinterpret_cast<double *>(cblock[0]), state_sz[0],
                reinterpret_cast<double *>(&x[0]), 1,
                reinterpret_cast<double *>(&beta),
                reinterpret_cast<double *>(&dxdt[0]), 1);

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

      cblas_zhemv(CblasRowMajor, CblasUpper, state_sz[L],
                  reinterpret_cast<double *>(&alpha),
                  reinterpret_cast<double *>(cblock[L]), state_sz[L],
                  reinterpret_cast<double *>(&x[offs[L]]), 1,
                  reinterpret_cast<double *>(&bt2),
                  reinterpret_cast<double *>(&dxdt[offs[L]]), 1);

      cblas_zgemv(CblasRowMajor, CblasTrans, state_sz[L + 1], state_sz[L],
                  reinterpret_cast<double *>(&alp_fl),
                  reinterpret_cast<double *>(cdipole[L]), state_sz[L],
                  reinterpret_cast<double *>(&x[offs[L + 1]]), 1,
                  reinterpret_cast<double *>(&bt2),
                  reinterpret_cast<double *>(&dxdt[offs[L]]), 1);
    }

    cblas_zgemv(CblasRowMajor, CblasNoTrans, state_sz[L_max],
                state_sz[L_max - 1], reinterpret_cast<double *>(&alp_fl),
                reinterpret_cast<double *>(cdipole[L_max - 1]), state_sz[L_max - 1],
                reinterpret_cast<double *>(&x[offs[L_max - 1]]), 1,
                reinterpret_cast<double *>(&beta),
                reinterpret_cast<double *>(&dxdt[offs[L_max]]), 1);

    cblas_zhemv(CblasRowMajor, CblasUpper, state_sz[L_max],
                reinterpret_cast<double *>(&alpha),
                reinterpret_cast<double *>(cblock[L_max]), state_sz[L_max],
                reinterpret_cast<double *>(&x[offs[L_max]]), 1,
                reinterpret_cast<double *>(&bt2),
                reinterpret_cast<double *>(&dxdt[offs[L_max]]), 1);
  }
};

int td2e::prop(std::string output, int L_max, double t, double dt, int steps,
               fieldInit fieldst, fieldFcn field, double w, double Io,
               double cepd, int cycles, int ct_sz, std::vector<int> &offs,
               std::vector<int> &state_sz, stvupt &blocks, stvupt &dipoles,
               std::vector<std::complex<double>> &ct) {
  double IoA, wA, Ao, cepds, Wenv;
  pulse::ToAU(Io, w, IoA, wA);

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

  std::vector<state_type> cblock, cdipole;
  for (auto L = 0; L < L_max; ++L) {
    auto L_sz = state_sz[L];
    auto L1_sz = state_sz[L + 1];
    cblock.push_back(state_type(L_sz * L_sz));
    cdipole.push_back(state_type(L_sz * L1_sz));

    for (auto i = 0; i < L_sz; ++i) {
      for (auto jb = 0; jb < L_sz; ++jb) {
        cblock[L][i * L_sz + jb] =
            std::complex<double>(blocks[L].get()->at(i * L_sz + jb), 0.0);
      }
      for (auto jd = 0; jd < L1_sz; ++jd) {
        cdipole[L][i * L1_sz + jd] =
            std::complex<double>(dipoles[L].get()->at(i * L1_sz + jd), 0.0);
      }
    }
  }

  auto L_sz = state_sz[L_max];
  cblock.push_back(state_type(L_sz * L_sz));
  for (auto i = 0; i < L_sz * L_sz; ++i) {
    cblock[L_max][i] = std::complex<double>(blocks[L_max].get()->at(i), 0.0);
  }

  MatVec MV;
  MV.L_max = L_max;
  MV.state_sz = state_sz;
  MV.offs = offs;

  for (auto L = 0; L < L_max; ++L) {
    MV.cblock.push_back(cblock[L].data());
    MV.cdipole.push_back(cdipole[L].data());
  }
  MV.cblock.push_back(cblock[L_max].data());

  boost::numeric::odeint::runge_kutta_fehlberg78<state_type> rkf;

  std::vector<std::complex<double>> c_ground(state_sz[0]);
  for (auto i = 0; i < state_sz[0]; ++i) {
    c_ground[i] = ct[i];
  }

  auto ctnrm =
      cblas_dznrm2(state_sz[0], reinterpret_cast<double *>(&c_ground[0]), 1);

  for (auto &n : c_ground) {
    n /= ctnrm;
  }

  field_fl.open(output + "_field.dat", std::ios::out);
  f_pop.open(output + "_pop.dat", std::ios::out);

  std::complex<double> cdiff;
  cblas_zdotc_sub(state_sz[0], reinterpret_cast<double *>(&c_ground[0]), 1,
                  reinterpret_cast<double *>(&ct[0]), 1,
                  reinterpret_cast<double *>(&cdiff));

  f_pop << t << " " << std::norm(cdiff) << " "
        << cblas_dznrm2(state_sz[0], reinterpret_cast<double *>(&ct[0]), 1)
        << "\n";

  for (auto st = 0; st < steps; ++st) {
    field_fl << t + dt << " " << field(Ao, wA, cepds, Wenv, t + dt) << "\n";
    MV.field = field(Ao, wA, cepds, Wenv, t + dt);

    rkf.do_step(MV, ct, t, dt);

    t += dt;

    ctnrm = cblas_dznrm2(ct_sz, reinterpret_cast<double *>(&ct[0]), 1);

    for (auto &n : ct) {
      n /= ctnrm;
    }

    cblas_zdotc_sub(state_sz[0], reinterpret_cast<double *>(&c_ground[0]), 1,
                    reinterpret_cast<double *>(&ct[0]), 1,
                    reinterpret_cast<double *>(&cdiff));

    auto L0n = cblas_dznrm2(2500, reinterpret_cast<double *>(&ct[0]), 1);
    auto L1n = cblas_dznrm2(2500, reinterpret_cast<double *>(&ct[2500]), 1);
    auto L2n = cblas_dznrm2(2500, reinterpret_cast<double *>(&ct[5000]), 1);

    f_pop << std::setprecision(16) << t << " " << std::norm(cdiff) << " "
          << L0n * L0n << " " << L1n * L1n << " " << L2n * L2n << " "
          << L0n * L0n + L1n * L1n + L2n * L2n << "\n";

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
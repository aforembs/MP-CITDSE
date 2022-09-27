#include "gr2e.hpp"

class MatVec {
public:
  int L_sz;
  double *block;

  void operator()(std::vector<double> &x, std::vector<double> &dxdt,
                  [[maybe_unused]] double t) const {
    auto alpha = -1.0;
    auto beta = 0.0;

    cblas_dsymv(CblasRowMajor, CblasUpper, L_sz, alpha, block, L_sz, x.data(),
                1, beta, dxdt.data(), 1);
  }
};

int gr2e::prop(std::string output, int L_sz, double t, double dt,
               std::vector<double> &block, std::vector<double> &c0) {
  std::vector<double> res(L_sz);
  MatVec MV;
  MV.L_sz = L_sz;
  MV.block = block.data();

  c0.resize(L_sz);

  for (auto i = 0; i < L_sz; ++i) {
    c0[i] = 1.0;
  }

  auto c0nrm = cblas_dnrm2(L_sz, c0.data(), 1);

  for (auto &n : c0) {
    n /= c0nrm;
  }

  boost::numeric::odeint::runge_kutta_fehlberg78<std::vector<double>> rkf;

  std::fstream f_en(output + "_gr_conv.dat", std::ios::out);

  double tot_en_last = 0.0;
  while (true) {
    rkf.do_step(MV, c0, t, dt);

    t += dt;

    c0nrm = cblas_dnrm2(L_sz, c0.data(), 1);

    for (auto &n : c0) {
      n /= c0nrm;
    }

    auto tot_en = 0.0;

    cblas_dsymv(CblasRowMajor, CblasUpper, L_sz, 1.0, block.data(), L_sz,
                c0.data(), 1, 0.0, res.data(), 1);

    tot_en = cblas_ddot(L_sz, res.data(), 1, c0.data(), 1);
    f_en << t << " " << tot_en << "\n";
    if (std::abs(tot_en - tot_en_last) <
        std::numeric_limits<double>::epsilon()) {
      break;
    }
    tot_en_last = tot_en;
  }
  f_en.close();

  std::fstream f_c0;
  f_c0.open(output + "_c0.dat", std::ios::out);
  for (auto k = 0; k < L_sz; ++k) {
    f_c0 << k << "  " << c0[k] << "\n";
  }
  f_c0.close();

  std::cout << "Ground State Energy (a.u.): " << tot_en_last << "\n";

  return 0;
}
#include "tdse.hpp"

typedef std::vector<double> state_type;

class MatVecV {
public:
  int L_max;
  double field;
  std::vector<int> state_sz;
  std::vector<int> offs;
  std::vector<double *> eig;
  std::vector<double *> dipole;

  void operator()(state_type &x, state_type &dxdt,
                  [[maybe_unused]] double t) const {
    constexpr double beta = 0.0;
    constexpr double bt2 = 1.0;
    int off2m1 = 0;
    int off2 = 0;
    int off2p1 = offs[1] * 2;

    for (auto i = 0; i < state_sz[0]; ++i) {
      auto i2 = i * 2;
      dxdt[i2 + 1] = -eig[0][i] * x[i2];
      dxdt[i2] = eig[0][i] * x[i2 + 1];
    }

    cblas_dgemv(CblasColMajor, CblasTrans, state_sz[1], state_sz[0], -field,
                dipole[0], state_sz[1], &x[off2p1], 2, bt2, &dxdt[0], 2);
    cblas_dgemv(CblasColMajor, CblasTrans, state_sz[1], state_sz[0], -field,
                dipole[0], state_sz[1], &x[off2p1 + 1], 2, bt2, &dxdt[1], 2);

    for (auto L = 1; L < L_max; ++L) {
      off2m1 = offs[L - 1] * 2;
      off2 = offs[L] * 2;
      off2p1 = offs[L + 1] * 2;
      cblas_dgemv(CblasColMajor, CblasNoTrans, state_sz[L], state_sz[L - 1],
                  field, dipole[L - 1], state_sz[L], &x[off2m1], 2, beta,
                  &dxdt[off2], 2);
      cblas_dgemv(CblasColMajor, CblasNoTrans, state_sz[L], state_sz[L - 1],
                  field, dipole[L - 1], state_sz[L], &x[off2m1 + 1], 2,
                  beta, &dxdt[off2 + 1], 2);

      for (auto i = 0; i < state_sz[L]; ++i) {
        auto i2 = off2 + i * 2;
        dxdt[i2 + 1] = -eig[L][i] * x[i2] + bt2 * dxdt[i2 + 1];
        dxdt[i2] = eig[L][i] * x[i2 + 1] + bt2 * dxdt[i2];
      }

      cblas_dgemv(CblasColMajor, CblasTrans, state_sz[L + 1], state_sz[L],
                  -field, dipole[L], state_sz[L+1], &x[off2p1], 2, bt2,
                  &dxdt[off2], 2);
      cblas_dgemv(CblasColMajor, CblasTrans, state_sz[L + 1], state_sz[L],
                  -field, dipole[L], state_sz[L+1], &x[off2p1 + 1], 2, bt2,
                  &dxdt[off2 + 1], 2);
    }
    off2m1 = offs[L_max - 1] * 2;
    off2 = offs[L_max] * 2;

    cblas_dgemv(CblasColMajor, CblasNoTrans, state_sz[L_max],
                state_sz[L_max - 1], field, dipole[L_max - 1],
                state_sz[L_max], &x[off2m1], 2, beta, &dxdt[off2], 2);
    cblas_dgemv(CblasColMajor, CblasNoTrans, state_sz[L_max],
                state_sz[L_max - 1], field, dipole[L_max - 1],
                state_sz[L_max], &x[off2m1 + 1], 2, beta, &dxdt[off2 + 1],
                2);

    for (auto i = 0; i < state_sz[L_max]; ++i) {
      auto i2 = off2 + i * 2;
      dxdt[i2 + 1] = -eig[L_max][i] * x[i2] + bt2 * dxdt[i2 + 1];
      dxdt[i2] = eig[L_max][i] * x[i2 + 1] + bt2 * dxdt[i2];
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
  std::vector<double *> dipole;

  void operator()(state_type &x, state_type &dxdt,
                  [[maybe_unused]] double t) const {
    constexpr double beta = 0.0;
    constexpr double bt2 = 1.0;
    int off2m1 = 0;
    int off2 = 0;
    int off2p1 = offs[1] * 2;

    for (auto i = 0; i < state_sz[0]; ++i) {
      auto i2 = i * 2;
      dxdt[i2 + 1] = -eig[0][i] * x[i2];
      dxdt[i2] = eig[0][i] * x[i2 + 1];
    }

    cblas_dgemv(CblasColMajor, CblasTrans, state_sz[1], state_sz[0], field,
                dipole[0], state_sz[1], &x[off2p1], 2, bt2, &dxdt[1], 2);
    cblas_dgemv(CblasColMajor, CblasTrans, state_sz[1], state_sz[0], field,
                dipole[0], state_sz[1], &x[off2p1 + 1], 2, bt2, &dxdt[0], 2);

    for (auto L = 1; L < L_max; ++L) {
      off2m1 = offs[L - 1] * 2;
      off2 = offs[L] * 2;
      off2p1 = offs[L + 1] * 2;
      cblas_dgemv(CblasColMajor, CblasNoTrans, state_sz[L], state_sz[L - 1],
                  field, dipole[L - 1], state_sz[L], &x[off2m1], 2, beta,
                  &dxdt[off2 + 1], 2);
      cblas_dgemv(CblasColMajor, CblasNoTrans, state_sz[L], state_sz[L - 1],
                  field, dipole[L - 1], state_sz[L], &x[off2m1 + 1], 2,
                  beta, &dxdt[off2], 2);

      for (auto i = 0; i < state_sz[L]; ++i) {
        auto i2 = off2 + i * 2;
        dxdt[i2 + 1] = -eig[L][i] * x[i2] + bt2 * dxdt[i2 + 1];
        dxdt[i2] = eig[L][i] * x[i2 + 1] + bt2 * dxdt[i2];
      }

      cblas_dgemv(CblasColMajor, CblasTrans, state_sz[L + 1], state_sz[L],
                  field, dipole[L], state_sz[L+1], &x[off2p1], 2, bt2,
                  &dxdt[off2 + 1], 2);
      cblas_dgemv(CblasColMajor, CblasTrans, state_sz[L + 1], state_sz[L],
                  field, dipole[L], state_sz[L+1], &x[off2p1 + 1], 2, bt2,
                  &dxdt[off2], 2);
    }
    off2m1 = offs[L_max - 1] * 2;
    off2 = offs[L_max] * 2;

    cblas_dgemv(CblasColMajor, CblasNoTrans, state_sz[L_max],
                state_sz[L_max - 1], field, dipole[L_max - 1],
                state_sz[L_max], &x[off2m1], 2, beta, &dxdt[off2 + 1], 2);
    cblas_dgemv(CblasColMajor, CblasNoTrans, state_sz[L_max],
                state_sz[L_max - 1], field, dipole[L_max - 1],
                state_sz[L_max], &x[off2m1 + 1], 2, beta, &dxdt[off2], 2);

    for (auto i = 0; i < state_sz[L_max]; ++i) {
      auto i2 = off2 + i * 2;
      dxdt[i2 + 1] = -eig[L_max][i] * x[i2] + bt2 * dxdt[i2 + 1];
      dxdt[i2] = eig[L_max][i] * x[i2 + 1] + bt2 * dxdt[i2];
    }
  }
};

int tdse::propV(std::string output, int L_max, double t, double dt, int steps,
                int pop_n, int pop_l, fieldFcn field, pulse::params &pars,
                int ct_sz, std::vector<int> &offs, std::vector<int> &state_sz,
                stvupt &eig, stvupt &dipoles, std::vector<double> &ct) {
  int print = steps / 10;

  std::fstream f_out, field_fl, f_pop;
  f_out.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
  f_out << "#Index, Re(c(t)), Im(c(t)), |c(t)|^2\n";
  for (auto k = 0; k < ct_sz * 2; k += 2) {
    f_out << k / 2 << "  " << ct[k] << "  " << ct[k + 1] << "\n";
  }
  f_out.close();

  MatVecV MV;
  MV.L_max = L_max;
  MV.state_sz = state_sz;
  MV.offs = offs;

  for (auto L = 0; L < L_max; ++L) {
    MV.eig.push_back(eig[L]->data());
    MV.dipole.push_back(dipoles[L]->data());
  }
  MV.eig.push_back(eig[L_max]->data());

  boost::numeric::odeint::runge_kutta_fehlberg78<state_type> rkf;

  field_fl.open(output + "_field.dat", std::ios::out);
  field_fl << "#time (a.u.), field (A(t) or E(t)) (a.u.)\n";

  f_pop.open(output + "_pop" + std::to_string(pop_n + 1) +
                 std::to_string(pop_l) + ".dat",
             std::ios::out);

  f_pop << "#time (a.u.), population\n";
  f_pop << std::setprecision(16) << t << " "
        << std::norm(std::complex<double>(ct[offs[pop_l] * 2 + pop_n * 2],
                                          ct[offs[pop_l] * 2 + pop_n * 2 + 1]))
        << "\n";

  for (auto st = 0; st < steps; ++st) {
    MV.field = field(pars, t + dt);
    field_fl << t + dt << " " << MV.field << "\n";

    rkf.do_step(MV, ct, t, dt);

    t += dt;

    auto ctnrm = cblas_dznrm2(ct_sz, &ct[0], 1);

    for (auto &n : ct) {
      n /= ctnrm;
    }

    f_pop << std::setprecision(16) << t << " "
          << std::norm(
                 std::complex<double>(ct[offs[pop_l] * 2 + pop_n * 2],
                                      ct[offs[pop_l] * 2 + pop_n * 2 + 1]))
          << "\n";

    if (st % print == 0) {
      std::cout << MV.field << "\n";
      std::fstream f_ct;
      f_ct.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
      f_ct << "#Index, Re(c(t)), Im(c(t)), |c(t)|^2\n";
      for (auto k = 0; k < ct_sz * 2; k += 2) {
        f_ct << k / 2 << "  " << ct[k] << "  " << ct[k + 1] << "\n";
      }
      f_ct.close();
    }
  }
  field_fl.close();
  f_pop.close();
  return 0;
}

int tdse::propL(std::string output, int L_max, double t, double dt, int steps,
                int pop_n, int pop_l, fieldFcn field, pulse::params &pars,
                int ct_sz, std::vector<int> &offs, std::vector<int> &state_sz,
                stvupt &eig, stvupt &dipoles, std::vector<double> &ct) {
  int print = steps / 10;

  std::fstream f_out, field_fl, f_pop;
  f_out.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
  f_out << "#Index, Re(c(t)), Im(c(t)), |c(t)|^2\n";
  for (auto k = 0; k < ct_sz * 2; k += 2) {
    f_out << k / 2 << "  " << ct[k] << "  " << ct[k + 1] << "\n";
  }
  f_out.close();

  MatVecL MV;
  MV.L_max = L_max;
  MV.state_sz = state_sz;
  MV.offs = offs;

  for (auto L = 0; L < L_max; ++L) {
    MV.eig.push_back(eig[L]->data());
    MV.dipole.push_back(dipoles[L]->data());
  }
  MV.eig.push_back(eig[L_max]->data());

  boost::numeric::odeint::runge_kutta_fehlberg78<state_type> rkf;

  field_fl.open(output + "_field.dat", std::ios::out);
  field_fl << "#time (a.u.), field (A(t) or E(t)) (a.u.)\n";

  f_pop.open(output + "_pop" + std::to_string(pop_n + 1) +
                 std::to_string(pop_l) + ".dat",
             std::ios::out);

  f_pop << "#time (a.u.), population\n";
  f_pop << std::setprecision(16) << t << " "
        << std::norm(std::complex<double>(ct[offs[pop_l] * 2 + pop_n * 2],
                                          ct[offs[pop_l] * 2 + pop_n * 2 + 1]))
        << "\n";

  for (auto st = 0; st < steps; ++st) {
    MV.field = field(pars, t + dt);
    field_fl << t + dt << " " << MV.field << "\n";

    rkf.do_step(MV, ct, t, dt);

    t += dt;

    auto ctnrm = cblas_dznrm2(ct_sz, &ct[0], 1);

    for (auto &n : ct) {
      n /= ctnrm;
    }

    f_pop << std::setprecision(16) << t << " "
          << std::norm(
                 std::complex<double>(ct[offs[pop_l] * 2 + pop_n * 2],
                                      ct[offs[pop_l] * 2 + pop_n * 2 + 1]))
          << "\n";

    if (st % print == 0) {
      std::cout << MV.field << "\n";
      std::fstream f_ct;
      f_ct.open(output + "_ct_" + std::to_string(t) + ".dat", std::ios::out);
      f_ct << "#Index, Re(c(t)), Im(c(t)), |c(t)|^2\n";
      for (auto k = 0; k < ct_sz * 2; k += 2) {
        f_ct << k / 2 << "  " << ct[k] << "  " << ct[k + 1] << "\n";
      }
      f_ct.close();
    }
  }
  field_fl.close();
  f_pop.close();
  return 0;
}

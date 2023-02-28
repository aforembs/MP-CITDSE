#include "integrator.hpp"

double dmx_int::tvelGL(int qsz, int lc_sz, int n1, int l1, int n2, int l2,
                       std::vector<double> &qx, std::vector<double> &qw,
                       std::vector<double> &wfn, std::vector<double> &wfnp) {
  double t_ab = 0.0;
  int p1 = l1 * lc_sz + n1 * qsz;
  int p2 = l2 * lc_sz + n2 * qsz;
  double v2 = (l1 * (l1 + 1) - l2 * (l2 + 1)) * 0.5;
  double pra;

  for (auto i = 0; i < qsz; ++i) {
    pra = wfn[p1 + i];

    t_ab += qw[i] * (pra * wfnp[n2 * qsz + i] - v2 * pra * wfn[p2 + i] / qx[i]);
  }
  return t_ab;
}

double dmx_int::tlenGL(int qsz, int lc_sz, int n1, int l1, int n2, int l2,
                       std::vector<double> &qx, std::vector<double> &qw,
                       std::vector<double> &wfn) {
  double t_ab = 0.0;
  int p1 = l1 * lc_sz + n1 * qsz;
  int p2 = l2 * lc_sz + n2 * qsz;

  for (auto i = 0; i < qsz; ++i) {
    t_ab += qw[i] * qx[i] * wfn[p1 + i] * wfn[p2 + i];
  }

  return t_ab;
}
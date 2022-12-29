#include "integrator_new.hpp"

namespace glw {
const double GLobw3[2] = {0.333333333333333333333333333333,
                          1.333333333333333333333333333333};
const double GLobw4[2] = {0.166666666666666666666666666667,
                          0.83333333333333333333333333333};
const double GLobw5[3] = {0.1, 0.544444444444444444444444444444,
                          0.71111111111111111111111111111};
const double GLobw6[3] = {0.0666666666666666666666666666667,
                          0.378474956297846980316612808212,
                          0.554858377035486353016720525121};
const double GLobw7[4] = {
    0.047619047619047619047619047619, 0.27682604736156594801070040629,
    0.431745381209862623417871022281, 0.487619047619047619047619047619};
const double GLobw8[4] = {
    0.0357142857142857142857142857143, 0.21070422714350603938299206578,
    0.34112269248350436476424067711, 0.4124587946587038815670529714};
const double GLobw9[5] = {
    0.0277777777777777777777777777778, 0.16549536156080552504633972003,
    0.274538712500161735280705618579, 0.34642851097304634511513153214,
    0.371519274376417233560090702948};
const double GLobw10[5] = {
    0.0222222222222222222222222222222, 0.133305990851070111126227170755,
    0.224889342063126452119457821731,  0.292042683679683757875582257374,
    0.327539761183897456656510527917};

const double *Globw[9] = {0, GLobw3, GLobw4, GLobw5,
                           GLobw6, GLobw7, GLobw8, GLobw9, GLobw10};
} // namespace glw

double fsltrLob(int k, int n, int bo, int glq_pt, int lc_sz, int n1, int l1,
                int n2, int l2, int n3, int l3, int n4, int l4,
                std::vector<double> &gl_w, std::vector<double> &gl_x,
                std::vector<double> &kkn, std::vector<double> &rk,
                std::vector<double> &rk_mid, std::vector<double> &wfn_o,
                std::vector<double> &wfn_i) {
  int kp1 = k + 1;
  int nqpt = n * glq_pt;
  int p2i = l2 * lc_sz + n2 * nqpt;
  int p2p = l4 * lc_sz + n4 * nqpt;
  int p1i = l1 * lc_sz + n1 * nqpt;
  int p1p = l3 * lc_sz + n3 * nqpt;
  int i1, i1bo;
  double dl, sl, loc_GL, r1, pr2, pr2a, chi;
  double rm1 = 0.0;
  double dlob;
  double fm1j = 0.0, fm1q = 0.0;
  double pr1k, pr1km;

  // first calculate Qk for all of 0->R
  double Qk = 0.0;

  for (auto i = bo - 1; i < n; ++i) {
    i1 = i + 1;
    i1bo = (i1 - bo) * glq_pt;
    dl = (kkn[i1] - kkn[i]) * 0.5;
    sl = (kkn[i1] + kkn[i]) * 0.5;

    // need to add inner lobatto points
    for (auto p = 0; p < glq_pt; ++p) {
      r1 = dl * gl_x[p] + sl;
      dlob = (r1 - rm1) * 0.5;
      pr2 = wfn_o[p2i + p + i1bo] * wfn_o[p2p + p + i1bo];

      Qk += dlob * (Lobo * fm1q +
                    Lobi * (1.0 / rk_mid[p + i1bo + kp1 * nqpt]) *
                        wfn_i[p2i + p + i1bo] * wfn_i[p2p + p + i1bo] +
                    Lobo * (1.0 / rk[p + i1bo + kp1 * nqpt]) * pr2);

      fm1q = (1.0 / rk[p + i1bo + kp1 * nqpt]) * pr2;
      rm1 = r1;
    }
  }

  double Jk = 0.0;
  double Fk = 0.0;
  rm1 = 0.0;
  fm1q = 0.0;

  for (auto i = bo - 1; i < n; ++i) {
    i1 = i + 1;
    i1bo = (i1 - bo) * glq_pt;
    dl = (kkn[i1] - kkn[i]) * 0.5;
    sl = (kkn[i1] + kkn[i]) * 0.5;
    loc_GL = 0.0;

    for (auto p = 0; p < glq_pt; ++p) {
      r1 = dl * gl_x[p] + sl;
      dlob = (r1 - rm1) * 0.5;

      pr2 = wfn_o[p2i + p + i1bo] * wfn_o[p2p + p + i1bo];
      pr2a = wfn_i[p2i + p + i1bo] * wfn_i[p2p + p + i1bo];

      // chi(r1)
      pr1k = rk[p + i1bo + k * nqpt];
      pr1km = 1.0 / rk[p + i1bo + kp1 * nqpt];
      Jk += dlob * (Lobo * fm1j + Lobi * rk_mid[p + i1bo + k * nqpt] * pr2a +
                    Lobo * pr1k * pr2);
      Qk -= dlob *
            (Lobo * fm1q + Lobi * (1.0 / rk_mid[p + i1bo + kp1 * nqpt]) * pr2a +
             Lobo * pr1km * pr2);
      chi = pr1km * Jk + pr1k * Qk;

      // Glq outer
      loc_GL += gl_w[p] * wfn_o[p1i + p + i1bo] * wfn_o[p1p + p + i1bo] * chi;

      rm1 = r1;
      fm1j = pr1k * pr2;
      fm1q = pr1km * pr2;
    }
    Fk += dl * loc_GL;
  }
  return Fk;
}

double fsltrTrap(int k, int n, int bo, int glq_pt, int lc_sz, int n1, int l1,
                 int n2, int l2, int n3, int l3, int n4, int l4,
                 std::vector<double> &gl_w, std::vector<double> &gl_x,
                 std::vector<double> &kkn, std::vector<double> &rk,
                 std::vector<double> &wfn_o) {
  int kp1 = k + 1;
  int nqpt = n * glq_pt;
  int p2i = l2 * lc_sz + n2 * nqpt;
  int p2p = l4 * lc_sz + n4 * nqpt;
  int p1i = l1 * lc_sz + n1 * nqpt;
  int p1p = l3 * lc_sz + n3 * nqpt;
  double dl, sl, loc_GL, r2, r1, pr2, chi;

  // first calculate Qk for all of 0->R
  double Qk = 0.0;
  double Jk = 0;
  double Fk = 0;
  double rm1 = 0.0;
  double dltr;
  double fm1j = 0.0, fm1q = 0.0;
  double pr1k, pr1km;
  auto i1 = 0, i1bo = 0;

  for (auto i = bo - 1; i < n; ++i) {
    i1 = i + 1;
    dl = (kkn[i1] - kkn[i]) * 0.5;
    sl = (kkn[i1] + kkn[i]) * 0.5;
    i1bo = (i1 - bo) * glq_pt;

    for (auto p = 0; p < glq_pt; ++p) {
      r2 = dl * gl_x[p] + sl;
      dltr = (r2 - rm1) * 0.5;
      pr2 = wfn_o[p2i + p + i1bo] * wfn_o[p2p + p + i1bo];

      Qk += dltr * (fm1q + (1.0 / rk[p + i1bo + kp1 * nqpt]) * pr2);

      fm1q = (1.0 / rk[p + i1bo + kp1 * nqpt]) * pr2;
      rm1 = r2;
    }
  }

  rm1 = 0.0;
  fm1q = 0.0;
  for (auto i = bo - 1; i < n; ++i) {
    i1 = i + 1;
    i1bo = (i1 - bo) * glq_pt;
    dl = (kkn[i + 1] - kkn[i]) * 0.5;
    sl = (kkn[i + 1] + kkn[i]) * 0.5;
    loc_GL = 0.0;

    for (auto p = 0; p < glq_pt; ++p) {
      r1 = dl * gl_x[p] + sl;
      dltr = (r1 - rm1) * 0.5;
      pr2 = wfn_o[p2i + p + i1bo] * wfn_o[p2p + p + i1bo];

      // chi(r1)
      pr1k = rk[p + i1bo + k * nqpt];
      pr1km = 1.0 / rk[p + i1bo + kp1 * nqpt];
      Jk += dltr * (fm1j + pr1k * pr2);
      Qk -= dltr * (fm1q + pr1km * pr2);
      chi = pr1km * Jk + pr1k * Qk;

      loc_GL += gl_w[p] * wfn_o[p1i + p + i1bo] * wfn_o[p1p + p + i1bo] * chi;

      rm1 = r1;
      fm1j = pr1k * pr2;
      fm1q = pr1km * pr2;
    }
    Fk += dl * loc_GL;
  }
  return Fk;
}
#include "integrator.hpp"

double FsltrGL2(int k, int n, int bo, int gl2,
                // int na, int la, int nb, int lb,
                // int nc, int lc, int nd, int ld,
                std::vector<double> &gl_ow, std::vector<double> &gl_ox,
                std::vector<double> &gl_iw, std::vector<double> &gl_ix,
                std::vector<double> &kkn, std::vector<double> &Bsp,
                std::vector<double> &Gsp,
                // std::vector<int> &offset,
                std::vector<double> &Cl1i_pt, std::vector<double> &Cl1p_pt,
                std::vector<double> &Cl2i_pt, std::vector<double> &Cl2p_pt) {
  int kp1 = k + 1;
  // double *Cl1i_pt=&Cf[offset[la]];
  // double *Cl1p_pt=&Cf[offset[lc]];
  // double *Cl2i_pt=&Cf[offset[lb]];
  // double *Cl2p_pt=&Cf[offset[ld]];
  double Pl1i = 0, Pl1p = 0, Pl2i = 0, Pl2p = 0;
  double Pl2igli = 0, Pl2pgli = 0;
  double dl, sl, loc_GL, r2, r1, gl2r, pr2, dli, chi;

  // int nna=na*n, nnb=nb*n, nnc=nc*n, nnd=nd*n;

  // first calculate Qk for all of 0->R
  double Qk = 0.0;
  double Jk = 0;
  double Fk = 0;

  double rm1 = 0.0;
  double jki = 0.0, qki = 0.0;
  auto sp = 0;
  auto ibo1j = 0, jbopi = 0;
  auto ai = 0;

  for (auto i = bo - 1; i < n; ++i) {
    dl = (kkn[i + 1] - kkn[i]) * 0.5;
    sl = (kkn[i + 1] + kkn[i]) * 0.5;
    loc_GL = 0.0;

    for (auto p = 0; p < bo; ++p) {
      r2 = dl * gl_ox[p] + sl;
      Pl2i = 0;
      Pl2p = 0;

      for (auto j = 0; j < bo; ++j) {
        jbopi = j + bo * (p + i * bo);
        ibo1j = i - bo + 1 + j;
        Pl2i += Cl2i_pt[ibo1j] * Bsp[jbopi];
        Pl2p += Cl2p_pt[ibo1j] * Bsp[jbopi];
      }
      loc_GL += gl_ow[p] * pow(r2, -kp1) * Pl2i * Pl2p;
    }
    Qk += dl * loc_GL;
  }

  std::ofstream outFile("dat/slt_test.dat", std::ofstream::out);

  for (auto i = bo - 1; i < n; ++i) {
    dl = (kkn[i + 1] - kkn[i]) * 0.5;
    sl = (kkn[i + 1] + kkn[i]) * 0.5;
    loc_GL = 0.0;

    sp = 0;
    for (auto p = 0; p < bo; ++p, sp += gl2) {
      r1 = dl * gl_ox[p] + sl;
      Pl1i = 0;
      Pl1p = 0;

      for (auto j = 0; j < bo; ++j) {
        ibo1j = i - bo + 1 + j;
        jbopi = j + bo * (p + i * bo);
        Pl1i += Cl1i_pt[ibo1j] * Bsp[jbopi];
        Pl1p += Cl1p_pt[ibo1j] * Bsp[jbopi];
      }

      // second, inner gaussian quadrature
      jki = 0;
      qki = 0;
      dli = (r1 - rm1) * 0.5;
      for (auto p2 = 0; p2 < gl2; p2 += 4) {
        Pl2igli = 0;
        Pl2pgli = 0;
        gl2r = dli * gl_ix[p2] + (r1 + rm1) * 0.5;
        for (auto j = 0; j < bo; ++j) {
          ibo1j = i - bo + 1 + j;
          jbopi = j + bo * (sp + p2 + i * gl2 * bo);
          ai = ibo1j - (kkn[i] > gl2r);
          Pl2igli += Cl2i_pt[ai] * Gsp[jbopi];
          Pl2pgli += Cl2p_pt[ai] * Gsp[jbopi];
        }
        pr2 = Pl2igli * Pl2pgli;
        jki += gl_iw[p2] * pow(gl2r, k) * pr2;
        qki += gl_iw[p2] * pow(gl2r, -kp1) * pr2;

        Pl2igli = 0;
        Pl2pgli = 0;
        gl2r = dli * gl_ix[p2 + 1] + (r1 + rm1) * 0.5;
        for (auto j = 0; j < bo; ++j) {
          ibo1j = i - bo + 1 + j;
          jbopi = j + bo * (sp + p2 + 1 + i * gl2 * bo);
          ai = ibo1j - (kkn[i] > gl2r);
          Pl2igli += Cl2i_pt[ai] * Gsp[jbopi];
          Pl2pgli += Cl2p_pt[ai] * Gsp[jbopi];
        }
        pr2 = Pl2igli * Pl2pgli;
        jki += gl_iw[p2 + 1] * pow(gl2r, k) * pr2;
        qki += gl_iw[p2 + 1] * pow(gl2r, -kp1) * pr2;

        Pl2igli = 0;
        Pl2pgli = 0;
        gl2r = dli * gl_ix[p2 + 2] + (r1 + rm1) * 0.5;
        for (auto j = 0; j < bo; ++j) {
          ibo1j = i - bo + 1 + j;
          jbopi = j + bo * (sp + p2 + 2 + i * gl2 * bo);
          ai = ibo1j - (kkn[i] > gl2r);
          Pl2igli += Cl2i_pt[ai] * Gsp[jbopi];
          Pl2pgli += Cl2p_pt[ai] * Gsp[jbopi];
        }
        pr2 = Pl2igli * Pl2pgli;
        jki += gl_iw[p2 + 2] * pow(gl2r, k) * pr2;
        qki += gl_iw[p2 + 2] * pow(gl2r, -kp1) * pr2;

        Pl2igli = 0;
        Pl2pgli = 0;
        gl2r = dli * gl_ix[p2 + 3] + (r1 + rm1) * 0.5;
        for (auto j = 0; j < bo; ++j) {
          ibo1j = i - bo + 1 + j;
          jbopi = j + bo * (sp + p2 + 3 + i * gl2 * bo);
          ai = ibo1j - (kkn[i] > gl2r);
          Pl2igli += Cl2i_pt[ai] * Gsp[jbopi];
          Pl2pgli += Cl2p_pt[ai] * Gsp[jbopi];
        }
        pr2 = Pl2igli * Pl2pgli;
        jki += gl_iw[p2 + 3] * pow(gl2r, k) * pr2;
        qki += gl_iw[p2 + 3] * pow(gl2r, -kp1) * pr2;
      }

      // chi(r1)
      Jk += dli * jki;
      Qk -= dli * qki;
      chi = pow(r1, -kp1) * Jk + pow(r1, k) * Qk;

      outFile << r1 << " " << Jk << " " << Qk << " " << chi << "\n";

      // Glq outer
      loc_GL += gl_ow[p] * Pl1i * Pl1p * chi;

      rm1 = r1;
    }
    Fk += dl * loc_GL;
  }

  return Fk;
}

double FsltrMM(int k, int n, int bo, int glq_pt, int lc_sz, int n1, int l1,
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
  int i1, i1bo, pt;
  double dl, sl, loc_GL, r2, r1, pr2, pr2a, chi;
  double dlob;

  // first calculate Qk for all of 0->R
  double Qk = 0.0;
  double Jk = 0;
  double Fk = 0;
  double rm1 = 0.0;
  double fm1q = 0.0;

  constexpr double Lob3o = 0.3333333333333333333333333e0;
  constexpr double Lob3i = 0.1333333333333333333333333e1;

  for (auto i = bo - 1; i < n; ++i) {
    i1 = i + 1;
    i1bo = (i1 - bo) * glq_pt;
    dl = (kkn[i1] - kkn[i]) * 0.5;
    sl = (kkn[i1] + kkn[i]) * 0.5;

    // first 2 intervals
    pt = 0;
    r2 = dl * gl_x[pt] + sl;
    dlob = (r2 - rm1) * 0.5;
    pr2 = wfn_o[p2i + pt + i1bo] * wfn_o[p2p + pt + i1bo];

    Qk += dlob * (fm1q + (1.0 / rk[pt + i1bo + kp1 * nqpt]) * pr2);

    fm1q = (1.0 / rk[pt + i1bo + kp1 * nqpt]) * pr2;
    rm1 = r2;

    pt = 1;
    r2 = dl * gl_x[pt] + sl;
    dlob = (r2 - rm1) * 0.5;
    pr2 = wfn_o[p2i + pt + i1bo] * wfn_o[p2p + pt + i1bo];

    Qk += dlob * (fm1q + (1.0 / rk[pt + i1bo + kp1 * nqpt]) * pr2);

    fm1q = (1.0 / rk[pt + i1bo + kp1 * nqpt]) * pr2;
    rm1 = r2;

    // middle intervals
    pt = glq_pt - 1;
    for (auto p = 2; p < pt; ++p) {
      r2 = dl * gl_x[p] + sl;
      dlob = (r2 - rm1) * 0.5;
      pr2 = wfn_o[p2i + p + i1bo] * wfn_o[p2p + p + i1bo];

      Qk += dlob * (Lob3o * fm1q +
                    Lob3i * (1.0 / rk_mid[p + i1bo + kp1 * nqpt]) *
                        wfn_i[p2i + p + i1bo] * wfn_i[p2p + p + i1bo] +
                    Lob3o * (1.0 / rk[p + i1bo + kp1 * nqpt]) * pr2);

      fm1q = (1.0 / rk[p + i1bo + kp1 * nqpt]) * pr2;
      rm1 = r2;
    }

    // last interval
    r2 = dl * gl_x[pt] + sl;
    dlob = (r2 - rm1) * 0.5;
    pr2 = wfn_o[p2i + pt + i1bo] * wfn_o[p2p + pt + i1bo];

    Qk += dlob * (fm1q + (1.0 / rk[pt + i1bo + kp1 * nqpt]) * pr2);

    fm1q = (1.0 / rk[pt + i1bo + kp1 * nqpt]) * pr2;
    rm1 = r2;
  }

  rm1 = 0.0;
  fm1q = 0.0;
  double fm1j = 0.0;
  double pr1k, pr1km;

  for (auto i = bo - 1; i < n; ++i) {
    i1 = i + 1;
    i1bo = (i1 - bo) * glq_pt;
    dl = (kkn[i + 1] - kkn[i]) * 0.5;
    sl = (kkn[i + 1] + kkn[i]) * 0.5;
    loc_GL = 0.0;

    pt = 0;
    r1 = dl * gl_x[pt] + sl;
    dlob = (r1 - rm1) * 0.5;
    pr2 = wfn_o[p2i + pt + i1bo] * wfn_o[p2p + pt + i1bo];

    // chi(r1)
    pr1k = rk[pt + i1bo + k * nqpt];
    pr1km = 1.0 / rk[pt + i1bo + kp1 * nqpt];
    Jk += dlob * (fm1j + pr1k * pr2);
    Qk -= dlob * (fm1q + pr1km * pr2);
    chi = pr1km * Jk + pr1k * Qk;

    loc_GL += gl_w[pt] * wfn_o[p1i + pt + i1bo] * wfn_o[p1p + pt + i1bo] * chi;

    rm1 = r1;
    fm1j = pr1k * pr2;
    fm1q = pr1km * pr2;

    pt = 1;
    r1 = dl * gl_x[pt] + sl;
    dlob = (r1 - rm1) * 0.5;
    pr2 = wfn_o[p2i + pt + i1bo] * wfn_o[p2p + pt + i1bo];

    // chi(r1)
    pr1k = rk[pt + i1bo + k * nqpt];
    pr1km = 1.0 / rk[pt + i1bo + kp1 * nqpt];
    Jk += dlob * (fm1j + pr1k * pr2);
    Qk -= dlob * (fm1q + pr1km * pr2);
    chi = pr1km * Jk + pr1k * Qk;

    loc_GL += gl_w[pt] * wfn_o[p1i + pt + i1bo] * wfn_o[p1p + pt + i1bo] * chi;

    rm1 = r1;
    fm1j = pr1k * pr2;
    fm1q = pr1km * pr2;

    pt = glq_pt - 1;
    for (auto p = 2; p < pt; ++p) {
      r1 = dl * gl_x[p] + sl;
      dlob = (r1 - rm1) * 0.5;

      pr2 = wfn_o[p2i + p + i1bo] * wfn_o[p2p + p + i1bo];
      pr2a = wfn_i[p2i + p + i1bo] * wfn_i[p2p + p + i1bo];

      // chi(r1)
      pr1k = rk[p + i1bo + k * nqpt];
      pr1km = 1.0 / rk[p + i1bo + kp1 * nqpt];
      Jk += dlob * (Lob3o * fm1j + Lob3i * rk_mid[p + i1bo + k * nqpt] * pr2a +
                    Lob3o * pr1k * pr2);
      Qk -= dlob * (Lob3o * fm1q +
                    Lob3i * (1.0 / rk_mid[p + i1bo + kp1 * nqpt]) * pr2a +
                    Lob3o * pr1km * pr2);
      chi = pr1km * Jk + pr1k * Qk;

      // Glq outer
      loc_GL += gl_w[p] * wfn_o[p1i + p + i1bo] * wfn_o[p1p + p + i1bo] * chi;

      rm1 = r1;
      fm1j = pr1k * pr2;
      fm1q = pr1km * pr2;
    }

    r1 = dl * gl_x[pt] + sl;
    dlob = (r1 - rm1) * 0.5;
    pr2 = wfn_o[p2i + pt + i1bo] * wfn_o[p2p + pt + i1bo];

    // chi(r1)
    pr1k = rk[pt + i1bo + k * nqpt];
    pr1km = 1.0 / rk[pt + i1bo + kp1 * nqpt];
    Jk += dlob * (fm1j + pr1k * pr2);
    Qk -= dlob * (fm1q + pr1km * pr2);
    chi = pr1km * Jk + pr1k * Qk;

    loc_GL += gl_w[pt] * wfn_o[p1i + pt + i1bo] * wfn_o[p1p + pt + i1bo] * chi;

    rm1 = r1;
    fm1j = pr1k * pr2;
    fm1q = pr1km * pr2;

    Fk += dl * loc_GL;
  }

  return Fk;
}

double FsltrLob4GL(int k, int n, int bo, int glq_pt, int lc_sz, int n1, int l1,
                   int n2, int l2, int n3, int l3, int n4, int l4,
                   std::vector<double> &gl_w, std::vector<double> &gl_x,
                   std::vector<double> &kkn, std::vector<double> &rk,
                   std::vector<double> &rk_mid, std::vector<double> &wfn_o,
                   std::vector<double> &wfn_i) {
  int kp1 = k + 1;
  int nqpt = n * glq_pt;
  int nqpt2 = 2 * nqpt;
  int p2i = l2 * lc_sz + n2 * nqpt;
  int p2p = l4 * lc_sz + n4 * nqpt;
  int p2ii = l2 * lc_sz * 2 + n2 * nqpt2;
  int p2pi = l4 * lc_sz * 2 + n4 * nqpt2;
  int p1i = l1 * lc_sz + n1 * nqpt;
  int p1p = l3 * lc_sz + n3 * nqpt;
  int i1, i1bo, i2bo;
  double dl, sl, loc_GL, r2, r1, pr2, pr2a, pr2b, chi;
  double dlob;

  // first calculate Qk for all of 0->R
  double Qk = 0.0;
  double Jk = 0;
  double Fk = 0;
  double rm1 = 0.0;
  double fm1q = 0.0;
  auto sp = 0;

  constexpr double Lob4o = 0.1666666666666666666666667e0;
  constexpr double Lob4i = 0.8333333333333333333333333e0;

  for (auto i = bo - 1; i < n; ++i) {
    i1 = i + 1;
    i1bo = (i1 - bo) * glq_pt;
    i2bo = i1bo * 2;
    dl = (kkn[i1] - kkn[i]) * 0.5;
    sl = (kkn[i1] + kkn[i]) * 0.5;

    sp = 0;
    for (auto p = 0; p < glq_pt; ++p, sp += 2) {
      r2 = dl * gl_x[p] + sl;
      dlob = (r2 - rm1) * 0.5;
      pr2 = wfn_o[p2i + p + i1bo] * wfn_o[p2p + p + i1bo];

      Qk += dlob *
            (fm1q +
             Lob4i * (1.0 / rk_mid[sp + i2bo + kp1 * nqpt2]) *
                 wfn_i[p2ii + sp + i2bo] * wfn_i[p2pi + sp + i2bo] +
             Lob4i * (1.0 / rk_mid[sp + 1 + i2bo + kp1 * nqpt2]) *
                 wfn_i[p2ii + sp + 1 + i2bo] * wfn_i[p2pi + sp + 1 + i2bo] +
             Lob4o * (1.0 / rk[p + i1bo + kp1 * nqpt]) * pr2);

      fm1q = Lob4o * (1.0 / rk[p + i1bo + kp1 * nqpt]) * pr2;
      rm1 = r2;
    }
  }

  rm1 = 0.0;
  fm1q = 0.0;
  double fm1j = 0.0;
  double pr1k, pr1km;

  for (auto i = bo - 1; i < n; ++i) {
    i1 = i + 1;
    i1bo = (i1 - bo) * glq_pt;
    i2bo = i1bo * 2;
    dl = (kkn[i + 1] - kkn[i]) * 0.5;
    sl = (kkn[i + 1] + kkn[i]) * 0.5;
    loc_GL = 0.0;

    sp = 0;
    for (auto p = 0; p < glq_pt; ++p, sp += 2) {
      r1 = dl * gl_x[p] + sl;
      dlob = (r1 - rm1) * 0.5;

      pr2 = wfn_o[p2i + p + i1bo] * wfn_o[p2p + p + i1bo];
      pr2a = wfn_i[p2ii + sp + i2bo] * wfn_i[p2pi + sp + i2bo];
      pr2b = wfn_i[p2ii + sp + 1 + i2bo] * wfn_i[p2pi + sp + 1 + i2bo];

      // chi(r1)
      pr1k = rk[p + i1bo + k * nqpt];
      pr1km = 1.0 / rk[p + i1bo + kp1 * nqpt];
      Jk += dlob * (fm1j + Lob4i * rk_mid[sp + i2bo + k * nqpt2] * pr2a +
                    Lob4i * rk_mid[sp + 1 + i2bo + k * nqpt2] * pr2b +
                    Lob4o * pr1k * pr2);
      Qk -= dlob *
            (fm1q + Lob4i * (1.0 / rk_mid[sp + i2bo + kp1 * nqpt2]) * pr2a +
             Lob4i * (1.0 / rk_mid[sp + 1 + i2bo + kp1 * nqpt2]) * pr2b +
             Lob4o * pr1km * pr2);
      chi = pr1km * Jk + pr1k * Qk;

      // Glq outer
      loc_GL += gl_w[p] * wfn_o[p1i + p + i1bo] * wfn_o[p1p + p + i1bo] * chi;

      rm1 = r1;
      fm1j = Lob4o * pr1k * pr2;
      fm1q = Lob4o * pr1km * pr2;
    }
    Fk += dl * loc_GL;
  }

  return Fk;
}

double FsltrLob3GL(int k, int n, int bo, int glq_pt, int lc_sz, int n1, int l1,
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

  constexpr double Lobo = 0.3333333333333333333333333e0;
  constexpr double Lobi = 0.1333333333333333333333333e1;

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

double FsltrTrapGL(int k, int n, int bo, int glq_pt, int lc_sz, int n1, int l1,
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
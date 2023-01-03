#include "integrator_new.hpp"

namespace glw {
const double GLobw2[2] = {1.0,1.0}; // Trapezoid rule
const double GLobw3[3] = {0.333333333333333333333333333333, 1.333333333333333333333333333333,
                          0.333333333333333333333333333333};
const double GLobw4[4] = {0.166666666666666666666666666667, 0.83333333333333333333333333333,
                          0.83333333333333333333333333333, 0.166666666666666666666666666667};
const double GLobw5[5] = {0.1, 0.544444444444444444444444444444, 0.71111111111111111111111111111,
                          0.544444444444444444444444444444, 0.1};
const double GLobw6[6] = {0.0666666666666666666666666666667, 0.378474956297846980316612808212,
                          0.554858377035486353016720525121, 0.554858377035486353016720525121,
                          0.378474956297846980316612808212, 0.0666666666666666666666666666667};
const double GLobw7[7] = {0.047619047619047619047619047619, 0.27682604736156594801070040629,
                          0.431745381209862623417871022281, 0.487619047619047619047619047619,
                          0.431745381209862623417871022281, 0.27682604736156594801070040629,
                          0.047619047619047619047619047619};
const double GLobw8[8] = {0.0357142857142857142857142857143, 0.21070422714350603938299206578,
                          0.34112269248350436476424067711, 0.4124587946587038815670529714,
                          0.4124587946587038815670529714, 0.34112269248350436476424067711,
                          0.21070422714350603938299206578, 0.0357142857142857142857142857143};
const double GLobw9[9] = {0.0277777777777777777777777777778, 0.16549536156080552504633972003,
                          0.274538712500161735280705618579, 0.34642851097304634511513153214,
                          0.371519274376417233560090702948, 0.34642851097304634511513153214,
                          0.274538712500161735280705618579, 0.16549536156080552504633972003,
                          0.0277777777777777777777777777778};
const double GLobw10[10] = {0.0222222222222222222222222222222, 0.133305990851070111126227170755,
                            0.224889342063126452119457821731,  0.292042683679683757875582257374,
                            0.327539761183897456656510527917, 0.327539761183897456656510527917,
                            0.292042683679683757875582257374, 0.224889342063126452119457821731,
                            0.133305990851070111126227170755, 0.0222222222222222222222222222222};

const double *Globw[9] = {GLobw2, GLobw3, GLobw4, GLobw5,
                           GLobw6, GLobw7, GLobw8, GLobw9, GLobw10};
} // namespace glw

double fsltrLob(int k, int qsz, int pti_sz, int lc_sz, int lci_sz,
                int n1, int l1, int n2, int l2, 
                int n3, int l3, int n4, int l4,
                std::vector<double> &q_w, std::vector<uint8_t> &pq_dx,
                std::vector<double> &rk, std::vector<double> &rk_in, 
                std::vector<double> &wfn_o, std::vector<double> &wfn_i) {
  int kp1 = k + 1;
  int p2i_o = l2 * lc_sz + n2 * qsz; // P(r2)
  int p2p_o = l4 * lc_sz + n4 * qsz; // P'(r2)
  int p1i_o = l1 * lc_sz + n1 * qsz; // P(r1)
  int p1p_o = l3 * lc_sz + n3 * qsz; // P'(r1)
  int p2i_i = l2 * lc_sz + n2 * pti_sz; // P(r2) inner points
  int p2p_i = l4 * lc_sz + n4 * pti_sz; // P'(r2)
  int p1i_i = l1 * lc_sz + n1 * pti_sz; // P(r1)
  int p1p_i = l3 * lc_sz + n3 * pti_sz; // P'(r1)
  int num_pti = 0;
  int in_idx = 0;
  double dl, sl, loc_GL, r1, pr2, pr2_in, chi;
  double rm1 = 0.0;
  double w_scale;
  double fm1j = 0.0, fm1q = 0.0;
  double pr1k, pr1km;

  // first calculate Qk for all of 0->R
  double Qk = 0.0, Qk_loc = 0.0;

  for (auto i = 0; i < qsz; ++i) {
    r1 = rk[qsz + i];
    num_pti = pq_dx[i];
    w_scale = (r1 - rm1) * 0.5;
    pr2 = (1.0 / rk[i + kp1 * qsz]) * wfn_o[p2i_o + i] * wfn_o[p2p_o + i];

    Qk_loc = glw::GLobw[num_pti][0] * fm1q;
    for (auto pt=0; pt<num_pti; ++pt) {
      Qk_loc += glw::GLobw[num_pti][pt] * (1.0 / rk_in[in_idx + kp1 * pti_sz])
              * wfn_i[p2i_i + in_idx] * wfn_i[p2p_i + in_idx];
      ++in_idx;
    }

    Qk_loc += glw::GLobw[num_pti][num_pti+1] * pr2;
    Qk += w_scale * Qk_loc;

    fm1q = pr2;
    rm1 = r1;
  }

  double Jk = 0.0, Jk_loc = 0.0;
  double Fk = 0.0;
  in_idx = 0;
  rm1 = 0.0;
  fm1q = 0.0;

  for (auto i = 0; i < qsz; ++i) {
    r1 = rk[qsz + i];
    num_pti = pq_dx[i];
    w_scale = (r1 - rm1) * 0.5;

    pr2 = wfn_o[p2i_o + i] * wfn_o[p2p_o + i];

    // chi(r1)
    pr1k = rk[i + k * qsz];
    pr1km = 1.0 / rk[i + kp1 * qsz];

    Jk_loc = glw::GLobw[num_pti][0] * fm1j;
    Qk_loc = glw::GLobw[num_pti][0] * fm1q;
    for (auto pt=0; pt<num_pti; ++pt) {
      pr2_in = wfn_i[p2i_i + in_idx] * wfn_i[p2p_i + in_idx];
      Jk_loc += glw::GLobw[num_pti][pt] * rk_in[in_idx + k * pti_sz] * pr2_in;
      Qk_loc += glw::GLobw[num_pti][pt] * (1.0 / rk_in[in_idx + kp1 * pti_sz]) * pr2_in;
      ++in_idx;
    }
    Jk_loc += glw::GLobw[num_pti][num_pti+1] * pr1k * pr2;
    Qk_loc += glw::GLobw[num_pti][num_pti+1] * pr1km * pr2;

    Jk += w_scale * Jk_loc;
    Qk -= w_scale * Qk_loc;

    chi = pr1km * Jk + pr1k * Qk;

    // Glq outer
    Fk += q_w[i] * wfn_o[p1i_o + i] * wfn_o[p1p_o + i] * chi;

    rm1 = r1;
    fm1j = pr1k * pr2;
    fm1q = pr1km * pr2;
  }
  
  return Fk;
}
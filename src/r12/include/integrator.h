#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>

// needs an update
double FsltrGL2(int k, int n, int bo, int gl2,
            std::vector<double> &gl_ow, 
            std::vector<double> &gl_ox,
            std::vector<double> &gl_iw, 
            std::vector<double> &gl_ix, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Gsp,
            std::vector<double> &Cl1i_pt,
            std::vector<double> &Cl1p_pt,
            std::vector<double> &Cl2i_pt,
            std::vector<double> &Cl2p_pt);

double FsltrMM(int k, int n, int bo, 
              int glq_pt, int lc_sz,
              int n1, int l1, int n2, int l2,
              int n3, int l3, int n4, int l4,
              std::vector<double> &gl_w, 
              std::vector<double> &gl_x, 
              std::vector<double> &kkn,
              std::vector<double> &rk,
              std::vector<double> &rk_mid,
              std::vector<double> &wfn_o,
              std::vector<double> &wfn_i);

double FsltrLob4GL(int k, int n, int bo, 
                  int glq_pt, int lc_sz,
                  int n1, int l1, int n2, int l2,
                  int n3, int l3, int n4, int l4,
                  std::vector<double> &gl_w, 
                  std::vector<double> &gl_x, 
                  std::vector<double> &kkn,
                  std::vector<double> &rk,
                  std::vector<double> &rk_mid,
                  std::vector<double> &wfn_o,
                  std::vector<double> &wfn_i);

double FsltrLob3GL(int k, int n, int bo, 
                  int glq_pt, int lc_sz,
                  int n1, int l1, int n2, int l2,
                  int n3, int l3, int n4, int l4,
                  std::vector<double> &gl_w, 
                  std::vector<double> &gl_x, 
                  std::vector<double> &kkn,
                  std::vector<double> &rk,
                  std::vector<double> &rk_mid,
                  std::vector<double> &wfn_o,
                  std::vector<double> &wfn_i);

double FsltrTrapGL(int k, int n, int bo,
                  int glq_pt, int lc_sz,
                  int n1, int l1, int n2, int l2,
                  int n3, int l3, int n4, int l4,
                  std::vector<double> &gl_w, 
                  std::vector<double> &gl_x, 
                  std::vector<double> &kkn,
                  std::vector<double> &rk,
                  std::vector<double> &wfn_o);

#endif // INTEGRATOR_H_
#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <vector>

double tvelGL(int n, int glq_pt, int bo, int lc_sz, int n1, int l1, int n2,
              int l2, std::vector<double> &gl_w, std::vector<double> &gl_x,
              std::vector<double> &kkn, std::vector<double> &wfn,
              std::vector<double> &wfnp);

double tlenGL(int n, int glq_pt, int bo, int lc_sz, int n1, int l1, int n2,
              int l2, std::vector<double> &gl_w, std::vector<double> &gl_x,
              std::vector<double> &kkn, std::vector<double> &wfn);

#endif // INTEGRATOR_H_
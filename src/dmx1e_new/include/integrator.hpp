#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <vector>

namespace dmx_int {
double tvelGL(int qsz, int lc_sz, int n1, int l1, int n2, int l2,
              std::vector<double> &qx, std::vector<double> &qw,
              std::vector<double> &wfn, std::vector<double> &wfnp);

double tlenGL(int qsz, int lc_sz, int n1, int l1, int n2, int l2,
              std::vector<double> &qx, std::vector<double> &qw,
              std::vector<double> &wfn);
} // namespace dmx_int

#endif // INTEGRATOR_H_
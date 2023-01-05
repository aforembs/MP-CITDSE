#ifndef INTEGRATOR_NEW_H_
#define INTEGRATOR_NEW_H_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// needs an update
namespace intfn {
double fsltrLob(int k, int qsz, int pti_sz, int lc_sz, int lci_sz, int n1,
                int l1, int n2, int l2, int n3, int l3, int n4, int l4,
                std::vector<double> &q_w, std::vector<uint8_t> &pq_dx,
                std::vector<double> &rk, std::vector<double> &rk_in,
                std::vector<double> &wfn_o, std::vector<double> &wfn_i);

} // namespace intfn

#endif // INTEGRATOR_H_
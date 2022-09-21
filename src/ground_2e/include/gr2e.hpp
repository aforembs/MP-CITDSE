#ifndef GR2E_HPP_
#define GR2E_HPP_

#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
extern "C" {
#include <cblas.h>
}

namespace gr2e {
int prop(std::string output, int L_sz, double t, double dt,
         std::vector<double> &ens, std::vector<double> &block,
         std::vector<double> &c0);
} // namespace gr2e

#endif // GR2E_HPP_
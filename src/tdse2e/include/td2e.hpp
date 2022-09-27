#ifndef TD2E_HPP_
#define TD2E_HPP_

#include "pulse.hpp"
#include <boost/numeric/odeint.hpp>
#include <complex>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>
extern "C" {
#include <cblas.h>
}

using stvupt = std::vector<std::unique_ptr<std::vector<double>>>;
using fieldInit = std::function<void(double, double, double, int, double &,
                                     double &, double &)>;
using fieldFcn = std::function<double(double, double, double, double, double)>;

namespace td2e {
int prop(std::string output, int L_max, double t, double dt, int steps,
         fieldInit fieldst, fieldFcn field, double w, double Io, double cepd,
         int cycles, int ct_sz, std::vector<int> &offs,
         std::vector<int> &state_sz, stvupt &blocks, stvupt &dipoles,
         std::vector<std::complex<double>> &ct);
} // namespace td2e

#endif // TD2E_HPP_
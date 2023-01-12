#ifndef TDSE_H_
#define TDSE_H_

#include <boost/numeric/odeint.hpp>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>
extern "C" {
#include <cblas.h>
}
#include "pulse.hpp"

using fieldInit = std::function<void(double, double, double, int, double &,
                                     double &, double &)>;
using fieldFcn = std::function<double(double, double, double, double, double)>;

namespace tdse {
int propV(std::string output, int l_max, double t, double dt, int steps,
          fieldInit fieldst, fieldFcn field, double Io, double w, double cepd,
          int cycles, int e_sz, std::vector<int> &offs, std::vector<int> &doffs,
          std::vector<int> &state_sz, std::vector<double *> &eps_off,
          std::vector<double *> &dip_off,
          std::vector<std::complex<double>> &ct);

int propL(std::string output, int l_max, double t, double dt, int steps,
          fieldInit fieldst, fieldFcn field, double Io, double w, double cepd,
          int cycles, int e_sz, std::vector<int> &offs, std::vector<int> &doffs,
          std::vector<int> &state_sz, std::vector<double *> &eps_off,
          std::vector<double *> &dip_off,
          std::vector<std::complex<double>> &ct);
} // namespace tdse

#endif // TDSE_H_
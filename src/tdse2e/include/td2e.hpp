#ifndef TD2E_HPP_
#define TD2E_HPP_

/**
 * @file tdse.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for propagating the 2-electron TDSE
 * @version 1.0
 */

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

/**
 * @brief Namespace for functions used in the propagation of the 2-electron TDSE
 */
namespace td2e {
/**
 * @brief Function for propagating the 2-electron TDSE in the velocity gauge
 *
 * @param output path prefix for the output files
 * @param L_max the maximum total angular momentum used
 * @param t the time at the start of the propagation
 * @param dt the time step
 * @param steps the total number of time steps
 * @param fieldst function pointer to the field initialisation fuction
 * @param field function pointer to the function used to generate the field at
 * time t
 * @param Io the peak intensity of the electric field in W/cm^2
 * @param w the central photon energy of the pulse in eV
 * @param cepd phase offset
 * @param cycles number of cycles in the pulse
 * @param ct_sz the size of the coefficient vector
 * @param offs vector of integers containing offsets to the start of each
 * angular momentum within the eigenenergy/coefficient vector
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each total angular momentum
 * @param eig a vector of pointers to the first eigenenergy in each
 * 2-electron total angular momentum
 * @param dipoles a vector of pointers to the first dipole element in each
 * 2-electron dipole transition matrix
 * @param ct vector containing the 2-electron wavefunction coefficients at time
 * t
 * @return int default '0' error otherwise
 */
int propV(std::string output, int L_max, double t, double dt, int steps,
          fieldInit fieldst, fieldFcn field, double w, double Io, double cepd,
          int cycles, int ct_sz, std::vector<int> &offs,
          std::vector<int> &state_sz, stvupt &eig, stvupt &dipoles,
          std::vector<std::complex<double>> &ct);

/**
 * @brief Function for propagating the 2-electron TDSE in the length gauge
 *
 * @param output path prefix for the output files
 * @param L_max the maximum total angular momentum used
 * @param t the time at the start of the propagation
 * @param dt the time step
 * @param steps the total number of time steps
 * @param fieldst function pointer to the field initialisation fuction
 * @param field function pointer to the function used to generate the field at
 * time t
 * @param Io the peak intensity of the electric field in W/cm^2
 * @param w the central photon energy of the pulse in eV
 * @param cepd phase offset
 * @param cycles number of cycles in the pulse
 * @param ct_sz the size of the coefficient vector
 * @param offs vector of integers containing offsets to the start of each
 * angular momentum within the eigenenergy/coefficient vector
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each total angular momentum
 * @param eig a vector of pointers to the first eigenenergy in each
 * 2-electron total angular momentum
 * @param dipoles a vector of pointers to the first dipole element in each
 * 2-electron dipole transition matrix
 * @param ct vector containing the 2-electron wavefunction coefficients at time
 * t
 * @return int default '0' error otherwise
 */
int propL(std::string output, int L_max, double t, double dt, int steps,
          fieldInit fieldst, fieldFcn field, double w, double Io, double cepd,
          int cycles, int ct_sz, std::vector<int> &offs,
          std::vector<int> &state_sz, stvupt &eig, stvupt &dipoles,
          std::vector<std::complex<double>> &ct);
} // namespace td2e

#endif // TD2E_HPP_
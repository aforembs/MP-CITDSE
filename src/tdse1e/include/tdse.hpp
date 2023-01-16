#ifndef TDSE_H_
#define TDSE_H_

/**
 * @file tdse.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for propagating the 1-electron TDSE
 * @version 1.0
 */

#include <boost/numeric/odeint.hpp>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>
extern "C" {
#include <cblas.h>
}
#include "pulse.hpp"

using stvupt = std::vector<std::unique_ptr<std::vector<double>>>;
using fieldInit = std::function<void(double, double, double, int, double &,
                                     double &, double &)>;
using fieldFcn = std::function<double(double, double, double, double, double)>;

/**
 * @brief Namespace for functions used in the propagation of the 1-electron TDSE
 */
namespace tdse {
/**
 * @brief Function for propagating the 1-electron TDSE in the velocity gauge
 *
 * @param output path prefix for the output files
 * @param l_max the maximum 1-electron angular momentum
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
 * @param e_sz the size of the coefficient vector
 * @param offs vector of integers containing offsets to the start of each
 * angular momentum within the eigenenergy/coefficient vector
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each angular momentum
 * @param eig a vector of pointers to the first eigenenergy in each
 * 1-electron angular momentum
 * @param dipoles a vector of pointers to the first dipole element in each
 * 1-electron dipole transition matrix
 * @param ct vector containing the 1-electron wavefunction coefficients at time
 * t
 * @return int default '0' error otherwise
 */
int propV(std::string output, int l_max, double t, double dt, int steps,
          fieldInit fieldst, fieldFcn field, double Io, double w, double cepd,
          int cycles, int e_sz, std::vector<int> &offs,
          std::vector<int> &state_sz, std::vector<double *> &eig,
          stvupt &dipoles, std::vector<std::complex<double>> &ct);

/**
 * @brief Function for propagating the 1-electron TDSE in the length gauge
 *
 * @param output path prefix for the output files
 * @param l_max the maximum 1-electron angular momentum
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
 * @param e_sz the size of the coefficient vector
 * @param offs vector of integers containing offsets to the start of each
 * angular momentum within the eigenenergy/coefficient vector
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each angular momentum
 * @param eig a vector of pointers to the first eigenenergy in each
 * 1-electron angular momentum
 * @param dipoles a vector of pointers to the first dipole element in each
 * 1-electron dipole transition matrix
 * @param ct vector containing the 1-electron wavefunction coefficients at time
 * t
 * @return int default '0' error otherwise
 */
int propL(std::string output, int l_max, double t, double dt, int steps,
          fieldInit fieldst, fieldFcn field, double Io, double w, double cepd,
          int cycles, int e_sz, std::vector<int> &offs,
          std::vector<int> &state_sz, std::vector<double *> &eig,
          stvupt &dipoles, std::vector<std::complex<double>> &ct);
} // namespace tdse

#endif // TDSE_H_
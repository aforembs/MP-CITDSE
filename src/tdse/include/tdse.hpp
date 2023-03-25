#ifndef TDSE_HPP_
#define TDSE_HPP_

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

/**
 * @brief Aliases for vectors of unique_ptr and the pulse function pointer.
 *
 */
using stvupt = std::vector<std::unique_ptr<std::vector<double>>>;
using fieldFcn = std::function<double(pulse::params &, double)>;

/**
 * @brief Namespace for functions used in the propagation of the 2-electron TDSE
 */
namespace tdse {

/**
 * @brief Function for propagating the 2-electron TDSE in the velocity gauge
 *
 * @param output path prefix for the output files
 * @param L_max the maximum total angular momentum used
 * @param t the time at the start of the propagation
 * @param dt the time step
 * @param steps the total number of time steps
 * @param pop_n primary quantum number of the state whose population will be
 * tracked
 * @param pop_l angular quantum number of the state
 * @param field function pointer to the function used to generate the field at
 * time t
 * @param pars structure containing static pulse parameters
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
          int pop_n, int pop_l, fieldFcn field, pulse::params &pars, int ct_sz,
          std::vector<int> &offs, std::vector<int> &state_sz, stvupt &eig,
          stvupt &dipoles, std::vector<double> &ct);

/**
 * @brief Function for propagating the 2-electron TDSE in the length gauge
 *
 * @param output path prefix for the output files
 * @param L_max the maximum total angular momentum used
 * @param t the time at the start of the propagation
 * @param dt the time step
 * @param steps the total number of time steps
 * @param pop_n primary quantum number of the state whose population will be
 * tracked
 * @param pop_l angular quantum number of the state
 * @param field function pointer to the function used to generate the field at
 * time t
 * @param pars structure containing static pulse parameters
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
          int pop_n, int pop_l, fieldFcn field, pulse::params &pars, int ct_sz,
          std::vector<int> &offs, std::vector<int> &state_sz, stvupt &eig,
          stvupt &dipoles, std::vector<double> &ct);
} // namespace tdse

#endif // TDSE_HPP_
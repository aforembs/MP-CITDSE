#ifndef TD_READ_H_
#define TD_READ_H_

/**
 * @file td_read.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for reading the input file and components of the
 * 1-electron Hamiltonian.
 * @version 1.0
 */

#include <H5Cpp.h>
#include <cassert>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <vector>
#include <yaml-cpp/yaml.h>

using stvupt = std::vector<std::unique_ptr<std::vector<double>>>;

/**
 * @brief Namespace for functions used to read the components of the 1-electron
 * time-dependent Hamiltonian.
 */
namespace tdrd {
/**
 * @brief Function for reading settings from a .yaml file
 *
 * @param file path to the .yaml file
 * @param pot the name of the selected potential
 * @param gauge the gauge of the dipole elements 'v' (velocity) / 'l' (length)
 * @param l_max the maximum 1-electron angular momentum used
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each angular momentum
 * @param timestep the dt parameter used in the time propagation
 * @param w the central photon energy of the pulse in eV
 * @param Io the peak intensity of the electric field in W/cm^2
 * @param cepd phase offset
 * @param cycles number of cycles in the pulse
 * @return int default '0' error otherwise
 */
int readConfig(std::string file, std::string &pot, char &gauge, int &l_max,
               std::vector<int> &state_sz, double &timestep, double &w,
               double &Io, double &cepd, int &cycles);

/**
 * @brief Function for reading the eigenenergies of the time independent
 * 1-electron Hamiltonian.
 *
 * @param pot the name of the selected potential
 * @param l_max the maximum 1-electron angular momentum used
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each angular momentum
 * @param eps vector containing the eigenenergies of the 1-electron Hamiltonian
 * @param offs vector of integers containing offsets to the start of each
 * angular momentum within the eigenenergy/coefficient vector
 * @param eps_sz the size of the eigenenergy vector
 * @param eig a vector of pointers to the first eigenenergy in each
 * 1-electron angular momentum
 * @return int default '0' error otherwise
 */
int readEnergies(std::string pot, int l_max, std::vector<int> &state_sz,
                 std::vector<double> &eps, std::vector<int> &offs, int &eps_sz,
                 std::vector<double *> &eig);

/**
 * @brief Function for reading the 1-electron dipole transition matrices.
 *
 * @param pot the name of the selected potential
 * @param gauge the gauge of the dipole elements 'v' (velocity) / 'l' (length)
 * @param l_max the maximum 1-electron angular momentum used
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each angular momentum
 * @param dipoles a vector of pointers to the first dipole element in each
 * 1-electron dipole transition matrix
 * @return int default '0' error otherwise
 */
int readDipoles(std::string pot, char gauge, int l_max,
                std::vector<int> &state_sz, stvupt &dipoles);
} // namespace tdrd

#endif // TD_READ_H_
#ifndef TD_READ_HPP_
#define TD_READ_HPP_

/**
 * @file td_read.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for reading the input file and components of the
 * 2-electron Hamiltonian.
 * @version 1.0
 */

#include "dmx_typ.hpp"
#include <H5Cpp.h>
#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>
#include <yaml-cpp/yaml.h>
extern "C" {
#include <cblas.h>
}

using stvupt = std::vector<std::unique_ptr<std::vector<double>>>;

/**
 * @brief Namespace for functions used to read the components of the 2-electron
 * time-dependent Hamiltonian.
 */
namespace tdrd {
/**
 * @brief Function for reading settings from a .yaml file
 *
 * @param file path to the .yaml file
 * @param pot the name of the selected potential
 * @param set_base sub-setting group for L_max, different for 1e and 2e cases
 * @param option setting name for L_max, different for 1e and 2e cases
 * @param L_max the maximum total angular momentum used
 * @param gauge the gauge of the dipole elements 'v' (velocity) / 'l' (length)
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each total angular momentum
 * @param timestep the dt parameter used in the time propagation
 * @param w the central photon energy of the pulse in eV
 * @param Io the peak intensity of the electric field in W/cm^2
 * @param cepd phase offset
 * @param cycles number of cycles in the pulse
 * @return int default '0' error otherwise
 */
int readConfig(std::string file, std::string &pot, std::string set_base,
               std::string option, int &L_max, char &gauge,
               std::vector<int> &state_sz, double &timestep, double &w,
               double &Io, double &cepd, int &cycles);

/**
 * @brief Function for reading the eigenenergies of the Configuration
 * Interaction form of the time independent 2-electron Hamiltonian.
 *
 * @param pot the name of the selected potential
 * @param L_max the maximum total angular momentum used
 * @param ct_sz the size of the eigenenergy/coefficient vector
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each total angular momentum
 * @param offs vector of integers containing offsets to the start of each
 * angular momentum within the eigenenergy/coefficient vector
 * @param ens a vector of pointers to the first eigenenergy in each
 * 2-electron total angular momentum
 * @return int default '0' error otherwise
 */
int readEnergies(std::string pot, std::string setname, int L_max, int &ct_sz,
                 std::vector<int> &state_sz, std::vector<int> &offs,
                 stvupt &eig);

/**
 * @brief Function for reading the 2-electron dipole transition matrices.
 *
 * @param pot the name of the selected potential
 * @param gauge the gauge of the dipole elements 'v' (velocity) / 'l' (length)
 * @param L_max the maximum total angular momentum used
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each total angular momentum
 * @param dipoles a vector of pointers to the first dipole element in each
 * 2-electron dipole transition matrix
 * @return int default '0' error otherwise
 */
int readDipoles(std::string pot, std::string setname, char gauge, int L_max,
                std::vector<int> &state_sz, stvupt &dipoles);

/**
 * @brief
 *
 * @param file
 * @param ct_sz
 * @param ct
 * @return int
 */
int readInitCt(std::string file, int ct_sz,
               std::vector<std::complex<double>> &ct);
} // namespace tdrd

#endif // TD_READ_HPP_
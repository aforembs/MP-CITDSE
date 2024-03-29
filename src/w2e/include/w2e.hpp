#ifndef W2E_HPP_
#define W2E_HPP_

/**
 * @file w2e.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for generating the Configuration Interaction 2-electron
 * basis.
 * @version 1.0
 */

#include "dmx_typ.hpp"
#include <H5Cpp.h>
#include <cassert>
#include <iostream>
#ifndef __INTEL_MKL__
#include <lapacke.h>
extern "C" {
#include <cblas.h>
}
#else 
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#endif
#include <numeric>
#include <vector>
#include <yaml-cpp/yaml.h>

using stvupt = std::vector<std::unique_ptr<std::vector<double>>>;

/**
 * @brief Namespace for CI basis functions
 */
namespace w2e {
/**
 * @brief Function for reading settings from a .yaml file
 *
 * @param file path to the .yaml file
 * @param pot the name of the selected potential
 * @param gauge the gauge of the dipole elements 'v' (velocity) / 'l' (length)
 * @param L_max the maximum total angular momentum used
 * @return int default '0' error otherwise
 */
int readConfig(std::string file, std::string &pot, char &gauge, int &L_max);

/**
 * @brief Function for calculating and saving the eigenergies and eigenvectors
 * of the CI basis.
 *
 * @param pot the name of the selected potential, used as a file prefix for
 * outputs
 * @param L_max the maximum total angular momentum used
 * @param vecs the eigenvectors of the CI basis
 * @return int default '0' error otherwise
 */
int formCIh0(std::string pot, int L_max, stvupt &vecs);

/**
 * @brief Function for projecting the dipole matrix elements onto the CI basis
 *
 * @param pot the name of the selected potential, used as a file prefix for
 * outputs
 * @param gauge the gauge of the dipole elements 'v' (velocity) / 'l' (length)
 * @param L_max the maximum total angular momentum used
 * @param vecs the eigenvectors of the CI basis
 * @return int default '0' error otherwise
 */
int formCIDipoles(std::string pot, char gauge, int L_max, stvupt &vecs);
} // namespace w2e

#endif // W2E_HPP_

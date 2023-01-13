#ifndef DMX1E_H_
#define DMX1E_H_

/**
 * @file dmx1e.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for calculating 1-electron dipole matrices
 * @version 1.0
 */

#include "dmx_typ.hpp"
#include "fastgl.hpp"
#include "integrator.hpp"
#include <H5Cpp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <yaml-cpp/yaml.h>

/**
 * @brief Namespace for the 1-electron dipole matrix funcitons
 */
namespace dmx1e {
/**
 * @brief Function for reading settings from a .yaml file
 *
 * @param file path to the .yaml file
 * @param pot the name of the selected potential
 * @param qsz the size of the quadrature used for integration
 * @param gauge the gauge of the dipole elements 'v' (velocity) / 'l' (length)
 * @param l_max the maximum 1-electron angular momentum used
 * @return int default '0' error otherwise
 */
int readConfig(std::string file, std::string &pot, int &qsz, char &gauge,
               int &l_max);

/**
 * @brief Function for generating and saving the 1-electron dipole matrices
 *
 * @param pot the name of the selected potential
 * @param qsz the size of the quadrature used for integration
 * @param gauge the gauge of the dipole elements 'v' (velocity) / 'l' (length)
 * @param l_max the maximum 1-electron angular momentum used
 * @return int default '0' error otherwise
 */
int genDipole(std::string pot, int qsz, char gauge, int l_max);
} // namespace dmx1e

#endif // DMX1E_H_
#ifndef DMX2E_H_
#define DMX2E_H_

/**
 * @file dmx2e.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for calculating 2-electron dipole matrices
 * @version 1.0
 */

#include "cfg_in.hpp"
#include "dmx_typ.hpp"
#include "wigxjpf.h"
#include <H5Cpp.h>
#include <algorithm>
#include <cmath>
#include <execution>
#include <fstream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

/**
 * @brief Namespace for the 2-electron dipole matrix funcitons
 */
namespace dmx2e {
/**
 * @brief Function for reading settings from a .yaml file
 *
 * @param file path to the .yaml file
 * @param pot the name of the selected potential
 * @param L_max the maximum total angular momentum used
 * @param l_max the maximum 1-electron angular momentum used
 * @param gauge the gauge of the dipole elements 'v' (velocity) / 'l' (length)
 * @return int default '0' error otherwise
 */
int readConfig(std::string file, std::string &pot, int &L_max, int &l_max,
               char &gauge);

/**
 * @brief Function for calculating the 2-electron dipole matrix elements
 *
 * @param pot the name of the selected potential
 * @param L_max the maximum total angular momentum used
 * @param l_max the maximum 1-electron angular momentum used
 * @param gauge the gauge of the dipole elements 'v' (velocity) / 'l' (length)
 * @param dir the directory containing the cfg-<L>.inp files which define what
 * 2-electron configuraitons are included
 * @return int default '0' error otherwise
 */
int genDipole(std::string pot, int L_max, int l_max, char gauge,
              std::string dir);
} // namespace dmx2e

#endif // DMX2E_H_
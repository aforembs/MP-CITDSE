#ifndef GENIDX_HPP_
#define GENIDX_HPP_

/**
 * @file genidx.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for generating indices for the 2-electron
 * configurations
 * @version 1.0
 */

#include "cfg_in.hpp"
#include "dmx_typ.hpp"
#include <H5Cpp.h>
#include <iostream>
#include <memory>
#include <yaml-cpp/yaml.h>

/**
 * @brief Namespace for functions used to generate indices for 2-electron
 * configurations
 */
namespace genidx {
/**
 * @brief Function for reading settings from a .yaml file
 *
 * @param file path to the .yaml file
 * @param pot the name of the selected potential
 * @param L_max the maximum total angular momentum used
 * @return int default '0' error otherwise
 */
int readConfig(std::string file, std::string &pot, int &L_max);

/**
 * @brief Function for saving the primary and angular momentum quantum numbers
 * of the 1-electron states used in each configuration
 *
 * @param pot the name of the selected potential
 * @param L_max the maximum total angular momentum used
 * @param dir the directory containing the cfg-<L>.inp files which define what
 * 2-electron configuraitons are included
 * @return int default '0' error otherwise
 */
int saveIdx(std::string pot, int L_max, std::string dir);
} // namespace genidx

#endif // GENIDX_HPP_
#ifndef PES_HPP_
#define PES_HPP_

/**
 * @file pes.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for calculating the photoelectron energy distribution
 * from 2-electron TDSE coefficients.
 * @version 1.0
 */

#include "cfg_in.hpp"
#include "dmx_typ.hpp"
#include <H5Cpp.h>
#include <cassert>
#include <complex>
#include <cstring>
#include <iostream>
#ifndef __INTEL_MKL__
#include <lapacke.h>
#endif
#include <regex>
#include <vector>
#include <yaml-cpp/yaml.h>
extern "C" {
#include <cblas.h>
}

/**
 * @brief Namespace for functions used to generate a photoelectron energy
 * spectrum from 2-electron coefficents.
 */
namespace pes {
/**
 * @brief Function for reading settings from a .yaml file
 *
 * @param file path to the .yaml file
 * @param pot the name of the selected potential
 * @param set_base sub-setting group for L_max, different for 1e and 2e cases
 * @param option setting name for L_max, different for 1e and 2e cases
 * @param L_max the maximum total angular momentum used
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each total angular momentum 'L'
 * @return int default '0' error otherwise
 */
int readConfig(std::string file, std::string &pot, std::string set_base,
               std::string option, int &L_max, std::vector<int> &state_sz);

/**
 * @brief Read the 2-electron state coefficients from a file
 *
 * @param file path to the text file containing the 1-electron state
 * coefficients
 * @param ct vector containing the 1-electron state coefficients
 * @return int default '0' error otherwise
 */
int readCt(std::string file, std::vector<std::complex<double>> &ct);

/**
 * @brief Generate the photoelectron energy spectrum from the 1-electron state
 * coefficients
 *
 * @param pot the name of the selected potential
 * @param s_flag
 * @param l_max the maximum 1-electron angular momentum used
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each angular momentum
 * @param ct vector containing the 1-electron state coefficients
 * @param output path to directory where the PES file will be saved
 * @return int default '0' error otherwise
 */
int genPES1e(std::string pot, bool s_flag, int l_max,
             std::vector<int> &state_sz, std::vector<std::complex<double>> &ct,
             std::string output);

/**
 * @brief Generate the photoelectron energy spectrum from the 2-electron
 * state coefficients
 *
 * @param pot the name of the selected potential
 * @param s_flag
 * @param L_max the maximum total angular momentum used
 * @param state_sz vector containing the sizes of the coefficient vectors
 * for each angular momentum
 * @param ct vector containing the 1-electron state coefficients
 * @param output path to directory where the PES file will be saved
 * @return int default '0' error otherwise
 */
int genPES2e(std::string pot, bool s_flag, int L_max,
             std::vector<int> &state_sz, std::vector<std::complex<double>> &ct,
             std::string output);
} // namespace pes

#endif // PES_HPP_
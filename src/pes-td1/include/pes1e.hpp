#ifndef PES_H_
#define PES_H_

/**
 * @file pes1e.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for calculating the photoelectron energy distribution
 * from 1-electron TDSE coefficients.
 * @version 1.0
 */

#include <H5Cpp.h>
#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <yaml-cpp/yaml.h>

/**
 * @brief Namespace for functions used to generate a photoelectron energy
 * spectrum from 1-electron coefficents.
 */
namespace pes {
/**
 * @brief Function for reading settings from a .yaml file
 *
 * @param file path to the .yaml file
 * @param pot the name of the selected potential
 * @param l_max the maximum 1-electron angular momentum used
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each angular momentum
 * @return int default '0' error otherwise
 */
int readConfig(std::string file, std::string &pot, int &l_max,
               std::vector<int> &state_sz);

/**
 * @brief Read the 1-electron state coefficients from a file
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
 * @param l_max the maximum 1-electron angular momentum used
 * @param state_sz vector containing the sizes of the coefficient vectors for
 * each angular momentum
 * @param ct vector containing the 1-electron state coefficients
 * @param output path to directory where the PES file will be saved
 * @return int default '0' error otherwise
 */
int genPES(std::string pot, int l_max, std::vector<int> &state_sz,
           std::vector<std::complex<double>> &ct, std::string output);
} // namespace pes

#endif // PES_H_
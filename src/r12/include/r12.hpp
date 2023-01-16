#ifndef V_12_H_
#define V_12_H_

/**
 * @file r12.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for calculating the photoelectron energy distribution
 * from 1-electron TDSE coefficients.
 * @version 1.0
 */

#include "cfg_in.hpp"
#include "dmx_typ.hpp"
#include "fastgl.hpp"
#include "integrator.hpp"
#include "wigxjpf.h"
#include <H5Cpp.h>
#include <algorithm>
#include <iomanip>
#include <omp.h>
#include <vector>
#include <yaml-cpp/yaml.h>

/**
 * @brief Namespace for functions used in calculating the <12|r_{12}|1'2'> two
 * electron correlations.
 */
namespace r_12 {
/**
 * @brief Function for reading settings from a .yaml file
 *
 * @param file path to the .yaml file
 * @param qsz the size of the quadrature used for integration
 * @param pot the name of the selected potential
 * @param L_max the maximum total angular momentum used
 * @param k_limit parameter used to define the k_max limit of the multipole
 * series
 * @param lim_flag flag indicating if a numeric limit is set if false k_max is
 * defined by the Wigner 6j sybol
 * @return int default '0' error otherwise
 */
int readConfig(std::string file, int &qsz, std::string &pot, int &L_max,
               std::string &k_limit, bool &lim_flag);

/**
 * @brief Generate and save the <12|r_{12}|1'2'> two electron correlations for
 * each of the total angular momenta 'L' up to and including L_max.
 *
 * @param pot the name of the selected potential
 * @param L_max the maximum total angular momentum used
 * @param qsz the size of the quadrature used for integration
 * @param dir path to directory containing the cfg-<L>.inp files
 * @param lim_flag flag indicating if a numeric limit is set if false k_max is
 * defined by the Wigner 6j sybol
 * @param k_limit parameter used to define the k_max limit of the multipole
 * series
 * @return int default '0' error otherwise
 */
int r12Glob(std::string pot, int L_max, int qsz, std::string dir, bool lim_flag,
            std::string k_limit);

} // namespace r_12

#endif // V_12_H_
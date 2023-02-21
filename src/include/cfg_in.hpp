#ifndef CFG_IN_H_
#define CFG_IN_H_

/**
 * @file cfg_in.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for calculating the B-spline basis.
 * @version 1.0
 */

#include "dmx_typ.hpp"
#include <H5Cpp.h>
#include <algorithm>
#include <execution>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Namespace for functions used in reading the "cfg-<L>.inp" files
 */
namespace cfg {
/**
 * @brief Structure for holding elemets of a single line of a cfg file
 *
 * @param n1 primary quantum number of the 1st electron
 * @param l1 angular quantum number of the 1st electron
 * @param l2 angular quantum number of the 2nd electron
 * @param n2min minimum primary quantum number of the 2nd electron
 * @param n2max maximum primary quantum number of the 2nd electron
 */
struct line {
  int n1;
  int l1;
  int l2;
  int n2min;
  int n2max;
};

/**
 * @brief Function for reading the "cfg-<L>.inp" files and saving the lines into
 * a vector of cfg::line structures
 *
 * @param dir path to directory containing the "cfg-<L>.inp" files
 * @param L the total angular momentum
 * @param sym symmetry '1'/'3' singlet or triplet (for now only singlet)
 * @param ncf number of lines in the cfg file
 * @param cfgs vector containing lines from a cfg file
 * @return int default '0' error otherwise
 */
int readCfg(std::string dir, int L, int &sym, int &ncf,
            std::vector<line> &cfgs);
} // namespace cfg

#endif // CFG_IN_H_
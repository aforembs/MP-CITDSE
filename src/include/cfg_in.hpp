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
 * @brief
 */
namespace cfg {
struct line {
  int n1;
  int l1;
  int l2;
  int n2min;
  int n2max;
};

/**
 * @brief
 *
 * @param dir
 * @param L
 * @param sym
 * @param ncf
 * @param cfgs
 * @return int
 */
int readCfg(std::string dir, int L, int &sym, int &ncf,
            std::vector<line> &cfgs);

/**
 * @brief
 *
 * @param pot
 * @param gauge
 * @param L_max
 * @param dir
 * @return int
 */
int genL_idx(std::string pot, char gauge, int L_max, std::string dir);
} // namespace cfg

#endif // CFG_IN_H_
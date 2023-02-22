#ifndef DMX_TYP_H_
#define DMX_TYP_H_

/**
 * @file dmx_typ.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Structures for storing 2-electron configurations
 * @version 1.0
 */

#include <vector>

/**
 * @brief Structures for storing the 1-electron information for each
 * configuration
 *
 */
namespace dmtp {

/**
 * @brief Structure containing the uncorrelated energy of the direct product of
 * 2 1-electron states as well as their primary and angular momentum quantum
 * numbers.
 */
struct en_data {
  double en;
  int n1;
  int l1;
  int n2;
  int l2;
};

/**
 * @brief The four integers defining the primary and angular quantum numbers of
 * a 2-electron configuration.
 *
 */
struct idx4 {
  int n1;
  int l1;
  int n2;
  int l2;
};

} // namespace dmtp

#endif // DMX_TYP_H_
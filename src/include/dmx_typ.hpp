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
 * @brief
 *
 */
namespace dmtp {

/**
 * @brief
 *
 */
struct l_ab {
  int l1;
  int l2;
};

/**
 * @brief
 *
 */
struct en_data {
  double en;
  int n1;
  int l1;
  int n2;
  int l2;
};

/**
 * @brief
 *
 */
struct idx4 {
  int n1;
  int l1;
  int n2;
  int l2;
};

/**
 * @brief
 *
 */
struct en_L {
  int L;
  std::vector<l_ab> l_pair;
  std::vector<en_data> en_dat;
};

/**
 * @brief
 *
 */
struct dmx_dim {
  int row;
  int col;
};
} // namespace dmtp

#endif // DMX_TYP_H_
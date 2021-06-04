#ifndef DMX_TYP_H_
#define DMX_TYP_H_

#include <vector>

struct l_ab {
  uint l1;
  uint l2;
};

struct en_data {
  double en;
  uint n1;
  uint l1;
  uint n2;
  uint l2;
};

struct idx4 {
  uint n1;
  uint l1;
  uint n2;
  uint l2;
};

struct en_L {
  uint L;
  std::vector<l_ab> l_pair;
  std::vector<en_data> en_dat;
};

#endif // DMX_TYP_H_
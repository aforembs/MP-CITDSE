#ifndef DMX_TYP_H_
#define DMX_TYP_H_

#include <vector>

struct l_ab {
  int l1;
  int l2;
};

struct en_data {
  double en;
  int n1;
  int l1;
  int n2;
  int l2;
};

struct idx4 {
  int n1;
  int l1;
  int n2;
  int l2;
};

struct en_L {
  int L;
  std::vector<l_ab> l_pair;
  std::vector<en_data> en_dat;
};

struct dmx_dim {
  int row;
  int col;
};

#endif // DMX_TYP_H_
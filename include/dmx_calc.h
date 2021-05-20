#ifndef DMX_CALC_H_
#define DMX_CALC_H_

#include <cmath>
#include <vector>
#include <string>
#include <execution>
#include <algorithm>
#include "H5Cpp.h"
#include "wigner6j.h"

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

struct en_L {
  uint L;
  uint max_n;
  std::vector<l_ab> l_pair;
  std::vector<en_data> en_dat;
};

// Create a 2eDMX class!
class DMX2e {
  private:
    /* Class Variables */
    char gauge;
    std::string pot;

    en_L make_enL(uint L);
    int sort_L(en_L &Lif, uint N_sz);

  public:
    DMX2e(std::string cpot, char gau, uint L_max, uint N_max);


};

// int calculate_2edmx(uint L_max, uint Ni_sz, uint Nf_sz);

#endif // DMX_CALC_H_
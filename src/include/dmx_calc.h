#ifndef DMX_CALC_H_
#define DMX_CALC_H_

#include <cmath>
#include <vector>
#include <string>
#include <execution>
#include <algorithm>
#include <fstream>
#include "H5Cpp.h"
#include "dmx_typ.h"
#include "wigner6j.h"

class DMX2e {
  private:
    /* Class Variables */
    char gauge;
    std::string pot;

    en_L make_enL(uint L, uint l1e_max);
    int sort_L(uint L_max, std::vector<uint> &N_sz);

  public:
    DMX2e(std::string cpot, char gau, uint L_max, std::vector<uint> &N_max);


};

// int calculate_2edmx(uint L_max, uint Ni_sz, uint Nf_sz);

#endif // DMX_CALC_H_
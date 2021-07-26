#ifndef DMX_CALC_WIGNER_SYM_H_
#define DMX_CALC_WIGNER_SYM_H_

#include <cmath>
#include <algorithm>
#include <iostream>

double wigner_3j0(uint j1, uint j2, uint j3);

double wigner_6j(uint j1, uint j2, uint j3, 
                uint j4, uint j5, uint j6);

double wigner_6j_2e(uint L, uint la, uint lb, uint lc);

#endif // DMX_CALC_WIGNER_SYM_H_
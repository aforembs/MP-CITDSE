#ifndef DMX_CALC_WIGNER_SYM_H_
#define DMX_CALC_WIGNER_SYM_H_

#include <cmath>
#include <algorithm>
#include <iostream>

double wigner_3j0(unsigned int j1, unsigned int j2, unsigned int j3);

double wigner_6j(unsigned int j1, unsigned int j2, unsigned int j3, 
                unsigned int j4, unsigned int j5, unsigned int j6);

double wigner_6j_2e(unsigned int L, unsigned int la, unsigned int lb, unsigned int lc);

#endif // DMX_CALC_WIGNER_SYM_H_
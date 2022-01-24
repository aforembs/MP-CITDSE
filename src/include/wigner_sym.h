#ifndef DMX_CALC_WIGNER_SYM_H_
#define DMX_CALC_WIGNER_SYM_H_

#include <cmath>
#include <algorithm>
#include <iostream>

double wigner_3j0(int j1, int j2, int j3);

double wigner_6j(int j1, int j2, int j3, 
                int j4, int j5, int j6);

double wigner_6j_2e(int L, int la, int lb, int lc);

#endif // DMX_CALC_WIGNER_SYM_H_
#include <iostream>
#include "V_12.h"

int main() {

  std::vector<uint> N_maxs = {2,2};

  // output file prefix, L_max, vector of numbers of states per L
  V12("dat/he", 1, N_maxs);

  return 0;
}
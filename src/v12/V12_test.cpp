#include <iostream>
#include "V_12.h"

int main() {

  std::vector<uint> N_maxs = {20, 20, 10, 10}; // has to match dmx input

  // output file prefix, L_max, vector of numbers of states per L
  V12("dat/he", 3, N_maxs);

  return 0;
}
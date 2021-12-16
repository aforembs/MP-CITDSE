#include <iostream>
#include "V_12.h"

int main() {

  std::vector<uint> N_maxs = {100, 100, 50, 50}; // has to match dmx input

  // output file prefix, L_max, vector of numbers of states per L
  V12("dat/he", 3, "inp");

  return 0;
}
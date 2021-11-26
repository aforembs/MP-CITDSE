#include <iostream>
#include "V_12.h"

int main() {

  std::vector<uint> N_maxs = {100, 100, 50, 50}; // has to match dmx input

  // output file prefix, L_max, vector of numbers of states per L
  V12_alt("dat/he", 0, "inp");

  return 0;
}
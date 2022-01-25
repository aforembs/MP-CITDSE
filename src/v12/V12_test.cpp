#include <iostream>
#include "r12.h"

int main() {

  std::vector<uint> N_maxs = {100, 100, 50, 50}; // has to match dmx input

  // output file prefix, L_max, vector of numbers of states per L
  r_12::R12("dat/he", 2, 9, "inp");

  return 0;
}
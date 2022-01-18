#include <iostream>
#include "dmx2e.h"

int main() {

  // Vector of the maximum number of states from (1e) l= 0,1,2, ...
  // std::vector<uint> N_maxs = {100, 100, 50, 50};

  // DMX2e DMXcalculate("dat/he", 'v', 3, N_maxs);

  DMX2e DMXcalc("dat/he", 'v', 2, 5, "inp");

  return 0;
}
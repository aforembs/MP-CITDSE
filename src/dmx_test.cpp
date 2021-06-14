#include <iostream>
#include "dmx_calc.h"

int main() {

  // Vector of the maximum number of states from (1e) l= 0,1,2, ...
  std::vector<uint> N_maxs = {100, 100, 60, 60};

  DMX2e DMXcalculate("he", 'v', 3, N_maxs);

  return 0;
}
#include <iostream>
#include "dmx_calc.h"

int main() {

  std::vector<uint> N_maxs = {100, 100, 60, 60};

  DMX2e DMXcalculate("he", 'v', 3, N_maxs);

  return 0;
}
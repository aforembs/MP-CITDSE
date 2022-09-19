#ifndef PES2E_HPP_
#define PES2E_HPP_

#include "dmx_typ.h"
#include <cassert>
#include <complex>
#include <iostream>
#include <vector>

namespace pes2e {
int readConfig();

int readCt();

int genPES();
} // namespace pes2e

#endif // PES2E_HPP_
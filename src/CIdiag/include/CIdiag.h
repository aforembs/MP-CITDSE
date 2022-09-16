#ifndef CIDIAG_H_
#define CIDIAG_H_

#include "dmx_typ.h"
#include <H5Cpp.h>
#include <fstream>
#include <iostream>
#include <lapacke.h>
#include <memory>
#include <vector>

int CalcCI(std::string pot, char gauge, int L_max);

#endif // CIDIAG_H_
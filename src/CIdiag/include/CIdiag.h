#ifndef CIDIAG_H_
#define CIDIAG_H_

#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <lapacke.h>
#include "H5Cpp.h"
#include "dmx_typ.h"

int CalcCI(std::string pot, char gauge, int L_max);

#endif // CIDIAG_H_
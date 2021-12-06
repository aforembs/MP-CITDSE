#ifndef CIDIAG_H_
#define CIDIAG_H_

#include <vector>
#include <memory>
#include <iostream>
#include <lapacke.h>
#include "H5Cpp.h"

int CalcCI(std::string pot, char gauge, int L_max);

#endif // CIDIAG_H_
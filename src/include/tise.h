#ifndef TISE_H_
#define TISE_H_

#include <vector>
#include <algorithm>
#include "ModelV.h"
extern "C" {
  #include <lapacke.h>
}

namespace tise {

int GenCoeff();

}

#endif // TISE_H_
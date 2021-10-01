#ifndef TISE_H_
#define TISE_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <execution>
#include "ModelV.h"
#include "bsplines.h"
#include "H5Cpp.h"
#include <lapacke.h>


namespace tise {

int GenCoeff(int n, int k, int l,
            double z, double mass,
            std::string pot,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x,
            std::vector<double> &kkn,
            std::vector<double> &spl,
            std::vector<double> &splp,
            std::string outFile);

}

#endif // TISE_H_
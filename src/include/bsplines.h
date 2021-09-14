#ifndef BSPLINES_H_
#define BSPLINES_H_

#include <vector>
#include "slatec_f.h"

int bsplines(int n, int k, 
             std::vector<double> &gl_x,
             std::vector<double> &knots,
             std::vector<double> &splines);

int bsplinesp(int n, int k, 
              std::vector<double> &gl_x,
              std::vector<double> &knots,
              std::vector<double> &splinesp);

#endif // BSPLINES_H_
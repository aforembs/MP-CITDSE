#ifndef BSPLINES_H_
#define BSPLINES_H_

#include <cmath>
#include <vector>
#include <ifstream>
#include "slatec_f.h"

namespace bsp {

int GenKnots(int n, int k, double r_max, 
             double fkn, char type,
             std::vector<double> &kkn);

int GenKnots(int n, int k, 
             std::string file, char type, 
             std::vector<double> &kkn);

int WrKnotsH5();

int Splines(int n, int k, 
             std::vector<double> &gl_x,
             std::vector<double> &knots,
             std::vector<double> &splines);

int SplinesP(int n, int k, 
              std::vector<double> &gl_x,
              std::vector<double> &knots,
              std::vector<double> &splinesp);

int SplineInt(int n, int k,
              std::vector<double> &gl_w,
              std::vector<double> &gl_x,
              std::vector<double> &ov,
              std::vector<double> &spl,
              std::vector<double> &kkn,
              ModelV &Vptr);

int SplineInt(int n, int k,
              std::vector<double> &gl_w,
              std::vector<double> &gl_x,
              std::vector<double> &ov,
              int di, int dj,
              std::vector<double> &spl,
              std::vector<double> &splp,
              std::vector<double> &kkn,
              ModelV &V)            
}

#endif // BSPLINES_H_
#ifndef BSPLINES_H_
#define BSPLINES_H_

#include <cmath>
#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include "slatec_f.h"
#include "ModelV.h"
#include "H5Cpp.h"

namespace bsp {

int GenKnots(int n, int k, double r_max, 
             double fkn, char type,
             std::vector<double> &kkn);

int GenKnots(int n, int k, double r_max,
             std::string file, char type, 
             std::vector<double> &kkn);

int WrKnotsH5(int n, int k, double r_max, 
              double fkn, char type,
              std::string file,
              std::vector<double> &kkn);

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
              std::unique_ptr<ModelV> &Vptr);

int SplineInt(int n, int k,
              std::vector<double> &gl_w,
              std::vector<double> &gl_x,
              std::vector<double> &ov,
              int di, int dj,
              std::vector<double> &spl,
              std::vector<double> &splp,
              std::vector<double> &kkn,
              ModelV *V);          
}

#endif // BSPLINES_H_
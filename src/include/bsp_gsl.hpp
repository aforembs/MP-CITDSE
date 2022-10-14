#ifndef BSP_GSL_H_
#define BSP_GSL_H_

#include "H5Cpp.h"
#include "ModelV.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <cstring>
#include <vector>
extern "C" {
#include <gsl/gsl_bspline.h>
}

namespace bsp {

int GenKnots(int n, int k, double r_max, double fkn, char type,
             std::vector<double> &kkn);

int GenKnots(int n, int k, double r_max, std::string file, char type,
             std::vector<double> &kkn);

int WrKnotsH5(int n, int k, double r_max, double fkn, char type,
              std::string file, std::vector<double> &kkn);

int Splines(int n, int k, int glq_pt, std::vector<double> &gl_x,
            std::vector<double> &knots, std::vector<double> &splines);

int SimpSplines(int n, int k, int glq_pt, std::vector<double> &gl_x,
                std::vector<double> &knots, std::vector<double> &sisplines);

int Lob3Splines(int n, int k, int glq_pt, std::vector<double> &gl_x,
                std::vector<double> &knots, std::vector<double> &lsplines);

int Lob4Splines(int n, int k, int glq_pt, std::vector<double> &gl_x,
                std::vector<double> &knots, std::vector<double> &lobsplines);

int GL2Splines(int n, int k, int glq_pt, int glq_pt2,
               std::vector<double> &gl_outer, std::vector<double> &gl_inner,
               std::vector<double> &knots, std::vector<double> &glsplines);

int TrapSplines(int n, int k, int pt, std::vector<double> &knots,
                std::vector<double> &splines);

int SplinesP(int n, int k, int glq_pt, std::vector<double> &gl_x,
             std::vector<double> &knots, std::vector<double> &splinesp);

int SplineInt(int n, int k, int glq_pt, std::vector<double> &gl_w,
              std::vector<double> &gl_x, std::vector<double> &ov,
              std::vector<double> &spl, std::vector<double> &kkn,
              std::unique_ptr<ModelV> &Vptr);

int SplineInt(int n, int k, std::vector<double> &gl_w,
              std::vector<double> &gl_x, std::vector<double> &ov, int di,
              int dj, std::vector<double> &spl, std::vector<double> &splp,
              std::vector<double> &kkn, ModelV *V);
} // namespace bsp

#endif // BSP_GSL_H_
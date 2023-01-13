#ifndef BSP_GSL_H_
#define BSP_GSL_H_

/**
 * @file bsp_gsl.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for calculating the B-spline basis.
 * @version 1.0
 */

#include "H5Cpp.h"
#include "ModelV.hpp"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <vector>
extern "C" {
#include <gsl/gsl_bspline.h>
}

/**
 * @brief Namespace for B-spline basis funcitons.
 */
namespace bsp {
/**
 * @brief Function for generating a knot sequence over [0,R] using a given knot
 * distribution 'linear/sine/exponential'.
 *
 * @param n the number of knot intervals
 * @param k the maximum B-spline order
 * @param R_max the bounding radius of the simulated 'box'
 * @param fkn fequency parameter for sine/exponential distribution
 * @param type the type of knot distribution 'l/s/e'
 * @param kkn vector containing the knots
 * @return int default '0' error otherwise
 */
int genKnots(int n, int k, int R_max, double fkn, char type,
             std::vector<double> &kkn);

/**
 * @brief Function for reading a custom knot sequence over [0,R] from a user
 * supplied file.
 *
 * @param n the number of knot intervals
 * @param k the maximum B-spline order
 * @param R_max the bounding radius of the simulated 'box'
 * @param file the path of the file containing the custom knot distribution
 * @param type the type of the supplied file text/binary
 * @param kkn vector containing the knots
 * @return int default '0' error otherwise
 */
int genKnots(int n, int k, int R_max, std::string file, char type,
             std::vector<double> &kkn);

/**
 * @brief Function for saving the knot distribution to a HDF5 file.
 *
 * @param n the number of knot intervals
 * @param k the maximum B-spline order
 * @param R_max the bounding radius of the simulated 'box'
 * @param fkn fequency parameter for sine/exponential distribution
 * @param type the type of knot distribution 'l/s/e'
 * @param file path of the output file
 * @param kkn vector containing the knots
 * @return int default '0' error otherwise
 */
int wrKnotsH5(int n, int k, int R_max, double fkn, char type, std::string file,
              std::vector<double> &kkn);

/**
 * @brief Function for generating B-splines bases and their first order
 * derivatives within the predefined knot intervals on Gaussian quadrature
 * points.
 *
 * @param n the number of knot intervals
 * @param k the maximum B-spline order
 * @param glq_pt the order of the Gauss-Lagendre quadrature within each knot
 * interval
 * @param gl_x the [-1,1] positions of the quadrature points
 * @param kkn vector containing the knots
 * @param splines vector containing the B-splines bases at the 'gl_x' points
 * @param splinesp vector containing the first order derivatives of the
 * B-splines bases at the 'gl_x' points
 * @return int default '0' error otherwise
 */
int splines(int n, int k, int glq_pt, std::vector<double> &gl_x,
            std::vector<double> &kkn, std::vector<double> &splines,
            std::vector<double> &splinesp);

/**
 * @brief Function for calculating the integral over B-splines
 * '<B_i|Vptr->V(x)|B_j>' used in the solution of the 1-electron TISE
 *
 * @param n the number of knot intervals
 * @param k the maximum B-spline order
 * @param glq_pt the order of the Gauss-Lagendre quadrature within each knot
 * interval
 * @param gl_w the wights of the quadrature over the interval [-1,1]
 * @param gl_x the [-1,1] positions of the quadrature points
 * @param ov vector containing the resulting integral over B-splines
 * @param spl vector containing the B-splines bases at the 'gl_x' points
 * @param kkn vector containing the knots
 * @param Vptr pointer to a class containing the definition of the potential
 * function to be integrated V(x)
 * @return int default '0' error otherwise
 */
int splineInt(int n, int k, int glq_pt, std::vector<double> &gl_w,
              std::vector<double> &gl_x, std::vector<double> &ov,
              std::vector<double> &spl, std::vector<double> &kkn,
              std::unique_ptr<ModelV> &Vptr);
} // namespace bsp

#endif // BSP_GSL_H_
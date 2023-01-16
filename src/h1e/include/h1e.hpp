#ifndef H1E_H_
#define H1E_H_

/**
 * @file h1e.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for solving the 1-electron tdse on a B-splines
 * basis
 * @version 1.0
 */

#include "ModelV.hpp"
#include "bsp_gsl.hpp"
#include <H5Cpp.h>
#include <algorithm>
#include <execution>
#include <iostream>
#include <lapacke.h>
#include <omp.h>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

/**
 * @brief Namespace for functions used in soving the 1-electron TISE
 */
namespace h1e {
/**
 * @brief Function for reading settings from a .yaml file
 *
 * @param file path to the .yaml file
 * @param n the number of knot intervals = no. of 1-electron states + 2
 * @param k the maximum B-spline order
 * @param glq_pt the order of the Gauss-Lagendre quadrature within each knot
 * interval
 * @param R_max the bounding radius of the simulated 'box'
 * @param grid the type of knot distribution 'linear/sine/exponential/custom'
 * @param k_file the path of the file containing the custom knot distribution
 * @param pot the name of the selected potential
 * @param l_max the maximum 1-electron angular momentum used
 * @param z the value of the nuclear charge
 * @param mass 0.5 for atom 1 for positronium
 * @return int default '0' error otherwise
 */
int readConfig(std::string file, int &n, int &k, int &glq_pt, int &R_max,
               std::string &grid, std::string &k_file, std::string &pot,
               int &l_max, int &z, double &mass);

/**
 * @brief Function for solving the 1-electron TISE on a B-splines basis
 *
 * @param n the number of knot intervals = no. of 1-electron states + 2
 * @param k the maximum B-spline order
 * @param glq_pt the order of the Gauss-Lagendre quadrature within each knot
 * interval
 * @param l_max the maximum 1-electron angular momentum used
 * @param z the value of the nuclear charge
 * @param mass 0.5 for atom 1 for positronium
 * @param pot the name of the selected potential
 * @param gl_w the wights of the quadrature over the interval [-1,1]
 * @param gl_x the [-1,1] positions of the quadrature points
 * @param kkn vector containing the knots
 * @param spl vector containing the B-splines bases at the 'gl_x' points
 * @param splp vector containing the first order derivatives of the
 * @param outFile path of the HDF5 output file
 * @return int default '0' error otherwise
 */
int genCoeff(int n, int k, int glq_pt, int l_max, double z, double mass,
             std::string pot, std::vector<double> &gl_w,
             std::vector<double> &gl_x, std::vector<double> &kkn,
             std::vector<double> &spl, std::vector<double> &splp,
             std::string outFile);

} // namespace h1e

#endif // H1E_H_
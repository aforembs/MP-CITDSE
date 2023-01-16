#ifndef W1E_H
#define W1E_H

/**
 * @file w1e.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for obtaining the 1-electron wave functions on points
 * defined by the inner and outer quadratures for integration.
 * @version 1.0
 */

#include "fastgl.hpp"
#include <H5Cpp.h>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <yaml-cpp/yaml.h>
extern "C" {
#include <gsl/gsl_bspline.h>
}

/**
 * @brief
 */
namespace w1e {
/**
 * @brief Function for reading settings from a .yaml file
 *
 * @param file path to the .yaml file
 * @param qsz the size of the outer quadrature
 * @param R_max the bounding radius of the simulated 'box'
 * @param l_max the maximum 1-electron angular momentum used
 * @param pot the name of the selected potential
 * @param quad_type
 * @param quad_file
 * @param in_quad_layout
 * @param pt_file
 * @return int default '0' error otherwise
 */
int readConfig(std::string file, int &qsz, int &R_max, int &l_max,
               std::string &pot, std::string &quad_type, std::string &quad_file,
               std::string &in_quad_layout, std::string &pt_file);

/**
 * @brief Generate Gauss-Legendre nodes and weights for the outer quadrature.
 *
 * @param qsz the size of the outer quadrature
 * @param R_max the bounding radius of the simulated 'box'
 * @param q_x the points of the outer quadrature scaled for
 * (0,R_max)
 * @param q_w the weights of the outer quadrature on the 'qx' points scaled for
 * (0,R_max)
 * @return int default '0' error otherwise
 */
int genGaussLegendre(int qsz, int R_max, std::vector<double> &q_x,
                     std::vector<double> &q_w);

/**
 * @brief Generate the default inner quadrature point distribution.
 *
 * @param qsz the size of the outer quadrature
 * @param pq_dx vector containing the distribution of points for the inner
 * quadrature
 * @param pti_sz the number of inner quadrature points
 * @return int default '0' error otherwise
 */
int defaultPointLayout(int qsz, std::vector<uint8_t> &pq_dx, int &pti_sz);

/**
 * @brief Read a user supplied inner quadrature point distreibution.
 *
 * @param ftype point distribution file type 't'/'b' text or binary
 * @param pt_file path to the file containing the inner point distribution
 * @param qsz the size of the outer quadrature
 * @param pq_dx vector containing the distribution of points for the inner
 * quadrature
 * @param pti_sz the number of inner quadrature points
 * @return int default '0' error otherwise
 */
int userPointLayout(char ftype, std::string pt_file, int qsz,
                    std::vector<uint8_t> &pq_dx, int &pti_sz);

/**
 * @brief Read a user supplied outer quadrature and wieghts defined on (-1,1)
 * and transform to (0,R_max).
 *
 * @param qsz the size of the outer quadrature
 * @param R_max the bounding radius of the simulated 'box'
 * @param quad_file
 * @param type
 * @param q_x the points of the outer quadrature scaled for
 * (0,R_max)
 * @param q_w the weights of the outer quadrature on the 'qx' points scaled for
 * (0,R_max)
 * @return int default '0' error otherwise
 */
int readQuad(int qsz, int R_max, std::string quad_file, char type,
             std::vector<double> &q_x, std::vector<double> &q_w);

/**
 * @brief Generate and save the quadrature points, weights and the corresponding
 * values of the 1-electron wave functions
 *
 * @param pot the name of the selected potential
 * @param qsz the size of the outer quadrature
 * @param pti_sz the number of inner quadrature points
 * @param l_max the maximum 1-electron angular momentum used
 * @param q_x the points of the outer quadrature scaled for
 * (0,R_max)
 * @param q_w the weights of the outer quadrature on the 'qx' points scaled for
 * (0,R_max)
 * @param pq_dx vector containing the distribution of points for the inner
 * quadrature
 * @return int default '0' error otherwise
 */
int genWfn(std::string pot, int qsz, int pti_sz, int l_max,
           std::vector<double> &q_x, std::vector<double> &q_w,
           std::vector<uint8_t> &pq_dx);

} // namespace w1e

#endif // W1E_H
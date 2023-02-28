#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

/**
 * @file integrator.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for calculating the 't_{ab}^{v/l}' overlap integrals
 * required by the 1-electron dipole matrix elements
 * @version 1.0
 */

#include <vector>

/**
 * @brief Namespace fot the 1-electron dipole matrix integration functions
 */
namespace dmx_int {
/**
 * @brief Funciton for calulating the 't_{ab}' overlap integral between 2
 * 1-electron states in the velocity gauge
 *
 * @param qsz the size of the quadrature used for integration
 * @param lc_sz the number of wave function points for each 1-electron angular
 * momentum, used as an offset.
 * @param n1 the primary quantum number of the initial electronic state
 * @param l1 the angular momentum quantum number of the initial electronic state
 * @param n2 the primary quantum number of the final electronic state
 * @param l2 the angular momentum quantum number of the final electronic state
 * @param qx the positions of the quadrature points on the interval [0,R_max]
 * @param qw the weights of the quadrature on the 'qx' points scaled for
 * (0,R_max)
 * @param wfn the values of the 1-electron wave functions on the quadrature
 * points 'qx'
 * @param wfnp the values of the first derivatives of the 1-electron wave
 * functions on the quadrature points 'qx'
 * @return double the value of the overlap integral 't_{ab}^v'
 */
double tvelGL(int qsz, int lc_sz, int n1, int l1, int n2, int l2,
              std::vector<double> &qx, std::vector<double> &qw,
              std::vector<double> &wfn, std::vector<double> &wfnp);

/**
 * @brief Funciton for calulating the 't_{ab}' overlap integral between 2
 * 1-electron states in the length gauge
 *
 * @param qsz the size of the quadrature used for integration
 * @param lc_sz the number of wave function points for each 1-electron angular
 * momentum, used as an offset.
 * @param n1 the primary quantum number of the initial electronic state
 * @param l1 the angular momentum quantum number of the initial electronic state
 * @param n2 the primary quantum number of the final electronic state
 * @param l2 the angular momentum quantum number of the final electronic state
 * @param qx the positions of the quadrature points on the interval [0,R_max]
 * @param qw the weights of the quadrature on the 'qx' points scaled for
 * (0,R_max)
 * @param wfn the values of the 1-electron wave functions on the quadrature
 * points 'qx'
 * @return double the value of the overlap integral 't_{ab}^l'
 */
double tlenGL(int qsz, int lc_sz, int n1, int l1, int n2, int l2,
              std::vector<double> &qx, std::vector<double> &qw,
              std::vector<double> &wfn);
} // namespace dmx_int

#endif // INTEGRATOR_H_
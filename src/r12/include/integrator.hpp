#ifndef INTEGRATOR_NEW_H_
#define INTEGRATOR_NEW_H_

/**
 * @file integrator.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for calculating the photoelectron energy distribution
 * from 1-electron TDSE coefficients.
 * @version 1.0
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

/**
 * @brief Namespace for functions used the calculate the radial slater integrals
 * used in the computiation of 2-electron correlations.
 */
namespace intfn {
/**
 * @brief
 *
 * @param k the order of the multipole series
 * @param qsz the size of the quadrature used for integration
 * @param pti_sz the number of inner quadrature points
 * @param lc_sz the number of wave function points in the outer integral for
 * each 1-electron angular momentum, used as an offset.
 * @param lci_sz the number of wave function points in the inner integral for
 * each 1-electron angular momentum, used as an offset.
 * @param n1 the primary quantum number of the 1st electron in the 1st
 * configuration
 * @param l1 the angular momentum quantum number of the 1st electron in the 1st
 * configuration
 * @param n2 the primary quantum number of the 2nd electron in the 1st
 * configuration
 * @param l2 the angular momentum quantum number of the 2nd electron in the 1st
 * configuration
 * @param n3 the primary quantum number of the 1st electron in the 2nd
 * configuration
 * @param l3 the angular momentum quantum number of the 1st electron in the 2nd
 * configuration
 * @param n4 the primary quantum number of the 2nd electron in the 2nd
 * configuration
 * @param l4 the angular momentum quantum number of the 2nd electron in the 2nd
 * configuration
 * @param q_w the weights of the outer quadrature on the 'qx' points scaled for
 * (0,R_max)
 * @param pq_dx vector containing the distribution of points for the inner
 * quadrature
 * @param rk vector containing the outer quadrature point positions r^k where
 * 0<=k<=k_max
 * @param rk_in vector containing the inner quadrature point positions r^k where
 * 0<=k<=k_max
 * @param wfn_o the values of the 1-electron wave functions on the outer
 * quadrature points
 * @param wfn_i the values of the 1-electron wave functions on the inner
 * quadrature points
 * @return double the value of the slater integral for the given set of
 * configurations
 */
double fsltrLob(int k, int qsz, int pti_sz, int lc_sz, int lci_sz, int n1,
                int l1, int n2, int l2, int n3, int l3, int n4, int l4,
                std::vector<double> &q_w, std::vector<uint8_t> &pq_dx,
                std::vector<double> &rk, std::vector<double> &rk_in,
                std::vector<double> &wfn_o, std::vector<double> &wfn_i);

} // namespace intfn

#endif // INTEGRATOR_H_
#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>

// needs an update
double FsltrGL2(int k, int n, int bo, int gl2,
            std::vector<double> &gl_ow, 
            std::vector<double> &gl_ox,
            std::vector<double> &gl_iw, 
            std::vector<double> &gl_ix, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Gsp,
            std::vector<double> &Cl1i_pt,
            std::vector<double> &Cl1p_pt,
            std::vector<double> &Cl2i_pt,
            std::vector<double> &Cl2p_pt);

// needs an update 
double FsltrLob4GL(int k, int n, int bo,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Ssp,
            std::vector<double> &Cl1i_pt,
            std::vector<double> &Cl1p_pt,
            std::vector<double> &Cl2i_pt,
            std::vector<double> &Cl2p_pt);

double FsltrLob3GL(int k, int n, int bo,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Ssp,
            std::vector<double> &rk,
            std::vector<double> &rk_mid,
            std::vector<double> &Cl1i_pt,
            std::vector<double> &Cl2i_pt,
            std::vector<double> &p1p_buff,
            std::vector<double> &p2p_buff,
            std::vector<double> &p2p_mid,
            std::vector<double> &p2is);

#endif // INTEGRATOR_H_
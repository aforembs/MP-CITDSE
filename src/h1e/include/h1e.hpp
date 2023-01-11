#ifndef H1E_H_
#define H1E_H_

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

namespace h1e {

int ReadConfig(std::string file, int &n, int &k, int &glq_pt, int &r_max,
               std::string &grid, std::string &k_file, std::string &pot,
               int &l_max, int &z, double &mass);

int GenCoeff(int n, int k, int glq_pt, int l, double z, double mass,
             std::string pot, std::vector<double> &gl_w,
             std::vector<double> &gl_x, std::vector<double> &kkn,
             std::vector<double> &spl, std::vector<double> &splp,
             std::string outFile);

} // namespace h1e

#endif // H1E_H_
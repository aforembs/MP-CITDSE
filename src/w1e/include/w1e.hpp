#ifndef W1E_H
#define W1E_H

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

namespace w1e {

int readConfig(std::string file, int &qsz, int &R_max, int &l_max,
               std::string &pot, std::string &quad_type, std::string &quad_file,
               std::string &in_quad_layout, std::string &pt_file);

int genGaussLegendre(int qsz, int R_max, std::vector<double> &q_x,
                     std::vector<double> &q_w);

int defaultPointLayout(int qsz, std::vector<uint8_t> &pq_dx, int &pti_sz);

int userPointLayout(char ftype, std::string pt_file, int qsz,
                    std::vector<uint8_t> &pq_dx, int &pti_sz);

int readQuad(int qsz, int R_max, std::string quad_file, char type,
             std::vector<double> &q_x, std::vector<double> &q_w);

int genWfn(std::string pot, int qsz, int pti_sz, int l_max,
           std::vector<double> &q_x, std::vector<double> &q_w,
           std::vector<uint8_t> &pq_dx);

} // namespace w1e

#endif // W1E_H
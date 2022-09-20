#include "H5Cpp.h"
#include "bsplines.hpp"
#include "fastgl.hpp"
#include "slatec_f.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

int bsptest(std::string cpot, uint l1e_max, std::vector<uint> &N_max) {
  int n = 0;   // no. of points
  int k = 0;   // max B-spline order
  int nkn = 0; // no. of knots
  int nSt = 0; // no. of states
  std::vector<double> kkn;
  std::vector<double> Cf;
  std::vector<double *> C(l1e_max + 1);

  std::string filename;
  H5::DataSpace cspace;
  hsize_t offset[2], count[2], stride[2], block[2];
  hsize_t dimms[2];
  offset[0] = 0;
  offset[1] = 0;
  stride[0] = 1;
  stride[1] = 1;
  block[0] = 1;
  block[1] = 1;
  H5::DataSpace memspace;

  int tot_states = 0;
  for (auto &b : N_max) {
    tot_states += b;
  }

  filename = cpot + std::to_string(0) + ".h5";
  auto file =
      std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  file->openAttribute("N").read(H5::PredType::NATIVE_INT32, &n);
  file->openAttribute("K").read(H5::PredType::NATIVE_INT32, &k);
  auto rset =
      std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Knots")));
  nkn = rset->getSpace().getSimpleExtentNpoints();

  kkn.reserve(nkn);
  rset->read(&kkn[0], H5::PredType::NATIVE_DOUBLE); // read knots

  rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("En")));
  nSt = rset->getSpace().getSimpleExtentNpoints();

  Cf.reserve(tot_states * n);
  C[0] = &Cf[0];
  int nst_prev = 0;
  for (uint i = 1; i <= l1e_max; ++i) {
    nst_prev += N_max[i - 1];
    C[i] = &Cf[nst_prev * n];
  }
  std::cout << n << " " << nSt << " " << nkn << "\n";
  // read coefficients for all l
  for (uint l = 0; l <= l1e_max; ++l) {
    count[0] = N_max[l];
    count[1] = n;
    dimms[0] = count[0];
    dimms[1] = count[1];
    memspace.setExtentSimple(2, dimms, NULL);

    filename = cpot + std::to_string(l) + ".h5";
    file =
        std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    rset = std::unique_ptr<H5::DataSet>(
        new H5::DataSet(file->openDataSet("Coeff")));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(C[l], H5::PredType::NATIVE_DOUBLE, memspace, cspace);
  }

  std::vector<double> gl_x(k);
  for (int i = 1; i <= k; ++i) {
    gl_x[k - i] =
        fastgl::GLPair(k, i).x(); // generate GL nodes over B-splines support
  }

  std::vector<double> Bsplines;
  bsp::Splines(n, k, gl_x, kkn, Bsplines);
  double dl, sl, x;

  std::cout << Cf[0] << " " << Cf[n] << " " << Cf[2 * n] << "\n";

  std::ofstream outFile("dat/wfn_n10l0.dat", std::ofstream::out);
  // Output ground state wfn
  double x_val = 0;
  for (auto i = k - 1; i < n; ++i) {
    dl = (kkn[i + 1] - kkn[i]) * 0.5;
    sl = (kkn[i + 1] + kkn[i]) * 0.5;

    for (int p = 0; p < k; ++p) {
      x = dl * gl_x[p] + sl; // x-transformation
      x_val = 0;

      for (int j = 0; j < k; ++j) {
        x_val += Cf[n + i - k + 1 + j] * Bsplines[j + k * (p + i * k)];
      }
      outFile << std::setiosflags(std::ios::scientific) << std::setprecision(12)
              << x << " " << x_val << "\n";
    }
  }
  outFile.close();
  return 0;
}

int main() {
  std::vector<uint> Nm;
  Nm.push_back(1);
  Nm.push_back(1);
  bsptest("dat/he", 1, Nm);
  return 0;
}
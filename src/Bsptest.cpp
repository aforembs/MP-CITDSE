#include <boost/math/special_functions/cardinal_b_spline.hpp> //?
#include "fastgl.h"

int bsptest(std::string cpot, uint l1e_max, std::vector<uint> &N_max) {
  int n=0;   // no. of points
  int k=0;   // max B-spline order
  int nkn=0; // no. of knots
  int nCf=0; // no. of coefficients
  int nSt=0; // no. of states
  std::vector<double> kkn;
  std::vector<double> Cf;
  std::vector<double*> C(l1e_max+1);

  std::string filename;
  H5::H5File    *file=nullptr;
  H5::DataSet   *rset=nullptr;
  H5::DataSpace cspace;
  hsize_t offset[1], count[1], stride[1], block[1];
  hsize_t dimms[1];
  offset[0]=0;
  count[0] =N_sz[0];
  stride[0]=1;
  block[0] =1;
  dimms[0] =N_sz[0];
  H5::DataSpace memspace;

  int tot_states = 0;
  for(auto &n : N_max) {
    tot_states += n;
  }

  filename = cpot + std::to_string(0) + ".h5";
  file = new H5::H5File(filename, H5F_ACC_RDONLY);
  file->openAttribute("N").read(H5::PredType::NATIVE_INT32, &n);
  file->openAttribute("K").read(H5::PredType::NATIVE_INT32, &k);
  rset = new H5::DataSet(file->openDataSet("Knots"));
  nkn = rset->getSpace().getSimpleExtentNpoints();

  kkn.reserve(nknots);
  rset->read(&kkn[0], H5::PredType::NATIVE_DOUBLE);
  delete rset;

  rset = new H5::DataSet(file->openDataSet("En"));
  nSt = rset->getSpace().getSimpleExtentNpoints();
  delete rset;
  delete file;

  Cf.reserve(tot_states*nSt);
  C[0]=&Cf[0];
  int nst_prev = 0;
  for(int i=1; i<=l1e_max; ++i) {
    nst_prev += N_max[i-1];
    C[i] = &Cf[nst_prev*nSt];
  }

  for(int l=0; l<=l1e_max; ++l) {
    count[0] = N_max[l]*nSt;
    dimms[0] = count[l];
    memspace.setExtentSimple(1, dimms, NULL);

    filename = cpot + std::to_string(l) + ".h5";
    file = new H5::H5File(filename, H5F_ACC_RDONLY);
    rset = new H5::DataSet(file->openDataSet("Coeff"));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(C[l], H5::PredType::NATIVE_DOUBLE, memspace, cspace);
    delete rset;
    delete file;
  }

  std::vector<fastgl::QuadPair> pvec(k);
  for(int i=1; i<=k; ++i) {
    pvec[i-1] = fastgl::GLPair(k, i);
  }

}
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "H5Cpp.h"
#include "fastgl.h"
#include "slatec_f.h"
#include "bsplines.h"

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
  hsize_t offset[2], count[2], stride[2], block[2];
  hsize_t dimms[2];
  offset[0]=0; offset[1]=0;
  count[0] =0;
  stride[0]=1; stride[1]=1;
  block[0] =1; block[1]=1;
  dimms[0] =0;
  H5::DataSpace memspace;

  int tot_states = 0;
  for(auto &b : N_max) {
    tot_states += b;
  }

  filename = cpot + std::to_string(0) + ".h5";
  file = new H5::H5File(filename, H5F_ACC_RDONLY);
  file->openAttribute("N").read(H5::PredType::NATIVE_INT32, &n);
  file->openAttribute("K").read(H5::PredType::NATIVE_INT32, &k);
  rset = new H5::DataSet(file->openDataSet("Knots"));
  nkn = rset->getSpace().getSimpleExtentNpoints();

  kkn.reserve(nkn);
  rset->read(&kkn[0], H5::PredType::NATIVE_DOUBLE); //read knots
  delete rset;

  rset = new H5::DataSet(file->openDataSet("En"));
  nSt = rset->getSpace().getSimpleExtentNpoints();
  delete rset;
  delete file;

  Cf.reserve(tot_states*n);
  C[0]=&Cf[1];
  Cf[0]=0.0; Cf[nSt+1] =0.0;
  int nst_prev = 0;
  for(int i=1; i<=l1e_max; ++i) {
    nst_prev += N_max[i-1];
    C[i] = &Cf[nst_prev*n+1];
    Cf[nst_prev*n]=0.0;
    Cf[i*n-1]=0.0;
  }
  std::cout << n << " " << nSt << " " << nkn <<"\n";
  // read coefficients for all l
  for(int l=0; l<=l1e_max; ++l) {
    count[0] = N_max[l]; count[1] = nSt;
    dimms[0] = count[0]; dimms[1] = count[1];
    memspace.setExtentSimple(2, dimms, NULL);

    filename = cpot + std::to_string(l) + ".h5";
    file = new H5::H5File(filename, H5F_ACC_RDONLY);
    rset = new H5::DataSet(file->openDataSet("Coeff"));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(C[l], H5::PredType::NATIVE_DOUBLE, memspace, cspace);
    delete rset;
    delete file;
  }

  std::vector<double> gl_x(k);
  for(int i=1; i<=k; ++i) {
    gl_x[k-i] = fastgl::GLPair(k, i).x(); // generate GL nodes over B-splines support
  }

  int len = (k+1)*(k+2)/2;
  std::vector<double> Db(k);
  std::vector<double> work(len);
  std::vector<double> Bsplines;

  Bsplines.reserve(nSt*k*k);
  bsplines(n, k, gl_x, kkn, Bsplines);
  int i1 ;
  double dl, sl, x;
             
  // for(auto i=k-1; i<n; ++i){
  //   dl = (kkn[i+1] - kkn[i])*0.5;
  //   sl = (kkn[i+1] + kkn[i])*0.5;

  //   for(int p=0; p<k; ++p){
  //     x = dl*gl_x[p] + sl;    //x-transformation
  //     i1 = i + 1 ;
  //     dbspvd_(&kkn[0], k, 1, x, i1, k, &Db[0], &work[0]);

  //     Bsplines.insert(std::end(Bsplines), std::begin(Db), std::end(Db));
  //   }
  // }

  std::ofstream outFile("wfn_l0.dat", std::ofstream::out);
  // Output ground state wfn
  int bidx=0;
  double x_val=0, x_part=0;
  for(auto i=k-1; i<n; ++i, ++bidx){
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;

    for(int p=0; p<k; ++p){
      x = dl*gl_x[p] + sl;    //x-transformation
      x_val = 0;

      for(int j=0; j<k; ++j) {
        x_part = Cf[i-k+1+j]*Bsplines[j+k*(p+bidx*k)];
        x_val += x_part;//Cf[i-k+1+j]*Bsplines[j+k*(p+bidx*k)];
        // if (p==4 &&j==3) {
        // outFile << std::setiosflags(std::ios::scientific)
        //         << std::setprecision(12) << x << " " << x_part << "\n";}
      }
      outFile << std::setiosflags(std::ios::scientific)
              << std::setprecision(12) << x << " " << x_val << "\n";
    }
  }
  outFile.close();
  return 0;
}

int main() {
  std::vector<uint> Nm;
  Nm.push_back(1);
  Nm.push_back(50);
  bsptest("he", 1, Nm);
  return 0;
}
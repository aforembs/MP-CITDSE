#include "tise.h"

int tise::WriteHdf5() {

  outFile = outFile + _potential + Global::app(l_) + _tail ;

   // Create and leave in define mode
  H5::H5File *file = new H5::H5File(outFile, H5F_ACC_TRUNC);

  // Check if the file was opened
  if (!file) {
    cerr << "# H5::H5File:: file couldn't opened: " << outFile.c_str() << "\n";
    exit(-1);
  }

  // Create dimensions
  hsize_t n_d[1]       = {n_};
  hsize_t m_d[1]       = {n_};
  hsize_t nKnots_d[1]  = {grid.Knots().size()};
  hsize_t att_space[1] = {1};
  hsize_t sqr_space[2] = {n_, n_};

  H5::Attribute Z = file->createAttribute("Z", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, att_space));
  H5::Attribute M = file->createAttribute("M", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, att_space));
  H5::Attribute N = file->createAttribute("N", H5::PredType::NATIVE_INT32, H5::DataSpace(1, att_space));
  H5::Attribute K = file->createAttribute("K", H5::PredType::NATIVE_INT32, H5::DataSpace(1, att_space));
  H5::Attribute R = file->createAttribute("R", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, att_space));
  H5::Attribute l = file->createAttribute("l", H5::PredType::NATIVE_INT32, H5::DataSpace(1, att_space));

  // Create variables
  H5::DataSet Knots = file->createDataSet("Knots", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, nKnots_d));
  H5::DataSet E_nl = file->createDataSet("En", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, n_d));
  H5::DataSet C_nl = file->createDataSet("Coeff", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, sqr_space));

  // Now write in the netCDF file
  double z, mass ;
  int k,n ;
  double r ;
  z = z_ ;
  mass = mass_ ;
  k = grid.k() ;
  n = grid.Knots().size() - grid.k() ;                    // Number of B-splines
  r = grid.Knots()[grid.Knots().size() - 1] ;

  double * kkn = new double[nKnots_d[0]];
  for(int i = 0 ; i < nKnots_d[0] ;i++)       kkn[i] = grid.Knots()[i];

  Z.write(H5::PredType::NATIVE_DOUBLE, &z);
  M.write(H5::PredType::NATIVE_DOUBLE, &mass);
  N.write(H5::PredType::NATIVE_INT32, &n);
  K.write(H5::PredType::NATIVE_INT32, &k);
  R.write(H5::PredType::NATIVE_DOUBLE, &r);
  Knots.write(kkn, H5::PredType::NATIVE_DOUBLE);
  l.write(H5::PredType::NATIVE_INT32, &l_);
  E_nl.write(Enl, H5::PredType::NATIVE_DOUBLE);
  C_nl.write(Cnl, H5::PredType::NATIVE_DOUBLE);

  std::cout << "# write:: HDF5 DATA FOR L =  "<< l_ << " STORED IN "<< outFile << '\n' << "\n";
  delete kkn;
  delete file;
  return 0;
}

int tise::GenCoeff(
                  std::vector<double> &gl_w, 
                  std::vector<double> &gl_x,
                  ) {
  int nm2=n-2;
  int nk=nm2*k;
  int llp1=0, nik=0, ni2=0;
  std::vector<double> ov_BB(nk), ov_dBdB(nk), ov_1_r2(nk), ov_V(nk), aa;
  std::vector<double> Enl, Cnl(n*nm2), Cnl_tmp;

  ModelV *v_1    = new V_c(1.0);
  ModelV *v_1_r2 = new V_c_r2(1.0);

  bsp::SplineInt(nm2, k, gl_w, gl_x, ov_BB, spl, v_1); // int B_iB_j dr
  bsp::SplineInt(nm2, k, gl_w, gl_x, ov_dBdB, splp, v_1); // int B_i d/dr^2 B_j dr
  bsp::SplineInt(nm2, k, gl_w, gl_x, ov_1_r2, spl, v_1_r2); // int B_iB_j/r^2 dr

  ModelV *v  = new  V_1_r(z);

  bsp::SplineInt(n, k, gl_w, gl_x, ov_V, spl, v); // int B_i V(r) B_j dr

  aa.reserve(nk);
  Enl.reserve(nm2);
  Cnl_tmp.reserve(nm2*nm2);

  for(int l=0, l<=l_max; ++l) {
    llp1=l*(l+1);
    for(int ni=0; ni<nm2; ++ni) {
      nik = ni*k;
      for(int j=0; j<k; ++j) {
        aa[j+nik] = mass*ov_dBdB[j+nik] - ov_V[j+nik] + mass*llp1*ov_1_r2[j+nik];
      }
    }

    int info = LAPACKE_dsbgv(LAPACK_COL_MAJOR, 'V', 'U', nm2, k-1, k-1, &aa[0],
                            k, &ov_BB[0], k, &Enl[0], &Cnl_tmp[0], nm2);
    // Reshape with zeros at r=0 & r=R
    for(int ni=0; ni<nm2; ++ni) {
      ni2=ni*nm2;
      std::copy(std::execution::par_unseq, Cnl_tmp.begin()+ni2, 
                Cnl_tmp.begin()+ni2+nm2, Cnl.begin()+1+ni*n);
    }
    // Write hdf5 file

  }

}
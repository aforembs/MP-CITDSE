#include "tise.h"

// for now same as Lampros', remove knots soon!!!
int WriteHdf5(int n, int k, int li,
              double z, double mass,
              std::string pot, 
              std::vector<double> &kkn,
              std::vector<double> &Enl,
              std::vector<double> &Cnl,
              std::string outFile) {
  auto nm2=n-2;
  outFile = outFile + pot + std::to_string(li) + ".h5" ;

   // Create and leave in define mode
  std::unique_ptr<H5::H5File> file(new H5::H5File(outFile, H5F_ACC_TRUNC));

  // Check if the file was opened
  if (!file) {
    std::cerr << "# H5::H5File:: file couldn't opened: " << outFile.c_str() << "\n";
    exit(-1);
  }

  // Create dimensions
  hsize_t n_d[1]       = {(hsize_t)nm2};
  hsize_t nKnots_d[1]  = {kkn.size()};
  hsize_t att_space[1] = {1};
  hsize_t sqr_space[2] = {(hsize_t)nm2, (hsize_t)n};

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
  double r = kkn[kkn.size() - 1];

  Z.write(H5::PredType::NATIVE_DOUBLE, &z);
  M.write(H5::PredType::NATIVE_DOUBLE, &mass);
  N.write(H5::PredType::NATIVE_INT32, &n);
  K.write(H5::PredType::NATIVE_INT32, &k);
  R.write(H5::PredType::NATIVE_DOUBLE, &r);
  Knots.write(&kkn[0], H5::PredType::NATIVE_DOUBLE);
  l.write(H5::PredType::NATIVE_INT32, &li);
  E_nl.write(&Enl[0], H5::PredType::NATIVE_DOUBLE);
  C_nl.write(&Cnl[0], H5::PredType::NATIVE_DOUBLE);

  std::cout << "# write:: HDF5 DATA FOR L =  "<< li << " STORED IN "<< outFile << "\n\n";
  return 0;
}

int tise::ReadConfig(std::string file,
                    int &n, int &k, int &r_max,
                    std::string &grid,
                    std::string &k_file,
                    std::string &pot,
                    int &l_max, int &z, 
                    double &mass) {

  YAML::Node settings = YAML::LoadFile(file);

  n = settings["Basis_Settings"]["state_no"].as<int>();
  std::cout << "Number of States:     "
            << n << std::endl;
  k = settings["Basis_Settings"]["max_spline_k"].as<int>();
  std::cout << "Maximum B-splines order:  "
            << k << std::endl;
  r_max = settings["Basis_Settings"]["R_max"].as<int>();
  std::cout << "Box radius:        "
            << r_max << std::endl;
  grid = settings["Basis_Settings"]["grid"].as<std::string>();
  std::cout << "Type of knot spacing:            "
            << grid << std::endl;
  if(grid.compare("custom")==0) {
    k_file = settings["Basis_Settings"]["grid"].as<std::string>();
    std::cout << "Custom knot sequence file:     "
              << k_file << std::endl;
  }
  pot = settings["Basis_Settings"]["potential"].as<std::string>();
  std::cout << "Core Potential:            "
            << pot << std::endl;
  l_max     = settings["Basis_Settings"]["l_max"].as<int>();
  std::cout << "Maximum l:                         "
            << l_max << std::endl;
  z = settings["Basis_Settings"]["atomic_no"].as<int>();
  std::cout << "Atomic number:  "
            << z << std::endl;
  mass = settings["Basis_Settings"]["mass"].as<double>();
  std::cout << "Atomic number:  "
            << mass << std::endl;
  return 0;
}

int tise::GenCoeff(int n, int k, int l_max,
                  double z, double mass,
                  std::string pot,
                  std::vector<double> &gl_w, 
                  std::vector<double> &gl_x,
                  std::vector<double> &kkn,
                  std::vector<double> &spl,
                  std::vector<double> &splp,
                  std::string outFile) {
  int nm2=n-2;
  int nk=nm2*k;
  int llp1=0, nik=0, ni2=0;
  std::vector<double> ov_BB(nk), ov_dBdB(nk), ov_1_r2(nk), ov_V(nk);//, aa;
  //std::vector<double> Enl, Cnl(n*nm2), Cnl_tmp;
  //std::vector<double> w_bb(nk);

  // Change these eventually
  ModelV *v_1    = new V_c(1.0);
  ModelV *v_1_r2 = new V_c_r2(1.0);

  bsp::SplineInt(nm2, k, gl_w, gl_x, ov_BB, spl, kkn, v_1); // int B_iB_j dr
  bsp::SplineInt(nm2, k, gl_w, gl_x, ov_dBdB, splp, kkn, v_1); // int B_i d/dr^2 B_j dr
  bsp::SplineInt(nm2, k, gl_w, gl_x, ov_1_r2, spl, kkn, v_1_r2); // int B_iB_j/r^2 dr

  ModelV *v  = new  V_1_r(z);

  bsp::SplineInt(nm2, k, gl_w, gl_x, ov_V, spl, kkn, v); // int B_i V(r) B_j dr

  delete v_1;
  delete v_1_r2;
  delete v;

  // aa.reserve(nk); // possibly *l and parallel
  // Enl.reserve(nm2);
  // Cnl_tmp.reserve(nm2*nm2);

  omp_set_num_threads(std::min(l_max+1, omp_get_max_threads()));
  #pragma omp parallel for private(llp1, nik, ni2)
  for(int l=0; l<=l_max; ++l) {
    std::vector<double> Enl, Cnl(n*nm2), Cnl_tmp, aa, w_bb;
    Enl.reserve(nm2);
    Cnl_tmp.reserve(nm2*nm2);
    aa.reserve(nk);
    w_bb.reserve(nk);
    std::copy( ov_BB.begin(), ov_BB.end(), w_bb.begin());
    llp1=l*(l+1);
    for(int ni=0; ni<nm2; ++ni) {
      nik = ni*k;
      for(int j=0; j<k; ++j) {
        aa[j+nik] = mass*ov_dBdB[j+nik] - ov_V[j+nik] + mass*llp1*ov_1_r2[j+nik];
      }
    }

    LAPACKE_dsbgvd(LAPACK_COL_MAJOR, 'V', 'U', nm2, k-1, k-1, &aa[0],
                  k, &w_bb[0], k, &Enl[0], &Cnl_tmp[0], nm2);

    // Reshape with zeros at r=0 & r=R
    for(int ni=0; ni<nm2; ++ni) {
      ni2=ni*nm2;
      std::copy( Cnl_tmp.begin()+ni2, 
                Cnl_tmp.begin()+ni2+nm2, Cnl.begin()+1+ni*n);
    }
    // Write hdf5 file
    WriteHdf5(n, k, l, z, mass, pot, kkn, Enl, Cnl, outFile);
  }

  return 0;
}
#include "h1e.hpp"

// for now same as Lampros'
int writeHdf5(int n, int k, int li, double z, double mass, std::string pot,
              std::vector<double> &kkn, std::vector<double> &Enl,
              std::vector<double> &Cnl, std::string outFile) {
  auto nm2 = n - 2;
  outFile = outFile + pot + std::to_string(li) + ".h5";

  // Create and leave in define mode
  auto file =
      std::unique_ptr<H5::H5File>(new H5::H5File(outFile, H5F_ACC_TRUNC));

  // Check if the file was opened
  if (!file) {
    std::cerr << "# H5::H5File:: file couldn't opened: " << outFile.c_str()
              << "\n";
    exit(-1);
  }

  // Create dimensions
  hsize_t n_d[1] = {(hsize_t)nm2};
  hsize_t nKnots_d[1] = {kkn.size()};
  hsize_t att_space[1] = {1};
  hsize_t sqr_space[2] = {(hsize_t)nm2, (hsize_t)n};

  H5::Attribute Z = file->createAttribute("Z", H5::PredType::NATIVE_DOUBLE,
                                          H5::DataSpace(1, att_space));
  H5::Attribute M = file->createAttribute("M", H5::PredType::NATIVE_DOUBLE,
                                          H5::DataSpace(1, att_space));
  H5::Attribute N = file->createAttribute("N", H5::PredType::NATIVE_INT32,
                                          H5::DataSpace(1, att_space));
  H5::Attribute K = file->createAttribute("K", H5::PredType::NATIVE_INT32,
                                          H5::DataSpace(1, att_space));
  H5::Attribute R = file->createAttribute("R", H5::PredType::NATIVE_DOUBLE,
                                          H5::DataSpace(1, att_space));
  H5::Attribute l = file->createAttribute("l", H5::PredType::NATIVE_INT32,
                                          H5::DataSpace(1, att_space));

  // Create variables
  H5::DataSet Knots = file->createDataSet("Knots", H5::PredType::NATIVE_DOUBLE,
                                          H5::DataSpace(1, nKnots_d));
  H5::DataSet E_nl = file->createDataSet("En", H5::PredType::NATIVE_DOUBLE,
                                         H5::DataSpace(1, n_d));
  H5::DataSet C_nl = file->createDataSet("Coeff", H5::PredType::NATIVE_DOUBLE,
                                         H5::DataSpace(2, sqr_space));

  // Now write in the netCDF file
  double r = kkn[kkn.size() - 1];

  Z.write(H5::PredType::NATIVE_DOUBLE, &z);
  M.write(H5::PredType::NATIVE_DOUBLE, &mass);
  N.write(H5::PredType::NATIVE_INT32, &n);
  K.write(H5::PredType::NATIVE_INT32, &k);
  R.write(H5::PredType::NATIVE_DOUBLE, &r);
  Knots.write(&kkn[0], H5::PredType::NATIVE_DOUBLE);
  l.write(H5::PredType::NATIVE_INT32, &li);
  E_nl.write(&Enl[li * nm2], H5::PredType::NATIVE_DOUBLE);
  C_nl.write(&Cnl[li * n * nm2], H5::PredType::NATIVE_DOUBLE);

  std::cout << "# write:: HDF5 DATA FOR L =  " << li << " STORED IN " << outFile
            << "\n\n";
  return 0;
}

int h1e::readConfig(std::string file, int &n, int &k, int &glq_pt, int &r_max,
                    std::string &grid, std::string &k_file, std::string &pot,
                    int &l_max, int &z, double &mass) {

  YAML::Node settings = YAML::LoadFile(file);

  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "Core Potential:                            " << pot
            << std::endl;

  n = settings["Basis_Settings"]["state_no"].as<int>();
  std::cout << "Number of States:                          " << n << std::endl;
  k = settings["Basis_Settings"]["max_spline_k"].as<int>();
  std::cout << "Maximum B-splines order:                   " << k << std::endl;
  glq_pt = settings["Basis_Settings"]["GL_quad_points"].as<int>();
  std::cout << "No. of quadrature points per knot:         " << glq_pt
            << std::endl;
  r_max = settings["Basis_Settings"]["R_max"].as<int>();
  std::cout << "Box radius:                                " << r_max
            << std::endl;
  grid = settings["Basis_Settings"]["grid"].as<std::string>();
  std::cout << "Type of knot spacing:                      " << grid
            << std::endl;
  if (grid.compare("user-defined") == 0) {
    k_file = settings["Basis_Settings"]["grid"].as<std::string>();
    std::cout << "Custom knot sequence file:     " << k_file << std::endl;
  }
  l_max = settings["Basis_Settings"]["l_max"].as<int>();
  std::cout << "Maximum l:                                 " << l_max
            << std::endl;
  z = settings["Basis_Settings"]["atomic_no"].as<int>();
  std::cout << "Atomic number:                             " << z << std::endl;
  mass = settings["Basis_Settings"]["mass"].as<double>();
  std::cout << "mass (0.5 atoms, 1 positronium):           " << mass
            << std::endl;
  return 0;
}

int h1e::genCoeff(int n, int k, int glq_pt, int l_max, double z, double mass,
                  std::string pot, std::vector<double> &gl_w,
                  std::vector<double> &gl_x, std::vector<double> &kkn,
                  std::vector<double> &spl, std::vector<double> &splp,
                  std::string outFile) {
  int nm2 = n - 2;
  int lm1 = l_max + 1;
  int nm22 = nm2 * nm2;
  int nk = nm2 * k;
  int llp1 = 0, nik = 0, lnk = 0;
  std::vector<double> ov_BB(nk), ov_dBdB(nk), ov_1_r2(nk), ov_V(nk);
  std::vector<double> Enl, Cnl(lm1 * n * nm2), Cnl_tmp, aa, w_bb;

  // Change these eventually
  auto v_1 = std::unique_ptr<ModelV>(new V_c(1.0));
  auto v_1_r2 = std::unique_ptr<ModelV>(new V_c_r2(1.0));

  // int B_iB_j dr
  bsp::splineInt(nm2, k, glq_pt, gl_w, gl_x, ov_BB, spl, kkn, v_1);
  // int B_i d/dr^2 B_j dr
  bsp::splineInt(nm2, k, glq_pt, gl_w, gl_x, ov_dBdB, splp, kkn, v_1);
  // int B_iB_j/r^2 dr
  bsp::splineInt(nm2, k, glq_pt, gl_w, gl_x, ov_1_r2, spl, kkn, v_1_r2);

  auto v = std::unique_ptr<ModelV>(new V_1_r(z));

  // int B_i V(r) B_j dr
  bsp::splineInt(nm2, k, glq_pt, gl_w, gl_x, ov_V, spl, kkn, v);

  Enl.resize(lm1 * nm2);
  Cnl_tmp.resize(lm1 * nm22);
  aa.resize(lm1 * nk);
  w_bb.resize(lm1 * nk);

  for (int l = 0; l <= l_max; ++l) {
    std::copy(std::execution::seq, ov_BB.begin(), ov_BB.end(),
              w_bb.begin() + l * nk);
  }

  for (int l = 0; l <= l_max; ++l) {
    llp1 = l * (l + 1);
    lnk = l * nk;
    for (int ni = 0; ni < nm2; ++ni) {
      nik = ni * k;
      for (int j = 0; j < k; ++j) {
        aa[lnk + j + nik] = mass * ov_dBdB[j + nik] - ov_V[j + nik] +
                            mass * llp1 * ov_1_r2[j + nik];
      }
    }

    LAPACKE_dsbgvd(LAPACK_COL_MAJOR, 'V', 'U', nm2, k - 1, k - 1, &aa[lnk], k,
                   &w_bb[lnk], k, &Enl[l * nm2], &Cnl_tmp[l * nm22], nm2);
  }

  for (int l = 0; l <= l_max; ++l) {
    // Reshape with zeros at r=0 & r=R
    for (int ni = 0; ni < nm2; ++ni) {
      auto idxh = l * nm22 + ni * nm2;
      auto st_it = Cnl_tmp.begin() + idxh;
      auto end_it = st_it + nm2;
      // make the wf have a positive derivative at r=0
      auto val = 0.0;
      for (auto i = 0; i < k; ++i)
        val += Cnl_tmp[idxh + i] * spl[i + (k - 1) * k * glq_pt];

      if (std::signbit(val)) {
        std::transform(std::execution::seq, st_it, end_it, st_it,
                       [](auto &a) { return a * -1.0; });
      }
      std::copy(std::execution::seq, st_it, end_it,
                Cnl.begin() + l * n * nm2 + 1 + ni * n);
    }

    // Write hdf5 file
    writeHdf5(n, k, l, z, mass, pot, kkn, Enl, Cnl, outFile);
  }

  return 0;
}
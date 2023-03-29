#include "dmx1e.hpp"

int dmx1e::readConfig(std::string file, std::string &pot, int &qsz, char &gauge,
                      int &l_max) {
  YAML::Node settings = YAML::LoadFile(file);

  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "Core Potential:                              " << pot
            << std::endl;
  qsz = settings["Global_Settings"]["Outer_quadrature_size"].as<int>();
  std::cout << "No. of outer quadrature points:              " << qsz
            << std::endl;
  gauge = settings["Global_Settings"]["gauge"].as<char>();
  std::cout << "Gauge type ('l' length/'v' velocity):        " << gauge
            << std::endl;

  l_max = settings["Basis_Settings"]["l_max"].as<int>();
  std::cout << "Maximum one electron angular momentum:       " << l_max
            << std::endl;

  return 0;
}

int dmx1e::genDipole(std::string pot, int qsz, char gauge, int l_max) {
  int nen = 0; // no. of basis states
  int lp1 = 0;
  double kl, t_ab;
  std::string filename;

  H5::DataSpace cspace;
  H5::DataSpace memspace;
  auto D_set = std::make_unique<H5::DataSet>();
  hsize_t d_dim[2];

  filename = pot + "_w1e0.h5";
  auto file =
      std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  auto rset =
      std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Pr_o")));
  auto lc_sz = rset->getSpace().getSimpleExtentNpoints();

  nen = lc_sz / qsz;
  std::vector<double> wfn(lc_sz * (l_max + 1));
  std::vector<double> D(nen * nen);
  std::vector<double> qx(qsz), qw(qsz);
  d_dim[0] = nen;
  d_dim[1] = nen;

  rset->read(wfn.data(), H5::PredType::NATIVE_DOUBLE);
  rset = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Qr_o")));
  rset->read(qx.data(), H5::PredType::NATIVE_DOUBLE);
  rset = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Qw_o")));
  rset->read(qw.data(), H5::PredType::NATIVE_DOUBLE);

  // read wavefunctions for all l
  for (int l = 1; l <= l_max; ++l) {
    filename = pot + "_w1e" + std::to_string(l) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    rset =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Pr_o")));
    rset->read(&wfn[l * lc_sz], H5::PredType::NATIVE_DOUBLE);
    file->close();
  }

  switch (gauge) {
  case 'v': {
    std::vector<double> wfnp;
    wfnp.reserve(lc_sz);
#pragma omp parallel
    {
      for (int l = 0; l < l_max; ++l) {
#pragma omp single
        {
          filename = pot + "_w1e" + std::to_string(l + 1) + ".h5";
          file = std::make_unique<H5::H5File>(
              H5::H5File(filename, H5F_ACC_RDONLY));
          rset = std::make_unique<H5::DataSet>(
              H5::DataSet(file->openDataSet("Pr_p")));
          rset->read(wfnp.data(), H5::PredType::NATIVE_DOUBLE);
          file->close();
        }
        lp1 = l + 1;
        kl = sqrt((double)(lp1 * lp1) / (4.0 * lp1 * lp1 - 1));
        for (int n1 = 0; n1 < nen; ++n1) {
#pragma omp for private(t_ab)
          for (int n2 = 0; n2 < nen; ++n2) {
            t_ab = kl * dmx_int::tvelGL(qsz, lc_sz, n1, l, n2, lp1, qx, qw, wfn,
                                        wfnp);
            D[n2 + nen * n1] = t_ab;
          }
        }
#pragma omp single
        {
          filename =
              pot + std::to_string(l) + std::to_string(l + 1) + gauge + ".h5";
          file =
              std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_TRUNC));
          D_set = std::make_unique<H5::DataSet>(H5::DataSet(file->createDataSet(
              "d_if", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, d_dim))));
          D_set->write(&D[0], H5::PredType::NATIVE_DOUBLE);
          file->close();
        }
      }
    }
    wfnp.clear();
    break;
  }
  case 'l': {
#pragma omp parallel
    {
      for (int l = 0; l < l_max; ++l) {
        lp1 = l + 1;
        kl = sqrt((double)(lp1 * lp1) / (4.0 * lp1 * lp1 - 1));
        for (int n1 = 0; n1 < nen; ++n1) {
#pragma omp for private(t_ab)
          for (int n2 = 0; n2 < nen; ++n2) {
            t_ab =
                kl * dmx_int::tlenGL(qsz, lc_sz, n1, l, n2, lp1, qx, qw, wfn);
            D[n2 + nen * n1] = t_ab;
          }
        }
#pragma omp single
        {
          filename =
              pot + std::to_string(l) + std::to_string(l + 1) + gauge + ".h5";
          file =
              std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_TRUNC));
          D_set = std::make_unique<H5::DataSet>(H5::DataSet(file->createDataSet(
              "d_if", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, d_dim))));
          D_set->write(&D[0], H5::PredType::NATIVE_DOUBLE);
          file->close();
        }
      }
    }
    break;
  }
  }

  return 0;
}
#include "w2e.hpp"

int w2e::readConfig(std::string file, std::string &pot, char &gauge,
                    int &L_max) {
  YAML::Node settings = YAML::LoadFile(file);
  std::cout << "Global Settings:" << std::endl;
  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "  Core Potential:                       " << pot << std::endl;
  L_max = settings["Global_Settings"]["L_max"].as<int>();
  std::cout << "  max l:                                " << L_max << std::endl;
  gauge = settings["Global_Settings"]["gauge"].as<char>();
  std::cout << "  Gauge type ('l' length/'v' velocity): " << gauge << std::endl;

  return 0;
}

int w2e::formCIh0(std::string pot, int L_max, stvupt &vecs) {
  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> data = nullptr;
  std::unique_ptr<H5::DataSet> diag_set = nullptr;
  std::vector<double> ens, v12, eig;
  int L_full_sz, v_sz;

  hsize_t dimms1[1], dimms2[2];

  for (auto L = 0; L <= L_max; ++L) {
    // read sum energies
    filename = pot + "2_" + std::to_string(L) + "En.h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    data =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("e_2e")));
    L_full_sz = data->getSpace().getSimpleExtentNpoints();
    ens.resize(L_full_sz);
    data->read(ens.data(), H5::PredType::NATIVE_DOUBLE);
    file->close();

    eig.resize(L_full_sz);
    vecs[L]->resize(L_full_sz * L_full_sz);
    v_sz = L_full_sz * (L_full_sz + 1) / 2;
    v12.resize(v_sz);

    // read V_12 triangular format
    filename = pot + "V12_" + std::to_string(L) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    data =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("V_12")));
    data->read(v12.data(), H5::PredType::NATIVE_DOUBLE);
    file->close();

    for (auto i = 0; i < L_full_sz; ++i) {
      v12[(2 * L_full_sz - i - 1) * i / 2 + i] += ens[i];
      vecs[L]->at(i * L_full_sz + i) = v12[(2 * L_full_sz - i - 1) * i / 2 + i];
      for (auto j = i + 1; j < L_full_sz; ++j) {
        vecs[L]->at(i * L_full_sz + j) =
            v12[(2 * L_full_sz - i - 1) * i / 2 + j];
      }
    }

    LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'L', L_full_sz, vecs[L]->data(),
                   L_full_sz, eig.data());

    dimms1[0] = L_full_sz;
    dimms2[0] = L_full_sz;
    dimms2[1] = L_full_sz;

    filename = pot + "CI" + std::to_string(L) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_TRUNC));
    data = std::make_unique<H5::DataSet>(H5::DataSet(file->createDataSet(
        "En_CI", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, dimms1))));
    data->write(eig.data(), H5::PredType::NATIVE_DOUBLE);
    data = std::make_unique<H5::DataSet>(H5::DataSet(file->createDataSet(
        "CI_vecs", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms2))));
    data->write(vecs[L]->data(), H5::PredType::NATIVE_DOUBLE);
    file->close();

    if (L == 0) {
      std::cout << "Ground state energy (a.u.): " << eig[0] << "\n";
    }
  }

  return 0;
}

int w2e::formCIDipoles(std::string pot, char gauge, int L_max, stvupt &vecs) {
  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> dl = nullptr;
  H5::DataSpace dl_space;
  hsize_t dimms[2];
  int dip_sz;
  std::vector<double> dipole, temp;

  for (auto L = 0; L < L_max; ++L) {
    filename =
        pot + "2_" + std::to_string(L) + std::to_string(L + 1) + gauge + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    dl = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("d_if")));
    dl_space = dl->getSpace();
    dl_space.getSimpleExtentDims(dimms, NULL);
    dip_sz = dl_space.getSimpleExtentNpoints();
    dipole.resize(dip_sz);
    temp.resize(dip_sz);
    dl->read(dipole.data(), H5::PredType::NATIVE_DOUBLE);
    file->close();

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, dimms[1], dimms[0],
                dimms[0], 1.0, dipole.data(), dimms[1], vecs[L]->data(),
                dimms[0], 0.0, temp.data(), dimms[1]);

    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, dimms[1], dimms[0],
                dimms[1], 1.0, vecs[L + 1]->data(), dimms[1], temp.data(),
                dimms[1], 0.0, dipole.data(), dimms[1]);

    filename =
        pot + "CI_" + std::to_string(L) + std::to_string(L + 1) + gauge + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_TRUNC));
    dl = std::make_unique<H5::DataSet>(H5::DataSet(file->createDataSet(
        "CI_dip", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms))));
    dl->write(dipole.data(), H5::PredType::NATIVE_DOUBLE);
    file->close();
  }

  return 0;
}
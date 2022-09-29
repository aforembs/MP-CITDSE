#include "conv.hpp"

int conv::readConfig(std::string file, std::string &pot, int &L_max,
                     char &gauge, std::vector<int> &state_sz) {

  YAML::Node settings = YAML::LoadFile(file);
  std::cout << "Global Settings:"
            << "\n";
  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "  Core Potential:                       " << pot << "\n";
  L_max = settings["Global_Settings"]["L_max"].as<int>();
  std::cout << "  max L:                                " << L_max << "\n";
  gauge = settings["Global_Settings"]["gauge"].as<char>();
  std::cout << "  Gauge type ('l' length/'v' velocity): " << gauge << "\n";

  std::cout << "Propagator Settings:"
            << "\n";
  int states_L = settings["Propagator_Settings"]["states_in_l"].size();

  assert(states_L == L_max + 1);

  std::cout << "  No. of states in each L:  ";
  for (auto i = 0; i < states_L; ++i) {
    state_sz.push_back(
        settings["Propagator_Settings"]["states_in_l"][i].as<int>());
    std::cout << " " << state_sz[i];
  }
  std::cout << "\n";

  return 0;
}

int conv::calcEvecs(std::string pot, char gauge, int L_max,
                    std::vector<int> &state_sz, stvupt &vecs) {
  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> edata = nullptr, v12data = nullptr,
                               diag_set = nullptr;
  std::unique_ptr<H5::DataSet> L_set = nullptr;
  H5::DataSpace memspace, espace, vspace;
  std::vector<double> ens, v12, eig;
  std::vector<int> ifail;
  int L_sz, L_full_sz, v_sz;

  // hsize_t offset[] = {0}, stride[] = {1}, block[] = {1};
  // hsize_t count[1], dimms[1];

  for (auto L = 0; L < L_max; ++L) {
    L_sz = state_sz[L];
    // count[0] = L_sz;
    // dimms[0] = L_sz;

    vecs[L].get()->resize(L_sz * L_sz);
    // memspace.setExtentSimple(1, dimms, NULL);

    // read sum energies
    filename =
        pot + "2_" + std::to_string(L) + std::to_string(L + 1) + gauge + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    edata =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("e_i")));
    L_full_sz = edata->getSpace().getSimpleExtentNpoints();
    ens.resize(L_full_sz);
    edata->read(ens.data(), H5::PredType::NATIVE_DOUBLE);
    file->close();

    v_sz = L_full_sz * (L_full_sz + 1) / 2;
    v12.reserve(v_sz);

    // read V_12
    filename = pot + "V12_" + std::to_string(L) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    v12data =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("V_12")));
    v12data->read(v12.data(), H5::PredType::NATIVE_DOUBLE);
    file->close();

    for (auto i = 0; i < L_sz; ++i) {
      v12[(2 * L_full_sz - i - 1) * i / 2 + i] += ens[i];
      vecs[L].get()->at(i * L_sz + i) =
          v12[(2 * L_full_sz - i - 1) * i / 2 + i];
      for (auto j = i + 1; j < L_sz; ++j) {
        vecs[L].get()->at(i * L_sz + j) =
            v12[(2 * L_full_sz - i - 1) * i / 2 + j];
      }
    }

    eig.resize(L_sz);

    std::cout << LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', L_sz,
                                vecs[L].get()->data(), L_sz, eig.data())
              << " " << L << "\n";
  }

  L_sz = state_sz[L_max];
  // count[0] = L_sz;
  // dimms[0] = L_sz;

  vecs[L_max].get()->resize(L_sz * L_sz);
  // memspace.setExtentSimple(1, dimms, NULL);

  filename = pot + "2_" + std::to_string(L_max - 1) + std::to_string(L_max) +
             gauge + ".h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  edata = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("e_f")));
  L_full_sz = edata->getSpace().getSimpleExtentNpoints();
  ens.resize(L_full_sz);
  edata->read(ens.data(), H5::PredType::NATIVE_DOUBLE);
  file->close();

  v_sz = L_full_sz * (L_full_sz + 1) / 2;
  v12.reserve(v_sz);

  // read V_12
  filename = pot + "V12_" + std::to_string(L_max) + ".h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  v12data =
      std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("V_12")));
  v12data->read(v12.data(), H5::PredType::NATIVE_DOUBLE);
  file->close();

  for (auto i = 0; i < L_sz; ++i) {
    v12[(2 * L_full_sz - i - 1) * i / 2 + i] += ens[i];
    vecs[L_max].get()->at(i * L_sz + i) =
        v12[(2 * L_full_sz - i - 1) * i / 2 + i];
    for (auto j = i + 1; j < L_sz; ++j) {
      vecs[L_max].get()->at(i * L_sz + j) =
          v12[(2 * L_full_sz - i - 1) * i / 2 + j];
    }
  }

  eig.resize(L_sz);

  std::cout << LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', L_sz,
                              vecs[L_max].get()->data(), L_sz, eig.data())
            << " " << L_max << "\n";

  return 0;
}

int conv::readDipoles(std::string pot, char gauge, int L_max,
                      std::vector<int> &state_sz, stvupt &dipoles) {
  H5::DataSpace memspace;
  hsize_t offset[] = {0, 0}, stride[] = {1, 1}, block[] = {1, 1};
  hsize_t count[2], dimms[2];
  int L_sz, L1_sz;

  for (auto L = 0; L < L_max; ++L) {
    L_sz = state_sz[L];
    L1_sz = state_sz[L + 1];
    count[0] = L1_sz;
    count[1] = L_sz;
    dimms[0] = L1_sz;
    dimms[1] = L_sz;

    dipoles[L].get()->resize(L_sz * L1_sz);
    memspace.setExtentSimple(2, dimms, NULL);
    auto filename =
        pot + "2_" + std::to_string(L) + std::to_string(L + 1) + gauge + ".h5";
    auto file =
        std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    auto dl =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("d_if")));
    auto dl_space = dl->getSpace();
    dl_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    dl->read(dipoles[L].get()->data(), H5::PredType::NATIVE_DOUBLE, memspace,
             dl_space);
    file->close();
  }

  return 0;
}

int conv::transDip(std::string pot, char gauge, int L_max,
                   std::vector<int> &state_sz, stvupt &vecs, stvupt &dipoles) {
  std::vector<double> Dd, Res;
  // double alpha = 1.0, beta = 0.0;
  hsize_t D_dimms[2];

  for (auto L = 0; L < L_max; ++L) {
    auto L_sz = state_sz[L];
    auto L1_sz = state_sz[L + 1];
    Dd.resize(L_sz * L1_sz);
    Res.resize(L_sz * L1_sz);

    std::cout << "L: " << L << " dip.size(): " << dipoles[L].get()->size()
              << " vecs.size(): " << vecs[L].get()->size() << "\n";

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, L1_sz, L_sz, L_sz,
                1.0, dipoles[L].get()->data(), L_sz, vecs[L].get()->data(),
                L_sz, 0.0, Res.data(), L_sz);

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, L1_sz, L_sz, L1_sz,
                1.0, vecs[L + 1].get()->data(), L1_sz, Res.data(), L_sz, 0.0,
                Dd.data(), L_sz);

    D_dimms[0] = L1_sz;
    D_dimms[1] = L_sz;

    for (auto m = 0; m < L1_sz; ++m) {
      for (auto n = 0; n < L_sz; ++n) {
        Res.at(n * L1_sz + m) = Dd.at(m * L_sz + n);
      }
    }

    auto outfile_name = pot + "2_" + std::to_string(L) + std::to_string(L + 1) +
                        gauge + "_diag.h5";
    auto outfile =
        std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_TRUNC));
    auto dmx = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
        "d_if", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, D_dimms))));
    dmx->write(Res.data(), H5::PredType::NATIVE_DOUBLE);
    outfile->close();
  }

  return 0;
}
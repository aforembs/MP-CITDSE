#include "conv.hpp"

int conv::readConfig(std::string file, std::string &pot, char &gauge,
                     int &L_max, std::vector<int> &state_sz) {
  YAML::Node settings = YAML::LoadFile(file);
  std::cout << "Global Settings:" << std::endl;
  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "  Core Potential:                       " << pot << std::endl;
  L_max = settings["Global_Settings"]["L_max"].as<int>();
  std::cout << "  max l:                                " << L_max << std::endl;
  gauge = settings["Global_Settings"]["gauge"].as<char>();
  std::cout << "  Gauge type ('l' length/'v' velocity): " << gauge << std::endl;

  std::cout << "Propagator Settings:" << std::endl;
  int states_L = settings["Propagator_Settings"]["states_in_l"].size();

  assert(states_L == L_max + 1);

  std::cout << "  No. of states in each l:  ";
  for (auto i = 0; i < states_L; ++i) {
    state_sz.push_back(
        settings["Propagator_Settings"]["states_in_l"][i].as<int>());
    std::cout << " " << state_sz[i];
  }

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
  std::vector<double> ens, v12, eig, vecs;
  std::vector<int> ifail;
  std::vector<idx4> L_idx;
  int L_sz, L_full_sz, v_sz;

  hsize_t offset[] = {0}, stride[] = {1}, block[] = {1};
  hsize_t count[1], dimms[1];

  for (auto L = 0; L < L_max; ++L) {
    L_sz = state_sz[L];
    count[0] = L_sz;
    dimms[0] = L_sz;

    ens.reserve(L_sz);
    vecs[L].get()->resize(L_sz * L_sz);
    L_idx.reserve(L_sz);
    memspace.setExtentSimple(1, dimms, NULL);

    // read sum energies
    filename =
        pot + "2_" + std::to_string(L) + std::to_string(L + 1) + gauge + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    edata =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("e_i")));
    L_full_sz = edata->getSpace().getSimpleExtentNpoints();
    espace = edata->getSpace();
    espace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    edata->read(&ens[0], H5::PredType::NATIVE_DOUBLE, memspace, espace);
    file->close();

    v_sz = L_sz * (L_sz + 1) / 2;
    v12.reserve(v_sz);

    // read V_12
    filename = pot + "V12_" + std::to_string(L) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    v12data =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("V_12")));
    v12data->read(&v12[0], H5::PredType::NATIVE_DOUBLE);
    file->close();

    for (auto i = 0; i < L_sz; ++i) {
      v12[(2 * L_full_sz - i - 1) * i / 2 + i] += ens[i];
      vecs[L].get()->at(i * L_sz + i) = v12[(2 * L_full_sz - i - 1) * i / 2 + i];
      for (auto j = i + 1; j < L_sz; ++j) {
        vecs[L].get()->at(i * L_sz + j) = v12[(2 * L_full_sz - i - 1) * i / 2 + j];
      }
    }

    eig.resize(L_sz);

    std::cout << LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', L_sz, vecs[L]->data(),
                                L_sz, &eig[0])
              << "\n";
  }

  L_sz = state_sz[L_max];
  count[0] = L_sz;
  dimms[0] = L_sz;

  ens.reserve(L_sz);
  vecs[L_max].get()->resize(L_sz * L_sz);
  L_idx.reserve(L_sz);
  memspace.setExtentSimple(1, dimms, NULL);

  filename = pot + "2_" + std::to_string(L_max - 1) + std::to_string(L_max) +
             gauge + ".h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  edata = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("e_f")));
  L_full_sz = edata->getSpace().getSimpleExtentNpoints();
  espace = edata->getSpace();
  espace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
  edata->read(&ens[0], H5::PredType::NATIVE_DOUBLE, memspace, espace);
  file->close();

  v_sz = L_sz * (L_sz + 1) / 2;
  v12.reserve(L_sz * L_sz);

  // read V_12
  filename = pot + "V12_" + std::to_string(L_max) + ".h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  v12data =
      std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("V_12")));
  v12data->read(&v12[0], H5::PredType::NATIVE_DOUBLE);
  file->close();

    for (auto i = 0; i < L_sz; ++i) {
      v12[(2 * L_full_sz - i - 1) * i / 2 + i] += ens[i];
      vecs[L_max].get()->at(i * L_sz + i) = v12[(2 * L_full_sz - i - 1) * i / 2 + i];
      for (auto j = i + 1; j < L_sz; ++j) {
        vecs[L_max].get()->at(i * L_sz + j) = v12[(2 * L_full_sz - i - 1) * i / 2 + j];
      }
    }

  eig.resize(L_sz);

  std::cout << LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', L_sz, vecs[L_max]->data(), L_sz,
                              &eig[0])
            << "\n";

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

int conv::transDip(std::vector<int> &state_sz, stvupt &vecs, stvupt &dipoles) {

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, state_sz[L], );

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,);

  return 0;
}
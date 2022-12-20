#include "pes2e.hpp"

int pes2e::readConfig(std::string file, std::string &pot, int &L_max,
                      int &l_max, std::vector<int> &state_sz) {
  YAML::Node settings = YAML::LoadFile(file);
  std::cout << "Global Settings:" << std::endl;
  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "  Core Potential:                       " << pot << std::endl;
  L_max = settings["Global_Settings"]["L_max"].as<int>();
  std::cout << "  max L:                                " << L_max << std::endl;

  std::cout << "Basis Settings:" << std::endl;
  l_max = settings["Basis_Settings"]["l_max"].as<int>();
  std::cout << "  max l:                                " << l_max << std::endl;

  std::cout << "Propagator Settings:" << std::endl;
  int states_L = settings["Propagator_Settings"]["states_in_l"].size();

  assert(states_L == L_max + 1);

  std::cout << "  No. of states in each l: ";
  for (auto i = 0; i < states_L; ++i) {
    state_sz.push_back(
        settings["Propagator_Settings"]["states_in_l"][i].as<int>());
    std::cout << " " << state_sz[i];
  }
  std::cout << std::endl;

  return 0;
}

int pes2e::readCt(std::string file, std::vector<std::complex<double>> &ct) {
  std::ifstream fl(file);
  std::string temp;
  while (std::getline(fl, temp)) {
    std::istringstream iss(temp);
    int idx;
    double real, imag;
    iss >> idx >> real >> imag;
    ct.push_back(std::complex<double>(real, imag));
  }

  return 0;
}

int pes2e::genPES2eb(std::string pot, int L_max, std::vector<int> &state_sz,
                     std::vector<std::complex<double>> &ct) {
  std::vector<idx4> ct_idx(ct.size());
  std::vector<int> offs;
  auto sum = 0;
  offs.push_back(sum);
  int L_sz;

  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> e_set = nullptr, L_set = nullptr;
  hsize_t offset[1] = {0}, stride[1] = {1}, block[1] = {1};
  hsize_t count[1], dimms[1];
  H5::DataSpace memspace, e_space, L_space;

  std::vector<double> en_2e(ct.size());

  // Read 2e energies
  for (auto L = 0; L <= L_max; ++L) {
    L_sz = state_sz[L];
    count[0] = L_sz;
    dimms[0] = L_sz;
    sum += L_sz;
    offs.push_back(sum);
    memspace.setExtentSimple(1, dimms, NULL);

    filename = pot + "2_" + std::to_string(L) + "En.h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    e_set =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("e_2e")));
    e_space = e_set->getSpace();
    e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    e_set->read(&en_2e[offs[L]], H5::PredType::NATIVE_DOUBLE, memspace,
                e_space);
    file->close();
  }

  // Read indices
  for (auto L = 0; L <= L_max; ++L) {
    L_sz = state_sz[L];
    count[0] = L_sz * 4;
    dimms[0] = L_sz * 4;
    memspace.setExtentSimple(1, dimms, NULL);

    filename = pot + "2_" + std::to_string(L) + "En.h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    L_set =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("idx")));
    L_space = L_set->getSpace();
    L_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    L_set->read(&ct_idx[offs[L]], H5::PredType::NATIVE_INT32, memspace,
                L_space);
    file->close();
  }

  std::vector<double> ens, v12, eig, vecs;
  std::vector<std::complex<double>> cvecs, c_proj;
  int L_full_sz, v_sz;
  double ion_yield = 0.0;
  // read sum energies
  int off = 0;
  for (auto L = 0; L <= L_max; ++L) {
    L_sz = state_sz[L];
    vecs.resize(L_sz * L_sz);
    cvecs.resize(L_sz * L_sz);
    c_proj.resize(L_sz);
    filename = pot + "2_" + std::to_string(L) + "En.h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    auto en_data =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("e_2e")));
    L_full_sz = en_data->getSpace().getSimpleExtentNpoints();
    ens.resize(L_full_sz);
    en_data->read(ens.data(), H5::PredType::NATIVE_DOUBLE);
    file->close();

    v_sz = L_full_sz * (L_full_sz + 1) / 2;
    v12.reserve(v_sz);

    // read V_12
    filename = pot + "V12_" + std::to_string(L) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    auto v12data =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("V_12")));
    v12data->read(v12.data(), H5::PredType::NATIVE_DOUBLE);
    file->close();

    for (auto i = 0; i < L_sz; ++i) {
      v12[(2 * L_full_sz - i - 1) * i / 2 + i] += ens[i];
      vecs.at(i * L_sz + i) = v12[(2 * L_full_sz - i - 1) * i / 2 + i];
      for (auto j = i + 1; j < L_sz; ++j) {
        vecs.at(i * L_sz + j) = v12[(2 * L_full_sz - i - 1) * i / 2 + j];
      }
    }

    eig.resize(L_sz);
    LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'N', 'U', L_sz, vecs.data(), L_sz,
                   eig.data());

    std::ofstream outfile(pot + "_pes" + std::to_string(L) + ".dat",
                          std::ios::out);
    for (auto i = 1; i < L_sz - 1; ++i) {
      outfile << std::setprecision(16) << eig[i] << " "
              << std::norm(
                     ct[off + i]) // * 2.0 / std::abs(eig[i + 1] - eig[i - 1])
              << "\n";
      if (eig[i] + 2.0 > 0.0) {
        ion_yield += std::norm(ct[off + i]);
      }
    }

    if (L == 0) {
      std::cout << std::setprecision(16)
                << " ground state pop: " << std::norm(ct[0]) << "\n";
    }
    outfile.close();
    off += L_sz;
  }

  // Get the norm of c(t)
  double nrm = 0;
  for (auto &v : ct) {
    nrm += std::norm(v);
  }

  // Print the ground population and the norm of c(t)
  std::cout << std::setprecision(16) << " norm: " << nrm
            << " yield: " << ion_yield << "\n";

  return 0;
}
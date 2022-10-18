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

int pes2e::genPES(std::string pot, std::string dir, int L_max, int l_max,
                  std::vector<int> &state_sz,
                  std::vector<std::complex<double>> &ct) {
  std::vector<idx4> ct_idx(ct.size());
  std::vector<int> offs;
  auto sum = 0;
  offs.push_back(sum);
  int L_sz;

  int ncf, sym;
  std::vector<cfg::line> cfgs;
  auto max_N = 0, max_nl = 0;
  cfg::line max_n2l;

  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> e_set = nullptr, L_set = nullptr;
  hsize_t offset[1] = {0}, stride[1] = {1}, block[1] = {1};
  hsize_t count[1], dimms[1];
  H5::DataSpace memspace, e_space, L_space;

  for (auto L = 0; L <= L_max; ++L) {
    cfg::readCfg(dir, L, sym, ncf, cfgs);

    max_n2l = *std::max_element(cfgs.begin(), cfgs.end(),
                                [](cfg::line const &a, cfg::line const &b) {
                                  return a.n2max < b.n2max;
                                });
    max_N = std::max(max_n2l.n2max, max_N);
    max_nl = std::max(ncf, max_nl);
  }

  std::vector<double> en_1e((l_max + 1) * max_N);

  count[0] = max_N;
  dimms[0] = max_N;
  memspace.setExtentSimple(1, dimms, NULL);

  // Read 1e energies
  for (auto l = 0; l <= l_max; ++l) {
    filename = pot + std::to_string(l) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    e_set = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
    e_space = e_set->getSpace();
    e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    e_set->read(&en_1e[l * max_N], H5::PredType::NATIVE_DOUBLE, memspace,
                e_space);
    file->close();
  }

  // Read indices
  for (auto L = 0; L <= L_max; ++L) {
    L_sz = state_sz[L];
    sum += L_sz;
    offs.push_back(sum);
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

  std::vector<double> PES_En;
  // int st_idx = l_max * max_N;
  for (auto i = 0; i < static_cast<int>(en_1e.size()); ++i) {
    PES_En.push_back(en_1e[i]);
  }

  std::sort(std::execution::par_unseq, PES_En.begin(), PES_En.end());
  PES_En.erase(
      std::unique(std::execution::par_unseq, PES_En.begin(), PES_En.end()),
      PES_En.end());
  std::cout << PES_En.size() << "\n";
  std::vector<double> PES(PES_En.size());

  // Loop for calculating dP/dEk = 2*|c_i|^2/|E_{i+1}-E_{i-1}|
  for (auto L = 0; L <= 0; ++L) {
    for (auto i = 1; i < state_sz[L]; ++i) {
      for (auto k = 1; k < static_cast<int>(PES_En.size()) - 1; ++k) {
        auto ndiff = std::abs(PES_En[k] - PES_En[k - 1]);
        auto pdiff = std::abs(PES_En[k + 1] - PES_En[k]);
        auto l1_idx = ct_idx[offs[L] + i].l1 * max_N;
        auto n1_idx = ct_idx[offs[L] + i].n1;

        bool cond1 = (en_1e[l1_idx + n1_idx] > PES_En[k] - 0.5 * ndiff) &&
                     (en_1e[l1_idx + n1_idx] < PES_En[k] + 0.5 * pdiff);
        if (cond1) {
          PES[k] += std::norm(ct[offs[L] + i]) * 2.0 / (ndiff + pdiff);
        }

        auto l2_idx = ct_idx[offs[L] + i].l2 * max_N;
        auto n2_idx = ct_idx[offs[L] + i].n2;

        bool cond2 = (en_1e[l2_idx + n2_idx] > PES_En[k] - 0.5 * ndiff) &&
                     (en_1e[l2_idx + n2_idx] < PES_En[k] + 0.5 * pdiff);
        if (cond2) {
          PES[k] += std::norm(ct[offs[L] + i]) * 2.0 / (ndiff + pdiff);
        }
      }
    }
  }

  // Write PES for energies > 0
  std::fstream outfile(pot + "_pes.dat", std::ios::out);
  for (auto i = 0; i < static_cast<int>(PES_En.size()) - 1; ++i) {
    if (PES_En[i] > 0.0) {
      outfile << PES_En[i] << " " << PES[i] << "\n";
    }
  }
  outfile.close();

  // Get the norm of c(t)
  double nrm = 0;
  for (auto &v : ct) {
    nrm += std::norm(v);
  }

  // Print the ground population and the norm of c(t)
  std::cout << std::setprecision(15) << " norm: " << nrm << "\n";

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

  std::vector<std::complex<double>> cg(state_sz[0]);
  // Read Ground State
  std::ifstream fl(pot + "_c0.dat");
  std::string temp;
  for (auto i = 0; i < state_sz[0]; ++i) {
    std::getline(fl, temp);
    std::istringstream iss(temp);
    int idx;
    double real;
    iss >> idx >> real;
    cg[i] = std::complex<double>(real, 0.0);
  }

  std::vector<double> ens, v12, eig, vecs;
  std::vector<std::complex<double>> cvecs, c_proj;
  int L_full_sz, v_sz;
  // read sum energies
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
    LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', L_sz, vecs.data(), L_sz,
                   eig.data());

    std::complex<double> alp(1.0, 0.0);
    std::complex<double> bet(0.0, 0.0);

    for (auto i = 0; i < static_cast<int>(vecs.size()); ++i) {
      cvecs[i] = std::complex<double>(vecs[i], 0.0);
    }

    cblas_zgemv(CblasRowMajor, CblasTrans, L_sz, L_sz,
                reinterpret_cast<double *>(&alp),
                reinterpret_cast<double *>(cvecs.data()), L_sz,
                reinterpret_cast<double *>(&ct[offs[L]]), 1,
                reinterpret_cast<double *>(&bet),
                reinterpret_cast<double *>(c_proj.data()), 1);

    std::ofstream outfile(pot + "_pes" + std::to_string(L) + ".dat",
                          std::ios::out);
    for (auto i = 0; i < L_sz; ++i) {
      outfile << std::setprecision(16) << eig[i] << " " << std::norm(c_proj[i])
              << "\n";
    }

    if (L == 0) {
      std::cout << std::setprecision(16)
                << " ground state pop: " << std::norm(c_proj[0]) << "\n";
    }
    outfile.close();
  }

  // Get the norm of c(t)
  double nrm = 0;
  for (auto &v : ct) {
    nrm += std::norm(v);
  }

  // Print the ground population and the norm of c(t)
  std::cout << std::setprecision(16) << " norm: " << nrm << "\n";

  return 0;
}
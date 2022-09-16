#include "td_read.hpp"

int tdrd::readConfig(std::string file, std::string &pot, char &gauge,
                     int &L_max, std::vector<int> &state_sz, int &Lanc_iter,
                     int &num_eval, double &timestep, double &w, double &Io,
                     double &cepd, int &cycles) {

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
  std::cout << std::endl;
  Lanc_iter = settings["Propagator_Settings"]["Lanczos_iterations"].as<int>();
  std::cout << "  No. of Lanczos Iterations:            " << Lanc_iter
            << std::endl;
  num_eval = settings["Propagator_Settings"]["min_eigenvals"].as<int>();
  std::cout << "  No. of Lowest Eigenvalues:            " << num_eval
            << std::endl;
  timestep = settings["Propagator_Settings"]["dt"].as<double>();
  std::cout << "  timestep dt:                          " << timestep
            << std::endl;

  std::cout << "Field Parameters:" << std::endl;
  w = settings["Field_Parameters"]["w"].as<double>();
  std::cout << "  w:                                    " << w << std::endl;
  Io = settings["Field_Parameters"]["Io"].as<double>();
  std::cout << "  Io:                                   " << Io << std::endl;
  cepd = settings["Field_Parameters"]["cepd"].as<double>();
  std::cout << "  cepd:                                 " << cepd << std::endl;
  cycles = settings["Field_Parameters"]["cycles"].as<int>();
  std::cout << "  cycles:                               " << cycles
            << std::endl;

  return 0;
}

template <class Container, class ElementType = typename Container::value_type>
constexpr auto element_of(const Container &, ElementType v = 0) {
  return v;
}

int tdrd::readStructure(std::string pot, char gauge, int L_max, int &ct_sz,
                        std::vector<int> &state_sz, std::vector<int> &offs,
                        stvupt &blocks) {
  std::string filename;
  std::unique_ptr<H5::H5File> efile = nullptr, file = nullptr;
  std::unique_ptr<H5::DataSet> edata = nullptr, v12data = nullptr;
  std::unique_ptr<H5::DataSet> diag_set = nullptr;
  std::unique_ptr<H5::DataSet> L_set = nullptr;
  H5::DataSpace memspace, espace, vspace, lspace;
  std::vector<double> ens, v12;
  std::vector<int> ifail;
  std::vector<idx4> L_idx;
  int L_sz, L_full_sz, v_sz;

  hsize_t offset1[] = {0}, stride1[] = {1}, block1[] = {1};
  hsize_t count1[1], dimms1[1];

  ct_sz = std::accumulate(state_sz.begin(), state_sz.end(),
                          element_of(state_sz, 0));
  auto sum = 0;
  offs.push_back(sum);

  for (auto L = 0; L < L_max; ++L) {
    L_sz = state_sz[L];
    count1[0] = L_sz;
    dimms1[0] = L_sz;
    sum += L_sz;
    offs.push_back(sum);

    ens.reserve(L_sz);
    blocks[L].get()->resize(L_sz * L_sz);
    L_idx.reserve(L_sz);
    memspace.setExtentSimple(1, dimms1, NULL);

    // read sum energies
    filename =
        pot + "2_" + std::to_string(L) + std::to_string(L + 1) + gauge + ".h5";
    efile = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    edata =
        std::make_unique<H5::DataSet>(H5::DataSet(efile->openDataSet("e_i")));
    L_full_sz = edata->getSpace().getSimpleExtentNpoints();
    espace = edata->getSpace();
    espace.selectHyperslab(H5S_SELECT_SET, count1, offset1, stride1, block1);
    edata->read(&ens[0], H5::PredType::NATIVE_DOUBLE, memspace, espace);
    efile->close();

    v_sz = L_full_sz * (L_full_sz + 1) / 2;
    v12.reserve(v_sz);

    // read V_12 triangular format
    filename = pot + "V12_" + std::to_string(L) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    v12data =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("V_12")));
    v12data->read(&v12[0], H5::PredType::NATIVE_DOUBLE);
    file->close();

    // get idx
    count1[0] = L_sz * 4;
    dimms1[0] = L_sz * 4;
    memspace.setExtentSimple(1, dimms1, NULL);
    filename = pot + std::to_string(L) + "idx.h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    L_set =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("idx")));
    lspace = L_set->getSpace();
    lspace.selectHyperslab(H5S_SELECT_SET, count1, offset1, stride1, block1);
    L_set->read(&L_idx[0], H5::PredType::NATIVE_INT32, memspace, lspace);
    file->close();

    // #pragma omp parallel
    {
      for (auto i = 0; i < L_sz; ++i) {
        v12[(2 * L_full_sz - i - 1) * i / 2 + i] += ens[i];
        blocks[L].get()->at(i * L_sz + i) =
            v12[(2 * L_full_sz - i - 1) * i / 2 + i];
        // #pragma omp for
        for (auto j = i + 1; j < L_sz; ++j) {
          blocks[L].get()->at(i * L_sz + j) =
              v12[(2 * L_full_sz - i - 1) * i / 2 + j];
        }
      }
    }
    v12.clear();
  }

  L_sz = state_sz[L_max];
  count1[0] = L_sz;
  dimms1[0] = L_sz;
  sum += L_sz;
  offs.push_back(sum);

  ens.reserve(L_sz);
  blocks[L_max].get()->resize(L_sz * L_sz);
  L_idx.reserve(L_sz);
  memspace.setExtentSimple(1, dimms1, NULL);

  filename = pot + "2_" + std::to_string(L_max - 1) + std::to_string(L_max) +
             gauge + ".h5";
  efile = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  edata = std::make_unique<H5::DataSet>(H5::DataSet(efile->openDataSet("e_f")));
  L_full_sz = edata->getSpace().getSimpleExtentNpoints();
  espace = edata->getSpace();
  espace.selectHyperslab(H5S_SELECT_SET, count1, offset1, stride1, block1);
  edata->read(&ens[0], H5::PredType::NATIVE_DOUBLE, memspace, espace);
  efile->close();

  v_sz = L_full_sz * (L_full_sz + 1) / 2;
  v12.reserve(v_sz);

  // read V_12
  filename = pot + "V12_" + std::to_string(L_max) + ".h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  v12data =
      std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("V_12")));
  v12data->read(&v12[0], H5::PredType::NATIVE_DOUBLE);
  file->close();

  // get idx
  count1[0] = L_sz * 4;
  dimms1[0] = L_sz * 4;
  memspace.setExtentSimple(1, dimms1, NULL);
  filename = pot + std::to_string(L_max) + "idx.h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  L_set = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("idx")));
  lspace = L_set->getSpace();
  lspace.selectHyperslab(H5S_SELECT_SET, count1, offset1, stride1, block1);
  L_set->read(&L_idx[0], H5::PredType::NATIVE_INT32, memspace, lspace);
  file->close();

  // #pragma omp parallel
  {
    for (auto i = 0; i < L_sz; ++i) {
      v12[(2 * L_full_sz - i - 1) * i / 2 + i] += ens[i];
      blocks[L_max].get()->at(i * L_sz + i) =
          v12[(2 * L_full_sz - i - 1) * i / 2 + i];
      // #pragma omp for
      for (auto j = i + 1; j < L_sz; ++j) {
        blocks[L_max].get()->at(i * L_sz + j) =
            v12[(2 * L_full_sz - i - 1) * i / 2 + j];
      }
    }
  }
  v12.clear();

  return 0;
}

int tdrd::readDipoles(std::string pot, char gauge, int L_max,
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
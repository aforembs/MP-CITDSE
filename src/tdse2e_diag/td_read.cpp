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

int tdrd::readEns(std::string pot, int L_max, int &ct_sz,
                  std::vector<int> &state_sz, std::vector<int> &offs,
                  stvupt &ens) {
  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> edata = nullptr;
  H5::DataSpace memspace, espace;
  std::vector<idx4> L_idx;
  int L_sz;

  hsize_t offset[] = {0}, stride[] = {1}, block[] = {1};
  hsize_t count[1], dimms[1];

  ct_sz = std::accumulate(state_sz.begin(), state_sz.end(),
                          element_of(state_sz, 0));
  auto sum = 0;
  offs.push_back(sum);

  for (auto L = 0; L <= L_max; ++L) {
    L_sz = state_sz[L];
    count[0] = L_sz;
    dimms[0] = L_sz;
    sum += L_sz;
    offs.push_back(sum);

    ens[L].get()->resize(L_sz);
    memspace.setExtentSimple(1, dimms, NULL);

    // read sum energies
    filename = pot + "_diag" + std::to_string(L) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    edata = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
    espace = edata->getSpace();
    espace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    edata->read(ens[L].get()->data(), H5::PredType::NATIVE_DOUBLE, memspace,
                espace);
    file->close();
  }

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
    auto filename = pot + "2_" + std::to_string(L) + std::to_string(L + 1) +
                    gauge + "_diag.h5";
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

int tdrd::readGrCt(std::string pot, std::vector<int> &state_sz,
                   std::vector<std::complex<double>> &ct) {
  std::ifstream fl(pot + "_c0.dat");
  std::string temp;
  for (auto i = 0; i < state_sz[0]; ++i) {
    std::getline(fl, temp);
    std::istringstream iss(temp);
    int idx;
    double real;
    iss >> idx >> real;
    ct[i] = std::complex<double>(real, 0.0);
  }

  return 0;
}
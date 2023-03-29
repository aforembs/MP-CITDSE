#include "td_read.hpp"

int tdrd::readConfig(std::string file, std::string &pot, std::string set_base,
                     std::string option, int &L_max, char &gauge,
                     std::vector<int> &state_sz, double &timestep,
                     std::string &shape, double &w, double &Io, double &cepd,
                     int &cycles) {

  YAML::Node settings = YAML::LoadFile(file);
  std::cout << "Global Settings:" << std::endl;
  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "  Core Potential:                       " << pot << std::endl;
  L_max = settings[set_base][option].as<int>();
  std::cout << "  max l/L:                              " << L_max << std::endl;
  gauge = settings["Global_Settings"]["gauge"].as<char>();
  std::cout << "  Gauge type ('l' length/'v' velocity): " << gauge << std::endl;

  std::cout << "Propagator Settings:" << std::endl;
  std::vector<int> loc_sz =
      settings["Propagator_Settings"]["states_in_l"].as<std::vector<int>>();
  assert(static_cast<int>(loc_sz.size()) == L_max + 1);
  std::cout << "  states per l/L:          ";
  for (auto i = 0; i < static_cast<int>(loc_sz.size()); ++i) {
    std::cout << loc_sz[i] << " ";
    state_sz.push_back(loc_sz[i]);
  }
  loc_sz.clear();

  std::cout << std::endl;
  timestep = settings["Propagator_Settings"]["dt"].as<double>();
  std::cout << "  timestep dt:                          " << timestep
            << std::endl;

  std::cout << "Field Parameters:" << std::endl;
  shape = settings["Field_Parameters"]["shape"].as<std::string>();
  std::cout << "  envelope shape:                       " << shape << std::endl;
  w = settings["Field_Parameters"]["w"].as<double>();
  std::cout << "  photon energy:                        " << w << std::endl;
  Io = settings["Field_Parameters"]["Io"].as<double>();
  std::cout << "  Peak intensity:                       " << Io << std::endl;
  cepd = settings["Field_Parameters"]["cepd"].as<double>();
  std::cout << "  phase shift:                          " << cepd << std::endl;
  cycles = settings["Field_Parameters"]["cycles"].as<int>();
  std::cout << "  number of cycles:                     " << cycles
            << std::endl;

  return 0;
}

int tdrd::readEnergies(std::string pot, std::string setname, int L_max,
                       int &ct_sz, std::vector<int> &state_sz,
                       std::vector<int> &offs, stvupt &eig) {
  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> edata = nullptr;
  H5::DataSpace e_space, memspace;
  int L_sz;
  hsize_t offset[] = {0}, stride[] = {1}, block[] = {1};
  hsize_t count[1], dimms[1];

  ct_sz = std::reduce(state_sz.begin(), state_sz.end(), 0);
  auto sum = 0;
  offs.push_back(sum);

  for (auto L = 0; L <= L_max; ++L) {
    L_sz = state_sz[L];
    sum += L_sz;
    offs.push_back(sum);
    count[0] = L_sz;
    dimms[0] = L_sz;
    memspace.setExtentSimple(1, dimms, NULL);

    eig[L]->resize(L_sz);

    filename = pot + std::to_string(L) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    edata =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet(setname)));
    e_space = edata->getSpace();
    e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    edata->read(eig[L]->data(), H5::PredType::NATIVE_DOUBLE, memspace, e_space);
    file->close();
  }

  return 0;
}

int tdrd::readDipoles(std::string pot, std::string setname, char gauge,
                      int L_max, std::vector<int> &state_sz, stvupt &dipoles) {
  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> dl = nullptr;
  H5::DataSpace dl_space, memspace;
  int L_sz, L1_sz;
  hsize_t offset[] = {0, 0}, stride[] = {1, 1}, block[] = {1, 1};
  hsize_t count[2], dimms[2];

  for (auto L = 0; L < L_max; ++L) {
    L_sz = state_sz[L];
    L1_sz = state_sz[L + 1];
    count[0] = L_sz;
    count[1] = L1_sz;
    dimms[0] = L_sz;
    dimms[1] = L1_sz;
    memspace.setExtentSimple(2, dimms, NULL);

    dipoles[L]->resize(L_sz * L1_sz);

    filename = pot + std::to_string(L) + std::to_string(L + 1) + gauge + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    dl = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet(setname)));
    dl_space = dl->getSpace();
    dl_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    dl->read(dipoles[L]->data(), H5::PredType::NATIVE_DOUBLE, memspace,
             dl_space);
    file->close();
  }

  return 0;
}

int tdrd::readInitCt(std::string file, int ct_sz, std::vector<double> &ct) {
  std::ifstream fl(file);
  std::string temp;
  int i = 0;
  while (i < ct_sz * 2) {
    std::getline(fl, temp);
    temp = std::regex_replace(temp, std::regex("^ +"), "");

    if (temp[0] == '#') {
      continue;
    }

    std::istringstream iss(temp);
    int idx;
    double real, complex;
    iss >> idx >> real >> complex;
    ct[i] = real;
    ct[i + 1] = complex;
    i += 2;
  }

  return 0;
}
#include "td_read.hpp"

int tdrd::readConfig(std::string file, std::string &pot, char &gauge,
                     int &l_max, std::vector<int> &state_sz, double &timestep,
                     double &w, double &Io, double &cepd, int &cycles) {

  YAML::Node settings = YAML::LoadFile(file);
  std::cout << "Global Settings:" << std::endl;
  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "  Core Potential:                       " << pot << std::endl;
  gauge = settings["Global_Settings"]["gauge"].as<char>();
  std::cout << "  Gauge type ('l' length/'v' velocity): " << gauge << std::endl;

  std::cout << "Basis Settings:" << std::endl;
  l_max = settings["Basis_Settings"]["l_max"].as<int>();
  std::cout << "  max l:                                " << l_max << std::endl;

  std::cout << "Propagator Settings:" << std::endl;
  int states_l = settings["Propagator_Settings"]["states_in_l"].size();

  assert(states_l == l_max + 1);

  std::cout << "  No. of states in each l:  ";
  for (auto i = 0; i < states_l; ++i) {
    state_sz.push_back(
        settings["Propagator_Settings"]["states_in_l"][i].as<int>());
    std::cout << " " << state_sz[i];
  }
  std::cout << std::endl;
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

int tdrd::readEnergies(std::string pot, int l_max, std::vector<int> &state_sz,
                       std::vector<double> &eps, std::vector<int> &offs,
                       int &eps_sz, std::vector<double *> &eig) {
  H5::DataSpace memspace;
  hsize_t offset[1], count[1], stride[1], block[1], dimms[1];
  offset[0] = 0;
  stride[0] = 1;
  block[0] = 1;

  eps_sz = std::reduce(state_sz.begin(), state_sz.end(), 0);
  eps.reserve(eps_sz);
  eig.push_back(eps.data());
  auto sum = 0;
  offs.push_back(sum);

  for (auto l = 0; l <= l_max; ++l) {
    count[0] = state_sz[l];
    dimms[0] = state_sz[l];
    memspace.setExtentSimple(1, dimms, NULL);
    auto filename = pot + std::to_string(l) + ".h5";
    auto file =
        std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    auto ei =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
    auto ei_space = ei->getSpace();
    ei_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    ei->read(eig[l], H5::PredType::NATIVE_DOUBLE, memspace, ei_space);
    file->close();
    sum += state_sz[l];
    offs.push_back(sum);
    eig.push_back(&eps[sum]);
  }

  return 0;
}

int tdrd::readDipoles(std::string pot, char gauge, int l_max,
                      std::vector<int> &state_sz, stvupt &dipoles) {
  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> dl = nullptr;
  H5::DataSpace dl_space, memspace;
  int l_sz, l1_sz;
  hsize_t offset[] = {0, 0}, stride[] = {1, 1}, block[] = {1, 1};
  hsize_t count[2], dimms[2];

  for (auto l = 0; l < l_max; ++l) {
    l_sz = state_sz[l];
    l1_sz = state_sz[l + 1];
    count[0] = l1_sz;
    count[1] = l_sz;
    dimms[0] = l1_sz;
    dimms[1] = l_sz;
    memspace.setExtentSimple(2, dimms, NULL);

    dipoles[l]->resize(l_sz * l1_sz);

    filename = pot + std::to_string(l) + std::to_string(l + 1) + gauge + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    dl = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("d_if")));
    dl_space = dl->getSpace();
    dl_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    dl->read(dipoles[l]->data(), H5::PredType::NATIVE_DOUBLE, memspace,
             dl_space);
    file->close();
  }

  return 0;
}
#include "gr_read.hpp"

int grrd::readConfig(std::string file, std::string &pot, char &gauge,
                     int &L0_sz, double &dt) {
  YAML::Node settings = YAML::LoadFile(file);
  std::cout << "Global Settings:"
            << "\n";
  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "  Core Potential:                       " << pot << "\n";
  gauge = settings["Global_Settings"]["gauge"].as<char>();
  std::cout << "  Gauge type ('l' length/'v' velocity): " << gauge << "\n";

  std::cout << "Propagator Settings:"
            << "\n";
  L0_sz = settings["Propagator_Settings"]["states_in_l"][0].as<int>();
  std::cout << "  No. of states in L = 0:               " << L0_sz << "\n";
  dt = settings["Propagator_Settings"]["dt"].as<double>();
  std::cout << "  dt:                                   " << dt << "\n";

  return 0;
}

int grrd::readStructure(std::string pot, int L0_sz,
                        std::vector<double> &ens, std::vector<double> &block) {
  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> edata = nullptr, v12data = nullptr;
  std::unique_ptr<H5::DataSet> L_set = nullptr;
  H5::DataSpace memspace, espace, lspace;
  std::vector<double> v12;
  std::vector<idx4> L_idx;
  int L_sz, v_sz;
  constexpr int L = 0;

  hsize_t offset[] = {0}, stride[] = {1}, hblock[] = {1};
  hsize_t count[] = {static_cast<hsize_t>(L0_sz * 4)},
          dimms[] = {static_cast<hsize_t>(L0_sz * 4)};
  memspace.setExtentSimple(1, dimms, NULL);

  L_idx.resize(L0_sz);
  block.resize(L0_sz * L0_sz);

  // get idx
  filename = pot + "2_" + std::to_string(L) + "En.h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  edata =
      std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("e_2e")));
  L_sz = edata->getSpace().getSimpleExtentNpoints();
  ens.resize(L_sz);
  espace = edata->getSpace();
  edata->read(ens.data(), H5::PredType::NATIVE_DOUBLE); 
  L_set =
      std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("idx")));
  lspace = L_set->getSpace();
  lspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, hblock);
  L_set->read(&L_idx[0], H5::PredType::NATIVE_INT32, memspace, lspace);
  file->close();

  v_sz = L_sz * (L_sz + 1) / 2;
  v12.resize(v_sz);

  // read V_12 triangular format
  filename = pot + "V12_" + std::to_string(L) + ".h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  v12data =
      std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("V_12")));
  v12data->read(v12.data(), H5::PredType::NATIVE_DOUBLE);
  file->close();

  // #pragma omp parallel
  {
    for (auto i = 0; i < L0_sz; ++i) {
      v12[(2 * L_sz - i - 1) * i / 2 + i] += ens[i];
      block.at(i * L0_sz + i) = v12[(2 * L_sz - i - 1) * i / 2 + i];
      // #pragma omp for
      for (auto j = i + 1; j < L0_sz; ++j) {
        block.at(i * L0_sz + j) = v12[(2 * L_sz - i - 1) * i / 2 + j];
      }
    }
  }
  v12.clear();

  return 0;
}
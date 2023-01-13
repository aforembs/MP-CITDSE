#include "pes.hpp"

int pes::readConfig(std::string file, std::string &pot, int &L_max, int &l_max,
                    std::vector<int> &state_sz) {
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

int pes::readCt(std::string file, std::vector<std::complex<double>> &ct) {
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

int pes::genPES(std::string pot, int L_max, std::vector<int> &state_sz,
                std::vector<std::complex<double>> &ct, std::string output) {
  std::vector<int> offs;
  auto sum = 0;
  offs.push_back(sum);
  int L_sz;

  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> edata = nullptr;
  std::vector<double> eig;
  int L_full_sz;
  double ion_yield = 0.0;
  double exitation = 0.0;
  // read sum energies
  int off = 0;
  for (auto L = 0; L <= L_max; ++L) {
    L_sz = state_sz[L];
    filename = pot + "CI" + std::to_string(L) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    edata =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En_CI")));
    L_full_sz = edata->getSpace().getSimpleExtentNpoints();
    eig.resize(L_full_sz);
    edata->read(eig.data(), H5::PredType::NATIVE_DOUBLE);
    file->close();

    std::ofstream outfile(output + "_pes" + std::to_string(L) + ".dat",
                          std::ios::out);
    for (auto i = 0; i < L_sz - 1; ++i) {
      outfile << std::setprecision(16) << eig[i] + 2.0 << " "
              << std::norm(
                     ct[off + i]) // * 2.0 / std::abs(eig[i + 1] - eig[i - 1])
              << "\n";
      if (eig[i] + 2.0 > 0.0) {
        ion_yield += std::norm(ct[off + i]);
      }
      if (off + i > 0) {
        exitation += std::norm(ct[off + i]);
      }
    }

    if (L == 0) {
      std::cout << std::setprecision(16)
                << "\nground state pop: " << std::norm(ct[0]) << "\n";
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
  std::cout << std::setprecision(16) << "norm: " << nrm
            << "\nyield: " << ion_yield << "\nexcited population: " << exitation
            << "\n";

  return 0;
}
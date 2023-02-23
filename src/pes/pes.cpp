#include "pes.hpp"

int pes::readConfig(std::string file, std::string &pot, std::string set_base,
                    std::string option, int &L_max,
                    std::vector<int> &state_sz) {
  YAML::Node settings = YAML::LoadFile(file);
  std::cout << "Global Settings:" << std::endl;
  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "  Core Potential:                       " << pot << std::endl;
  L_max = settings[set_base][option].as<int>();
  std::cout << "  max l/L:                              " << L_max << std::endl;

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

  return 0;
}

int pes::readCt(std::string file, std::vector<std::complex<double>> &ct) {
  std::ifstream fl(file);
  std::string temp;
  while (std::getline(fl, temp)) {
    temp = std::regex_replace(temp, std::regex("^ +"), "");

    if (temp[0] == '#') {
      continue;
    }

    std::istringstream iss(temp);
    int idx;
    double real, imag;
    iss >> idx >> real >> imag;
    ct.push_back(std::complex<double>(real, imag));
  }

  return 0;
}

int pes::genPES1e(std::string pot, bool s_flag, int l_max,
                  std::vector<int> &state_sz,
                  std::vector<std::complex<double>> &ct, std::string output) {
  int n, l_sz = state_sz[0];
  hsize_t offset[] = {0}, stride[] = {1}, block[] = {1};
  hsize_t count[] = {static_cast<hsize_t>(l_sz)};
  H5::DataSpace memspace(1, count, NULL);
  // Read the energies of l=0
  auto filename = pot + std::to_string(0) + ".h5";
  auto file =
      std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  file->openAttribute("N").read(H5::PredType::NATIVE_INT32, &n);
  std::vector<double> En(n * (l_max + 1));
  auto E_set =
      std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
  auto espace = E_set->getSpace();
  espace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
  E_set->read(&En[0], H5::PredType::NATIVE_DOUBLE, memspace, espace);
  file->close();

  // offsets for the different l within c(t)=|..l=0..|..l=1..|..l=...|
  std::vector<int> offs;
  auto sum = 0;
  offs.push_back(sum);
  sum += l_sz;
  offs.push_back(sum);

  for (auto l = 1; l <= l_max; ++l) {
    l_sz = state_sz[l];
    count[0] = l_sz;
    memspace.setExtentSimple(1, count, NULL);
    filename = pot + std::to_string(l) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    E_set = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
    espace = E_set->getSpace();
    espace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    E_set->read(&En[offs[l]], H5::PredType::NATIVE_DOUBLE, memspace, espace);
    file->close();
    sum += l_sz;
    offs.push_back(sum);
  }

  std::vector<double> PES_En;
  PES_En.insert(PES_En.end(), En.begin(), En.begin() + state_sz[0]);
  double ion_yield = 0.0;
  // Loop for calculating dP/dEk = |c_i|^2
  std::vector<double> PES(state_sz[l_max]);
  std::vector<double> PES_l(state_sz[l_max]);
  for (auto l = 0; l <= l_max; ++l) {
    for (auto i = 1; i < state_sz[l]; ++i) {
      for (auto k = 1; k < state_sz[0] - 1; ++k) {
        auto ndiff = std::abs(PES_En[k] - PES_En[k - 1]);
        auto pdiff = std::abs(PES_En[k + 1] - PES_En[k]);

        bool cond = (En[offs[l] + i] > PES_En[k] - 0.5 * ndiff) &&
                    (En[offs[l] + i] < PES_En[k] + 0.5 * pdiff);
        if (cond) {
          auto cf_en = std::norm(ct[offs[l] + i]);
          PES_l[k] += cf_en;
          break;
        }
      }
    }
    std::fstream lfile(output + "_pes" + std::to_string(l) + ".dat",
                       std::ios::out);
    for (auto i = 0; i < state_sz[0] - 1; ++i) {
      if (PES_En[i] > 0.0) {
        lfile << PES_En[i] << " " << PES_l[i] << "\n";
        ion_yield += PES_l[i];
      }
    }
    lfile.close();
    if (s_flag) {
      for (auto k = 1; k < state_sz[0] - 1; ++k) {
        PES[k] += PES_l[k];
      }
    }
    std::fill(PES_l.begin(), PES_l.end(), 0);
  }

  if (s_flag) {
    std::fstream outfile(output + "_pes.dat", std::ios::out);
    for (auto i = 0; i < state_sz[0] - 1; ++i) {
      if (PES_En[i] > 0.0) {
        outfile << PES_En[i] << " " << PES[i] << "\n";
      }
    }
    outfile.close();
  }

  // Get the norm of c(t)
  double nrm = 0;
  for (auto &v : ct) {
    nrm += std::norm(v);
  }

  // Print the ground population and the norm of c(t)
  std::cout << std::setprecision(16) << "ground state pop: " << std::norm(ct[0])
            << "\nnorm: " << nrm << "\nyield: " << ion_yield << "\n";

  return 0;
}

int pes::genPES2e(std::string pot, bool s_flag, int L_max,
                  std::vector<int> &state_sz,
                  std::vector<std::complex<double>> &ct, std::string output) {
  int L_sz;
  std::string filename;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::DataSet> edata = nullptr;
  hsize_t count[1] = {1}, offset[1] = {0}, stride[1] = {1}, block[1] = {1};
  std::vector<double> eig, eig_full;
  std::vector<int> offs;
  int sum = 0;
  offs.push_back(sum);
  int L_full_sz;
  double ion_yield = 0.0;
  double exitation = 0.0;
  double threshold = 0.0;

  // read first ionisation threshold
  filename = pot + "0.h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  edata = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
  auto set1 = edata->getSpace();
  set1.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
  edata->read(&threshold, H5::PredType::NATIVE_DOUBLE,
              H5::DataSpace(1, count, NULL), set1);
  file->close();

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

    if (s_flag) {
      sum += L_sz;
      eig_full.insert(eig_full.end(), eig.begin(), eig.begin() + L_sz);
      offs.push_back(sum);
    }

    std::ofstream outfile(output + "_pes" + std::to_string(L) + ".dat",
                          std::ios::out);
    for (auto i = 1; i < L_sz - 1; ++i) {
      outfile << std::setprecision(16) << eig[i] - threshold << " "
              << std::norm(ct[off + i]) << "\n";
      if (eig[i] - threshold > 0.0) {
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

  if (s_flag) {
    std::vector<double> PES(state_sz[L_max]);
    for (auto L = 0; L <= L_max; ++L) {
      for (auto i = 1; i < state_sz[L]; ++i) {
        for (auto k = 1; k < state_sz[L_max]; ++k) {
          auto ndiff = std::abs(eig_full[offs[L_max] + k] -
                                eig_full[offs[L_max] + k - 1]);
          auto pdiff = std::abs(eig_full[offs[L_max] + k + 1] -
                                eig_full[offs[L_max] + k]);

          bool cond =
              (eig_full[offs[L] + i] >
               eig_full[offs[L_max] + k] - 0.5 * ndiff) &&
              (eig_full[offs[L] + i] < eig_full[offs[L_max] + k] + 0.5 * pdiff);
          if (cond) {
            PES[k] += std::norm(ct[offs[L] + i]);
            break;
          }
        }
      }
    }
    std::fstream outfile(output + "_pes.dat", std::ios::out);
    for (auto i = 0; i < state_sz[L_max] - 1; ++i) {
      auto energy = eig_full[offs[L_max] + i] - threshold;
      if (energy > 0.0) {
        outfile << energy << " " << PES[i] << "\n";
      }
    }
    outfile.close();
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
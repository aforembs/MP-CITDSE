#include "genidx.hpp"

int genidx::readConfig(std::string file, std::string &pot, int &L_max) {
  YAML::Node settings = YAML::LoadFile(file);

  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "Core Potential:                              " << pot
            << std::endl;
  L_max = settings["Global_Settings"]["L_max"].as<int>();
  std::cout << "Maximum total two electron angular momentum: " << L_max
            << std::endl;

  return 0;
}

int genidx::sortEn(std::string pot, int L_max, std::string dir) {
  int l1 = 0, last_l1 = 0;
  int l2 = 0, last_l2 = 0;
  std::string filename;
  std::string outfile_name;
  std::unique_ptr<H5::DataSet> e1 = nullptr, ei = nullptr;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::H5File> outfile = nullptr;
  std::vector<double> en12, en;
  std::vector<idx4> idx;

  int ncf, sym, t_sz;
  std::vector<cfg::line> cfgs;
  en_data en_d;
  idx4 idx_elm;

  H5::DataSpace e_space;
  hsize_t offset[1] = {0}, stride[1] = {1}, block[1] = {1};
  hsize_t count[1], dimms[1], write_sz[1], idx_sz[1];
  H5::DataSpace memspace_l;

  // Loop
  for (int L = 0; L <= L_max; ++L) {
    cfg::readCfg(dir, L, sym, ncf, cfgs);

    t_sz = 0;
    for (auto i = 0; i < ncf; ++i) {
      t_sz += cfgs[i].n2max - cfgs[i].n2min;
    }

    auto max_n2l = *std::max_element(
        cfgs.begin(), cfgs.end(), [](cfg::line const &a, cfg::line const &b) {
          return a.n2max < b.n2max;
        });
    auto max_Nsz = max_n2l.n2max;
    en12.reserve(max_Nsz * 2);
    count[0] = max_Nsz;
    dimms[0] = max_Nsz;
    memspace_l.setExtentSimple(1, dimms, NULL);

    en.reserve(t_sz);
    idx.reserve(t_sz);
    write_sz[0] = t_sz;
    idx_sz[0] = t_sz * 4;

    filename = pot + "0.h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    e1 = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
    e_space = e1->getSpace();
    e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);

    e_space = e1->getSpace();
    e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
    file->close();

    last_l1 = 0;
    last_l2 = 0;
    for (const auto &line : cfgs) {
      l1 = line.l1;
      l2 = line.l2;

      if (l1 != last_l1) {
        filename = pot + std::to_string(l1) + ".h5";
        file =
            std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
        e1 =
            std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
        e_space = e1->getSpace();
        e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
        file->close();
      }

      if (l2 != last_l2) {
        filename = pot + std::to_string(l2) + ".h5";
        file =
            std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
        e1 =
            std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
        e_space = e1->getSpace();
        e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l,
                 e_space);
        file->close();
      }

      double ent, en1 = en12[line.n1];
      for (int n2 = line.n2min; n2 < line.n2max; ++n2) {
        ent = en1 + en12[max_Nsz + n2];
        en_d = {ent, line.n1, l1, n2, l2};
        idx_elm = {line.n1, l1, n2, l2};
        en.emplace_back(ent);
        idx.emplace_back(idx_elm);
      }
      last_l1 = l1;
      last_l2 = l2;
    }

    // write the energies to a file
    outfile_name = pot + "2_" + std::to_string(L) + "En.h5";
    outfile =
        std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_TRUNC));
    ei = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
        "e_2e", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, write_sz))));
    ei->write(&en[0], H5::PredType::NATIVE_DOUBLE);

    en.clear();

    // write indices to file
    ei = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
        "idx", H5::PredType::NATIVE_INT32, H5::DataSpace(1, idx_sz))));
    ei->write(&idx[0], H5::PredType::NATIVE_INT32);
    outfile->close();

    idx.clear();
    cfgs.clear();
  }

  return 0;
}
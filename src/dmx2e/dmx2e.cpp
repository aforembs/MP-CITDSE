#include "dmx2e.hpp"

int dmx2e::ReadConfig(std::string file, std::string &pot, int &L_max,
                      int &l_max, char &gauge) {
  YAML::Node settings = YAML::LoadFile(file);

  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "Core Potential:                              " << pot
            << std::endl;
  L_max = settings["Global_Settings"]["L_max"].as<int>();
  std::cout << "Maximum total two electron angular momentum: " << L_max
            << std::endl;
  gauge = settings["Global_Settings"]["gauge"].as<char>();
  std::cout << "Gauge type ('l' length/'v' velocity):        " << gauge
            << std::endl;

  l_max = settings["Basis_Settings"]["l_max"].as<int>();
  std::cout << "Maximum one electron angular momentum:       " << l_max
            << std::endl;
  return 0;
}

/* Function for calculating which l1,l2 combinations belong to which L
 * i.e. for L=0 l1,l2 = 0,0 1,1 2,2 ...
 * @param in L       - The 2e total orbital angular momentum L
 * @param in l1e_max - The maximum 1e angular momentum l1 and l2
 */
en_L make_enL(int L, int l1e_max) {
  int sz = 0;
  int st = 0;
  int l1 = 0;
  int l2 = 0;
  std::vector<l_ab> lp;
  l_ab l_pair;

  if (L == 0) { // For L=0, l1=l2
    sz = l1e_max + 1;
    lp.reserve(sz);
    for (int i = 0; i < sz; ++i) {
      l_pair = {i, i};
      lp.emplace_back(l_pair);
    }
  } else if (L == 1) { // For L=1, l2=l1+1
    sz = l1e_max;
    lp.reserve(sz);
    for (int i = 0; i < sz; ++i) {
      l_pair = {i, i + 1};
      lp.emplace_back(l_pair);
    }
  } else if ((L ^ 1) == (L + 1) && L != 0) { // L>=2, L even
    sz = l1e_max * 2;
    lp.reserve(sz);
    l_pair = {0, L};
    lp.emplace_back(l_pair);
    for (l1 = 1; l1 <= l1e_max; ++l1) {
      if (l1 * 2 >= L) {
        st = l1;
      } else {
        st = L - l1;
      }
      for (l2 = st; l2 <= l1e_max; l2 += 2) {
        l_pair = {l1, l2};
        lp.emplace_back(l_pair);
      }
    }
  } else if (L & 1 && L != 1) { // L>=3, L odd
    sz = l1e_max * 2;
    lp.reserve(sz);
    l_pair = {0, L};
    lp.emplace_back(l_pair);
    for (l1 = 1; l1 <= l1e_max; ++l1) {
      if (l1 * 2 + 1 >= L) {
        st = l1 + 1;
      } else {
        st = L - l1;
      }
      for (l2 = st; l2 <= l1e_max; l2 += 2) {
        l_pair = {l1, l2};
        lp.emplace_back(l_pair);
      }
    }
  }

  auto en_vec = std::vector<en_data>();
  en_L Len = {L, lp, en_vec};
  return Len;
}

int dmx2e::SortL(std::string cpot, int L_max, char gauge, std::string dir) {
  int l1 = 0;
  int l2 = 0;
  std::string filename;
  std::string outfile_name;
  std::unique_ptr<H5::DataSet> e1 = nullptr, ei = nullptr;
  std::unique_ptr<H5::H5File> file = nullptr;
  std::unique_ptr<H5::H5File> outfile = nullptr;
  std::vector<double> en12, en_srtd;
  std::vector<idx4> idx_data;
  std::vector<en_data> Li_dat, Lf_dat;

  int ncf, sym;
  std::vector<cfg::line> cfgs;
  cfg::readCfg(dir, 0, sym, ncf, cfgs);

  en_data en_d;
  idx4 idx_elm;

  H5::DataSpace e_space;
  hsize_t offset[1], count[1], stride[1], block[1], dimms[1];
  offset[0] = 0;
  stride[0] = 1;
  block[0] = 1;
  H5::DataSpace memspace_l;

  hsize_t t_sz = 0;
  for (auto i = 0; i < ncf; ++i) {
    t_sz += cfgs[i].n2max - cfgs[i].n2min;
  }

  auto max_n2l = *std::max_element(
      cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
  auto max_Nsz = max_n2l.n2max;
  en12.reserve(max_Nsz * 2);
  count[0] = max_Nsz;
  dimms[0] = max_Nsz;
  memspace_l.setExtentSimple(1, dimms, NULL);

  // Sort the energies in parallel using C++17 built in parallel sort
  std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(),
            [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });

  std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(),
            [](cfg::line const &a, cfg::line const &b) {
              return (a.l2 < b.l2) && (a.l1 == b.l1);
            });

  // L initial sort & save
  Li_dat.reserve(t_sz);
  en_srtd.reserve(t_sz);
  idx_data.reserve(t_sz);
  hsize_t write_sz[] = {t_sz};
  hsize_t idx_sz[] = {t_sz * 4};

  filename = cpot + "0.h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  e1 = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
  e_space = e1->getSpace();
  e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
  e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);

  e_space = e1->getSpace();
  e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
  e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);

  int last_l1 = 0, last_l2 = 0;
  for (const auto &line : cfgs) {
    l1 = line.l1;
    l2 = line.l2;

    if (l1 != last_l1) {
      filename = cpot + std::to_string(l1) + ".h5";
      file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
    }

    if (l2 != last_l2) {
      filename = cpot + std::to_string(l2) + ".h5";
      file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l,
               e_space);
    }

    double ent, en1 = en12[line.n1];
    for (int n2 = line.n2min; n2 < line.n2max; ++n2) {
      ent = en1 + en12[max_Nsz + n2];
      en_d = {ent, line.n1, l1, n2, l2};
      Li_dat.emplace_back(en_d);
    }
    last_l1 = l1;
    last_l2 = l2;
  }

  // Sort the energies in parallel using C++17 built in parallel sort
  std::sort(std::execution::par_unseq, Li_dat.begin(), Li_dat.end(),
            [](en_data const &a, en_data const &b) { return a.en < b.en; });

  // Separate the sorted energies from the sorted indices
  for (auto i : Li_dat) {
    en_srtd.emplace_back(i.en);
    idx_elm = {i.n1, i.l1, i.n2, i.l2};
    idx_data.emplace_back(idx_elm);
  }

  // write the energies to a file
  outfile_name = cpot + "2_01" + gauge + ".h5";
  outfile =
      std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
      "e_i", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, write_sz))));
  ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);
  outfile->close();

  en_srtd.clear();

  // write indices to file
  outfile_name = cpot + "0idx.h5";
  outfile =
      std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
      "idx", H5::PredType::NATIVE_INT32, H5::DataSpace(1, idx_sz))));
  ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);
  outfile->close();

  idx_data.clear();
  Li_dat.clear();
  cfgs.clear();

  // L final sort & save
  cfg::readCfg(dir, 1, sym, ncf, cfgs);

  t_sz = 0;
  for (auto i = 0; i < ncf; ++i) {
    t_sz += cfgs[i].n2max - cfgs[i].n2min;
  }

  max_n2l = *std::max_element(
      cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
  max_Nsz = max_n2l.n2max;
  en12.reserve(max_Nsz * 2);
  count[0] = max_Nsz;
  dimms[0] = max_Nsz;
  memspace_l.setExtentSimple(1, dimms, NULL);

  // Sort the energies in parallel using C++17 built in parallel sort
  std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(),
            [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });

  std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(),
            [](cfg::line const &a, cfg::line const &b) {
              return (a.l2 < b.l2) && (a.l1 == b.l1);
            });

  Lf_dat.reserve(t_sz);
  en_srtd.reserve(t_sz);
  idx_data.reserve(t_sz);
  write_sz[0] = t_sz;
  idx_sz[0] = t_sz * 4;

  filename = cpot + "0.h5";
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
      filename = cpot + std::to_string(l1) + ".h5";
      file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
      file->close();
    }

    if (l2 != last_l2) {
      filename = cpot + std::to_string(l2) + ".h5";
      file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
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
      Lf_dat.emplace_back(en_d);
    }
    last_l1 = l1;
    last_l2 = l2;
  }

  std::sort(std::execution::par_unseq, Lf_dat.begin(), Lf_dat.end(),
            [](en_data const &a, en_data const &b) { return a.en < b.en; });

  for (auto i : Lf_dat) {
    en_srtd.emplace_back(i.en);
    idx_elm = {i.n1, i.l1, i.n2, i.l2};
    idx_data.emplace_back(idx_elm);
  }

  outfile_name = cpot + "2_01" + gauge + ".h5";
  outfile =
      std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_RDWR));
  ei = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
      "e_f", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, write_sz))));
  ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);
  outfile->close();

  // write indices to file
  outfile_name = cpot + "1idx.h5";
  outfile =
      std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
      "idx", H5::PredType::NATIVE_INT32, H5::DataSpace(1, idx_sz))));
  ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);
  outfile->close();

  Lf_dat.clear();

  // Loop
  for (int L_itr = 2; L_itr <= L_max; ++L_itr) {
    outfile_name = cpot + "2_" + std::to_string(L_itr - 1) +
                   std::to_string(L_itr) + gauge + ".h5";
    outfile =
        std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_TRUNC));
    ei = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
        "e_i", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, write_sz))));
    ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);

    en_srtd.clear();
    idx_data.clear();
    Lf_dat.clear();
    cfgs.clear();

    cfg::readCfg(dir, L_itr, sym, ncf, cfgs);

    t_sz = 0;
    for (auto i = 0; i < ncf; ++i) {
      t_sz += cfgs[i].n2max - cfgs[i].n2min;
    }

    max_n2l = *std::max_element(cfgs.begin(), cfgs.end(),
                                [](cfg::line const &a, cfg::line const &b) {
                                  return a.n2max < b.n2max;
                                });
    max_Nsz = max_n2l.n2max;
    en12.reserve(max_Nsz * 2);
    count[0] = max_Nsz;
    dimms[0] = max_Nsz;
    memspace_l.setExtentSimple(1, dimms, NULL);

    // Sort the energies in parallel using C++17 built in parallel sort
    std::sort(
        std::execution::par_unseq, cfgs.begin(), cfgs.end(),
        [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });

    std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(),
              [](cfg::line const &a, cfg::line const &b) {
                return (a.l2 < b.l2) && (a.l1 == b.l1);
              });

    Lf_dat.reserve(t_sz);
    en_srtd.reserve(t_sz);
    idx_data.reserve(t_sz);
    write_sz[0] = t_sz;
    idx_sz[0] = t_sz * 4;

    filename = cpot + "0.h5";
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
        filename = cpot + std::to_string(l1) + ".h5";
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
        filename = cpot + std::to_string(l2) + ".h5";
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
        Lf_dat.emplace_back(en_d);
      }
      last_l1 = l1;
      last_l2 = l2;
    }

    std::sort(std::execution::par_unseq, Lf_dat.begin(), Lf_dat.end(),
              [](en_data const &a, en_data const &b) { return a.en < b.en; });

    for (auto i : Lf_dat) {
      en_srtd.emplace_back(i.en);
      idx_elm = {i.n1, i.l1, i.n2, i.l2};
      idx_data.emplace_back(idx_elm);
    }

    ei = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
        "e_f", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, write_sz))));
    ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);
    outfile->close();

    // write indices to file
    outfile_name = cpot + std::to_string(L_itr) + "idx.h5";
    outfile =
        std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_TRUNC));
    ei = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
        "idx", H5::PredType::NATIVE_INT32, H5::DataSpace(1, idx_sz))));
    ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);
    outfile->close();
  }

  return 0;
}

int dmx2e::GenDipole(std::string cpot, int L_max, int l_m, char gauge,
                     std::string dir) {
  int l1i = 0, l2i = 0, l1f = 0, l2f = 0;
  int Li_sz = 0;
  int Lf_sz = 0;
  double Lsq = 0.0;
  std::string filename;
  std::string outfile_name;
  std::unique_ptr<H5::H5File> file = nullptr, outfile = nullptr;
  std::unique_ptr<H5::DataSet> dmx = nullptr, L_set = nullptr;
  H5::DataSpace dmx_space;
  hsize_t offset[2], count[2], stride[2], block[2];
  hsize_t dimms[2], T_dimms[2];
  offset[0] = 0;
  offset[1] = 0;
  stride[0] = 1;
  stride[1] = 1;
  block[0] = 1;
  block[1] = 1;
  H5::DataSpace memspace;

  std::vector<int> dmx_sz(L_max);
  std::vector<double> D_data;
  std::vector<double> T;

  int ncf, sym;
  std::vector<cfg::line> cfgs;
  std::vector<idx4 *> buffs(2);
  auto max_N = 0, max_nl = 0;
  cfg::line max_n2l;

  for (auto Li = 0; Li <= L_max; ++Li) {
    cfg::readCfg(dir, Li, sym, ncf, cfgs);

    max_n2l = *std::max_element(cfgs.begin(), cfgs.end(),
                                [](cfg::line const &a, cfg::line const &b) {
                                  return a.n2max < b.n2max;
                                });
    max_N = std::max(max_n2l.n2max, max_N);
    max_nl = std::max(ncf, max_nl);
  }
  auto max_Nsz = max_nl * max_N;
  std::vector<idx4> Li_idx(max_Nsz);
  std::vector<idx4> Lf_idx(max_Nsz);

  count[0] = max_N;
  count[1] = max_N;
  dimms[0] = max_N;
  dimms[1] = max_N;
  D_data.reserve(max_N * max_N * l_m);
  auto max2 = max_N * max_N;

  // Read all 1e dipoles
  for (auto i = 0; i < l_m; ++i) {
    filename = cpot + std::to_string(i) + std::to_string(i + 1) + gauge + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    dmx = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("d_if")));
    dmx_space = dmx->getSpace();
    dmx_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    memspace.setExtentSimple(2, dimms, NULL);
    dmx->read(&D_data[i * max2], H5::PredType::NATIVE_DOUBLE, memspace,
              dmx_space);
    file->close();
  }

  // for loop over L's
  int L = 0;
  int Lf_i = 1;

  // read L indices
  filename = cpot + std::to_string(L) + "idx.h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  L_set = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("idx")));
  Li_sz = L_set->getSpace().getSimpleExtentNpoints() / 4; // get the size of L
  L_set->read(&Li_idx[0], H5::PredType::NATIVE_INT32);
  file->close();

  filename = cpot + std::to_string(Lf_i) + "idx.h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  L_set = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("idx")));
  Lf_sz = L_set->getSpace().getSimpleExtentNpoints() / 4; // get the size of L+1
  L_set->read(&Lf_idx[0], H5::PredType::NATIVE_INT32);
  file->close();

  T.reserve(Li_sz * Lf_sz);
  T_dimms[0] = Lf_sz;
  T_dimms[1] = Li_sz;

  wig_table_init(2 * l_m, 6);

  // calculate 2e dipoles
  Lsq = sqrt(Lf_i);
#pragma omp parallel
  {
    wig_thread_temp_init(2 * l_m);

    for (int idLf = 0; idLf < Lf_sz; ++idLf) { // up to sz Lf
#pragma omp for private(l1i, l2i, l1f, l2f)
      for (int idLi = 0; idLi < Li_sz; ++idLi) { // up to sz Li
        l1i = Li_idx[idLi].l1;
        l2i = Li_idx[idLi].l2;
        l1f = Lf_idx[idLf].l1;
        l2f = Lf_idx[idLf].l2;
        auto Tif = 0.0;
        if (l1i + 1 == l1f && l2i == l2f) {
          Tif = Lsq * pow(-1, l2i + l1f) * sqrt(4.0 * l1f * l1f - 1.0) *
                wig6jj(0, 2, 2, 2 * l1f, 2 * l1i,
                       2 * l2i) * // wigner_6j_2e(L,l1f,l1i,l2i)*
                D_data[l1i * max2 + Li_idx[idLi].n1 + max_N * Lf_idx[idLf].n1] /
                sqrt(static_cast<double>(l1f));
        } else if (l2i + 1 == l2f && l1i == l1f) {
          Tif = Lsq * pow(-1, l1i + l2f) * sqrt(4.0 * l2f * l2f - 1.0) *
                wig6jj(0, 2, 2, 2 * l2f, 2 * l2i,
                       2 * l1i) * // wigner_6j_2e(L,l2f,l2i,l1i)*
                D_data[l2i * max2 + Li_idx[idLi].n2 + max_N * Lf_idx[idLf].n2] /
                sqrt(static_cast<double>(l2f));
        }
        if (l1i == l2i) {
          Tif *= 0.7071067811865475244008444e0;
        }
        if (l1f == l2f) {
          Tif *= 0.7071067811865475244008444e0;
        }
        T[idLi + idLf * Li_sz] = Tif;
      }
    }

    wig_temp_free();
  }
  // Save TL1;L2
  outfile_name =
      cpot + "2_" + std::to_string(L) + std::to_string(Lf_i) + gauge + ".h5";
  outfile =
      std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_RDWR));
  dmx = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
      "d_if", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, T_dimms))));
  dmx->write(&T[0], H5::PredType::NATIVE_DOUBLE);
  outfile->close();

  T.clear();

  // Set L=1, L+1=2 and change index buffer pointers
  int L2 = 0;
  int buf_Li = 1;
  int buf_Lf = 0;
  Li_sz = Lf_sz;
  buffs[0] = &Li_idx[0];
  buffs[1] = &Lf_idx[0];

  // Loop for L>0
  for (Lf_i = 2; Lf_i <= L_max; ++Lf_i) {
    L = Lf_i - 1;
    L2 = 2 * L;

    // read L indices
    filename = cpot + std::to_string(Lf_i) + "idx.h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    L_set =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("idx")));
    Lf_sz = L_set->getSpace().getSimpleExtentNpoints() / 4;
    L_set->read(buffs.at(buf_Lf), H5::PredType::NATIVE_INT32);
    file->close();

    T.reserve(Li_sz * Lf_sz);
    T_dimms[0] = Lf_sz;
    T_dimms[1] = Li_sz;

    // calculate 2e dipoles
    Lsq = sqrt(Lf_i);
#pragma omp parallel
    {
      wig_thread_temp_init(2 * l_m);
      for (int idLf = 0; idLf < Lf_sz; ++idLf) { // up to sz Lf
#pragma omp for private(l1i, l2i, l1f, l2f)
        for (int idLi = 0; idLi < Li_sz; ++idLi) { // up to sz Li
          l1i = buffs.at(buf_Li)[idLi].l1;
          l2i = buffs.at(buf_Li)[idLi].l2;
          l1f = buffs.at(buf_Lf)[idLf].l1;
          l2f = buffs.at(buf_Lf)[idLf].l2;
          auto Tif = 0.0;
          if (l1i + 1 == l1f && l2i == l2f) {
            Tif = Lsq * pow(-1, l2i + l1f) * sqrt(4.0 * l1f * l1f - 1.0) *
                  wig6jj(L2, L2 + 2, 2, 2 * l1f, 2 * l1i,
                         2 * l2i) * // wigner_6j_2e(L,l1f,l1i,l2i)*
                  D_data[l1i * max2 + buffs.at(buf_Li)[idLi].n1 +
                         max_N * buffs.at(buf_Lf)[idLf].n1] /
                  sqrt(static_cast<double>(l1f));
          } else if (l2i + 1 == l2f && l1i == l1f) {
            Tif = Lsq * pow(-1, l1i + l2f) * sqrt(4.0 * l2f * l2f - 1.0) *
                  wig6jj(L2, L2 + 2, 2, 2 * l2f, 2 * l2i,
                         2 * l1i) * // wigner_6j_2e(L,l2f,l2i,l1i)*
                  D_data[l2i * max2 + buffs.at(buf_Li)[idLi].n2 +
                         max_N * buffs.at(buf_Lf)[idLf].n2] /
                  sqrt(static_cast<double>(l2f));
          }
          if (l1i == l2i) {
            Tif *= 0.7071067811865475244008444e0;
          }
          if (l1f == l2f) {
            Tif *= 0.7071067811865475244008444e0;
          }
          T[idLi + idLf * Li_sz] = Tif;
        }
      }
      wig_temp_free();
    }

    outfile_name =
        cpot + "2_" + std::to_string(L) + std::to_string(Lf_i) + gauge + ".h5";
    outfile =
        std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_RDWR));
    dmx = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
        "d_if", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, T_dimms))));
    dmx->write(&T[0], H5::PredType::NATIVE_DOUBLE);
    outfile->close();

    T.clear();

    buf_Li = buf_Lf;
    buf_Lf = 1 - buf_Lf;
    Li_sz = Lf_sz;
  }

  wig_table_free();
  return 0;
}
#include "dmx2e.hpp"

int dmx2e::readConfig(std::string file, std::string &pot, int &L_max,
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

int dmx2e::genDipole(std::string pot, int L_max, int l_max, char gauge,
                     std::string dir) {
  constexpr double sq2 = 0.7071067811865475244008444e0;
  int l1i = 0, l2i = 0, l1f = 0, l2f = 0;
  int n1i = 0, n2i = 0, n1f = 0, n2f = 0;
  int Li_sz = 0;
  int Lf_sz = 0;
  double Lsq = 0.0;
  double dipT = 0.0;
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
  std::vector<dmtp::idx4 *> buffs(2);
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
  std::vector<dmtp::idx4> Li_idx(max_Nsz);
  std::vector<dmtp::idx4> Lf_idx(max_Nsz);

  count[0] = max_N;
  count[1] = max_N;
  dimms[0] = max_N;
  dimms[1] = max_N;
  D_data.resize(max_N * max_N * l_max);
  auto max2 = max_N * max_N;

  // Read all 1e dipoles
  for (auto i = 0; i < l_max; ++i) {
    filename = pot + std::to_string(i) + std::to_string(i + 1) + gauge + ".h5";
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
  filename = pot + "2_" + std::to_string(L) + "En.h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  L_set = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("idx")));
  Li_sz = L_set->getSpace().getSimpleExtentNpoints() / 4; // get the size of L
  L_set->read(&Li_idx[0], H5::PredType::NATIVE_INT32);
  file->close();

  filename = pot + "2_" + std::to_string(Lf_i) + "En.h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  L_set = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("idx")));
  Lf_sz = L_set->getSpace().getSimpleExtentNpoints() / 4; // get the size of L+1
  L_set->read(&Lf_idx[0], H5::PredType::NATIVE_INT32);
  file->close();

  T.resize(Li_sz * Lf_sz);
  T_dimms[0] = Li_sz;
  T_dimms[1] = Lf_sz;

  switch (gauge) {
  case 'v':
    dipT = -1.0;
    break;
  case 'l':
    dipT = 1.0;
    break;
  }

  wig_table_init(2 * l_max, 6);

  // calculate 2e dipoles
  Lsq = pow(-1, Lf_i) * sqrt(static_cast<double>(Lf_i));
#pragma omp parallel
  {
    wig_thread_temp_init(2 * l_max);

    for (int idLi = 0; idLi < Li_sz; ++idLi) { // up to sz Li
#pragma omp for private(l1i, l2i, l1f, l2f, n1i, n2i, n1f, n2f)
      for (int idLf = 0; idLf < Lf_sz; ++idLf) { // up to sz Lf
        l1i = Li_idx[idLi].l1;
        n1i = Li_idx[idLi].n1;
        l2i = Li_idx[idLi].l2;
        n2i = Li_idx[idLi].n2;
        l1f = Lf_idx[idLf].l1;
        n1f = Lf_idx[idLf].n1;
        l2f = Lf_idx[idLf].l2;
        n2f = Lf_idx[idLf].n2;
        auto Tif = 0.0;
        if (l1i + 1 == l1f && l2i == l2f && n2i == n2f) {
          Tif += pow(-1, l2i + l1f) * sqrt(4.0 * l1f * l1f - 1.0) *
                 wig6jj(2 * l1i, 2 * l2i, 0, 2, 2, 2 * l1f) *
                 D_data[l1i * max2 + n1f + max_N * n1i] /
                 sqrt(static_cast<double>(l1f));
        }
        if (l2i + 1 == l2f && l1i == l1f && n1i == n1f) {
          Tif += pow(-1, l1i + l2f) * sqrt(4.0 * l2f * l2f - 1.0) *
                 wig6jj(2 * l2i, 2 * l1i, 0, 2, 2, 2 * l2f) *
                 D_data[l2i * max2 + n2f + max_N * n2i] /
                 sqrt(static_cast<double>(l2f));
        }
        if (l1i + 1 == l2f && l2i == l1f && n2i == n1f) {
          Tif += pow(-1, l2i + l2f) * sqrt(4.0 * l2f * l2f - 1.0) *
                 wig6jj(2 * l1i, 2 * l2i, 0, 2, 2, 2 * l2f) *
                 D_data[l1i * max2 + n2f + max_N * n1i] /
                 sqrt(static_cast<double>(l2f));
        }
        if (l2i + 1 == l1f && l1i == l2f && n1i == n2f) {
          Tif += pow(-1, l1i + l1f) * sqrt(4.0 * l1f * l1f - 1.0) *
                 wig6jj(2 * l2i, 2 * l1i, 0, 2, 2, 2 * l1f) *
                 D_data[l2i * max2 + n1f + max_N * n2i] /
                 sqrt(static_cast<double>(l1f));
        }
        if (l1i - 1 == l1f && l2i == l2f && n2i == n2f) {
          Tif += dipT * pow(-1, l2i + l1i) * sqrt(4.0 * l1i * l1i - 1.0) *
                 wig6jj(2 * l1i, 2 * l2i, 0, 2, 2, 2 * l1f) *
                 D_data[l1f * max2 + n1i + max_N * n1f] /
                 sqrt(static_cast<double>(l1i));
        }
        if (l2i - 1 == l2f && l1i == l1f && n1i == n1f) {
          Tif += dipT * pow(-1, l1i + l2i) * sqrt(4.0 * l2i * l2i - 1.0) *
                 wig6jj(2 * l2i, 2 * l1i, 0, 2, 2, 2 * l2f) *
                 D_data[l2f * max2 + n2i + max_N * n2f] /
                 sqrt(static_cast<double>(l2i));
        }
        if (l1i - 1 == l2f && l2i == l1f && n2i == n1f) {
          Tif += dipT * pow(-1, l2i + l1i) * sqrt(4.0 * l1i * l1i - 1.0) *
                 wig6jj(2 * l1i, 2 * l2i, 0, 2, 2, 2 * l2f) *
                 D_data[l2f * max2 + n1i + max_N * n2f] /
                 sqrt(static_cast<double>(l1i));
        }
        if (l2i - 1 == l1f && l1i == l2f && n1i == n2f) {
          Tif += dipT * pow(-1, l1i + l2i) * sqrt(4.0 * l2i * l2i - 1.0) *
                 wig6jj(2 * l2i, 2 * l1i, 0, 2, 2, 2 * l1f) *
                 D_data[l1f * max2 + n2i + max_N * n1f] /
                 sqrt(static_cast<double>(l2i));
        }
        if (l1i == l2i && n1i == n2i) {
          Tif *= sq2;
        }
        if (l1f == l2f && n1f == n2f) {
          Tif *= sq2;
        }
        T[idLf + idLi * Lf_sz] = Tif * Lsq;
      }
    }

    wig_temp_free();
  }
  // Save TL1;L2
  outfile_name =
      pot + "2_" + std::to_string(L) + std::to_string(Lf_i) + gauge + ".h5";
  outfile =
      std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_TRUNC));
  dmx = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
      "d_if", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, T_dimms))));
  dmx->write(T.data(), H5::PredType::NATIVE_DOUBLE);
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
    filename = pot + "2_" + std::to_string(Lf_i) + "En.h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    L_set =
        std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("idx")));
    Lf_sz = L_set->getSpace().getSimpleExtentNpoints() / 4;
    L_set->read(buffs.at(buf_Lf), H5::PredType::NATIVE_INT32);
    file->close();

    T.resize(Li_sz * Lf_sz);
    T_dimms[0] = Li_sz;
    T_dimms[1] = Lf_sz;

    // calculate 2e dipoles
    Lsq = pow(-1, Lf_i) * sqrt(static_cast<double>(Lf_i));
#pragma omp parallel
    {
      wig_thread_temp_init(2 * l_max);
      for (int idLi = 0; idLi < Li_sz; ++idLi) { // up to sz Li
#pragma omp for private(l1i, l2i, l1f, l2f, n1i, n2i, n1f, n2f)
        for (int idLf = 0; idLf < Lf_sz; ++idLf) { // up to sz Lf
          l1i = buffs.at(buf_Li)[idLi].l1;
          n1i = buffs.at(buf_Li)[idLi].n1;
          l2i = buffs.at(buf_Li)[idLi].l2;
          n2i = buffs.at(buf_Li)[idLi].n2;
          l1f = buffs.at(buf_Lf)[idLf].l1;
          n1f = buffs.at(buf_Lf)[idLf].n1;
          l2f = buffs.at(buf_Lf)[idLf].l2;
          n2f = buffs.at(buf_Lf)[idLf].n2;
          auto Tif = 0.0;
          if (l1i + 1 == l1f && l2i == l2f && n2i == n2f) {
            Tif += pow(-1, l2i + l1f) * sqrt(4.0 * l1f * l1f - 1.0) *
                   wig6jj(2 * l1i, 2 * l2i, L2, L2 + 2, 2, 2 * l1f) *
                   D_data[l1i * max2 + n1f + max_N * n1i] /
                   sqrt(static_cast<double>(l1f));
          }
          if (l2i + 1 == l2f && l1i == l1f && n1i == n1f) {
            Tif += pow(-1, l1i + l2f) * sqrt(4.0 * l2f * l2f - 1.0) *
                   wig6jj(2 * l2i, 2 * l1i, L2, L2 + 2, 2, 2 * l2f) *
                   D_data[l2i * max2 + n2f + max_N * n2i] /
                   sqrt(static_cast<double>(l2f));
          }
          if (l1i + 1 == l2f && l2i == l1f && n2i == n1f) {
            Tif += pow(-1, l2i + l2f) * sqrt(4.0 * l2f * l2f - 1.0) *
                   wig6jj(2 * l1i, 2 * l2i, L2, L2 + 2, 2, 2 * l2f) *
                   D_data[l1i * max2 + n2f + max_N * n1i] /
                   sqrt(static_cast<double>(l2f));
          }
          if (l2i + 1 == l1f && l1i == l2f && n1i == n2f) {
            Tif += pow(-1, l1i + l1f) * sqrt(4.0 * l1f * l1f - 1.0) *
                   wig6jj(2 * l2i, 2 * l1i, L2, L2 + 2, 2, 2 * l1f) *
                   D_data[l2i * max2 + n1f + max_N * n2i] /
                   sqrt(static_cast<double>(l1f));
          }
          if (l1i - 1 == l1f && l2i == l2f && n2i == n2f) {
            Tif += dipT * pow(-1, l2i + l1i) * sqrt(4.0 * l1i * l1i - 1.0) *
                   wig6jj(2 * l1i, 2 * l2i, L2, L2 + 2, 2, 2 * l1f) *
                   D_data[l1f * max2 + n1i + max_N * n1f] /
                   sqrt(static_cast<double>(l1i));
          }
          if (l2i - 1 == l2f && l1i == l1f && n1i == n1f) {
            Tif += dipT * pow(-1, l1i + l2i) * sqrt(4.0 * l2i * l2i - 1.0) *
                   wig6jj(2 * l2i, 2 * l1i, L2, L2 + 2, 2, 2 * l2f) *
                   D_data[l2f * max2 + n2i + max_N * n2f] /
                   sqrt(static_cast<double>(l2i));
          }
          if (l1i - 1 == l2f && l2i == l1f && n2i == n1f) {
            Tif += dipT * pow(-1, l2i + l1i) * sqrt(4.0 * l1i * l1i - 1.0) *
                   wig6jj(2 * l1i, 2 * l2i, L2, L2 + 2, 2, 2 * l2f) *
                   D_data[l2f * max2 + n1i + max_N * n2f] /
                   sqrt(static_cast<double>(l1i));
          }
          if (l2i - 1 == l1f && l1i == l2f && n1i == n2f) {
            Tif += dipT * pow(-1, l1i + l2i) * sqrt(4.0 * l2i * l2i - 1.0) *
                   wig6jj(2 * l2i, 2 * l1i, L2, L2 + 2, 2, 2 * l1f) *
                   D_data[l1f * max2 + n2i + max_N * n1f] /
                   sqrt(static_cast<double>(l2i));
          }
          if (l1i == l2i && n1i == n2i) {
            Tif *= sq2;
          }
          if (l1f == l2f && n1f == n2f) {
            Tif *= sq2;
          }
          T[idLf + idLi * Lf_sz] = Tif * Lsq;
        }
      }
      wig_temp_free();
    }

    outfile_name =
        pot + "2_" + std::to_string(L) + std::to_string(Lf_i) + gauge + ".h5";
    outfile =
        std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_TRUNC));
    dmx = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
        "d_if", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, T_dimms))));
    dmx->write(T.data(), H5::PredType::NATIVE_DOUBLE);
    outfile->close();

    T.clear();

    buf_Li = buf_Lf;
    buf_Lf = 1 - buf_Lf;
    Li_sz = Lf_sz;
  }

  wig_table_free();
  return 0;
}
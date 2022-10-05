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
  filename = cpot + "2_" + std::to_string(L) + "En.h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  L_set = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("idx")));
  Li_sz = L_set->getSpace().getSimpleExtentNpoints() / 4; // get the size of L
  L_set->read(&Li_idx[0], H5::PredType::NATIVE_INT32);
  file->close();

  filename = cpot + "2_" + std::to_string(Lf_i) + "En.h5";
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
      std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_TRUNC));
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
    filename = cpot + "2_" + std::to_string(Lf_i) + "En.h5";
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
        std::make_unique<H5::H5File>(H5::H5File(outfile_name, H5F_ACC_TRUNC));
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
#include "r12.hpp"
#include "time_tst.hpp"

int r_12::readConfig(std::string file, int &qsz, std::string &pot, int &L_max,
                     char &gauge) {
  YAML::Node settings = YAML::LoadFile(file);

  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "Core Potential:                              " << pot
            << std::endl;
  qsz = settings["Global_Settings"]["Outer_quadrature_size"].as<int>();
  std::cout << "No. of GL-quadrature points between knots:   " << qsz
            << std::endl;
  L_max = settings["Global_Settings"]["L_max"].as<int>();
  std::cout << "Maximum total two electron angular momentum: " << L_max
            << std::endl;
  gauge = settings["Global_Settings"]["gauge"].as<char>();
  std::cout << "Gauge type ('l' length/'v' velocity):        " << gauge
            << std::endl;

  return 0;
}

int Rpowk(int qsz, int pti_sz, int k_max, std::vector<double> &qx_o,
          std::vector<double> &qx_i, std::vector<double> &r_out,
          std::vector<double> &r_in) {
  r_out.resize((k_max + 1) * qsz);
  r_in.resize((k_max + 1) * pti_sz);
  auto km1 = 0;

  std::fill(r_out.begin(), r_out.begin() + qsz, 1.0);
  std::fill(r_in.begin(), r_in.begin() + pti_sz, 1.0);

  std::copy(qx_o.begin(), qx_o.end(), r_out.begin() + qsz);
  std::copy(qx_i.begin(), qx_i.end(), r_in.begin() + pti_sz);

  for (auto k = 2; k <= k_max; ++k) {
    km1 = k - 1;
    for (int i = 0; i < qsz; ++i) {
      r_out[i + k * qsz] = r_out[i + qsz] * r_out[i + km1 * qsz];
    }
    for (int i = 0; i < pti_sz; ++i) {
      r_in[i + k * pti_sz] = r_in[i + pti_sz] * r_in[i + km1 * pti_sz];
    }
  }

  return 0;
}

int r_12::r12Glob(std::string cpot, int L_max, int qsz, std::string dir) {
  bool min_dir = 0, min_exc = 0;
  double Y_norm = 0.0;
  std::vector<double> v_mat;
  double sum_k = 0.0;
  idx4 e12, e12p;

  int L_sz = 0, v_sz = 0;
  std::string filename;
  std::string outfile_name;

  // HDF5 defines
  std::unique_ptr<H5::H5File> outfile = nullptr, file = nullptr;
  std::unique_ptr<H5::DataSet> Po = nullptr, Pi = nullptr, Pidx = nullptr,
                               qr = nullptr, qw = nullptr, qri = nullptr;
  std::unique_ptr<H5::DataSet> V_set = nullptr;
  std::unique_ptr<H5::DataSet> L_set = nullptr;
  std::vector<idx4> L_idx;
  H5::DataSpace cspace;
  hsize_t offset[2] = {0, 0}, stride[2] = {1, 1}, block[2] = {1, 1};
  hsize_t count[2], counti[2], dimms[2], dimmsi[2], v_dim[1];
  H5::DataSpace memspace, memspacei;

  int ncf, sym, max_N = 0, l1_m = 0, l2_m = 0;
  std::vector<cfg::line> cfgs;
  cfg::line max_line;

  for (auto li = 0; li <= L_max; ++li) {
    cfg::readCfg(dir, li, sym, ncf, cfgs);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
                                 [](cfg::line const &a, cfg::line const &b) {
                                   return a.n2max < b.n2max;
                                 });
    max_N = std::max(max_N, max_line.n2max);

    max_line = *std::max_element(
        cfgs.begin(), cfgs.end(),
        [](cfg::line const &a, cfg::line const &b) { return a.l2 < b.l2; });
    l2_m = std::max(l2_m, max_line.l2);

    max_line = *std::max_element(
        cfgs.begin(), cfgs.end(),
        [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });
    l1_m = std::max(l1_m, max_line.l1);
  }
  auto k_max = l1_m + l2_m;
  auto lc_sz = max_N * qsz;

  count[0] = max_N;
  counti[0] = max_N;
  count[1] = qsz;
  dimms[0] = count[0];
  dimmsi[0] = counti[0];
  dimms[1] = count[1];

  memspace.setExtentSimple(2, dimms, NULL);

  // reserve space for wavefunctions
  auto e1_lm = std::max(l1_m, l2_m);
  std::vector<double> wfn_o(lc_sz * (e1_lm + 1));

  // read data from l=0
  filename = cpot + "_w1e0.h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  Po = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Pr_o")));
  cspace = Po->getSpace();
  int nen = cspace.getSimpleExtentNpoints() / qsz;
  cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
  Po->read(wfn_o.data(), H5::PredType::NATIVE_DOUBLE, memspace, cspace);
  Pi = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Pr_i")));
  cspace = Pi->getSpace();

  int pti_sz = cspace.getSimpleExtentNpoints() / nen;
  auto lci_sz = max_N * pti_sz;
  std::vector<double> wfn_i(lci_sz * (e1_lm + 1));
  std::vector<uint8_t> pq_dx(qsz);
  std::vector<double> qx_o(qsz), qw_o(qsz), qx_i(pti_sz);
  counti[1] = pti_sz;
  dimmsi[1] = counti[1];
  memspacei.setExtentSimple(2, dimmsi, NULL);

  cspace.selectHyperslab(H5S_SELECT_SET, counti, offset, stride, block);
  Pi->read(wfn_i.data(), H5::PredType::NATIVE_DOUBLE, memspacei, cspace);
  Pidx =
      std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Pr_idx")));
  Pidx->read(pq_dx.data(), H5::PredType::NATIVE_UCHAR);
  qr = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Qr_o")));
  qr->read(qx_o.data(), H5::PredType::NATIVE_DOUBLE);
  qw = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Qw_o")));
  qw->read(qw_o.data(), H5::PredType::NATIVE_DOUBLE);
  qri = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Qr_i")));
  qri->read(qx_i.data(), H5::PredType::NATIVE_DOUBLE);
  file->close();

  // read wavefunctions for all other l
  for (int l = 1; l <= e1_lm; ++l) {
    filename = cpot + "_w1e" + std::to_string(l) + ".h5";
    file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
    Po = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Pr_o")));
    cspace = Po->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    Po->read(&wfn_o[l * lc_sz], H5::PredType::NATIVE_DOUBLE, memspace, cspace);

    Pi = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Pr_i")));
    cspace = Pi->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, counti, offset, stride, block);
    Pi->read(&wfn_i[l * lci_sz], H5::PredType::NATIVE_DOUBLE, memspacei,
             cspace);
    file->close();
  }

  // Precalculate r^k
  std::vector<double> rk, rk_in;
  Rpowk(qsz, pti_sz, k_max, qx_o, qx_i, rk, rk_in);

  omp_lock_t copylock;
  omp_init_lock(&copylock);
  uint64_t st_time;

  wig_table_init(2 * (k_max + 2), 6);

#pragma omp parallel private(e12p)
  {
    wig_thread_temp_init(2 * (k_max + 2));

    for (int L = 0; L <= L_max; ++L) {

// Read n1l1;n2l2 indices for NL states
#pragma omp single
      {
        filename = cpot + "2_" + std::to_string(L) + "En.h5";
        file =
            std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
        L_set = std::make_unique<H5::DataSet>(
            H5::DataSet(file->openDataSet("idx")));
        L_sz = L_set->getSpace().getSimpleExtentNpoints() / 4;
        L_idx.reserve(L_sz);
        L_set->read(&L_idx[0], H5::PredType::NATIVE_INT32);
        file->close();

        std::cout << "L: " << L << " L_sz: " << L_sz << "\n";

        v_sz = L_sz * (L_sz + 1) / 2;
        v_mat.reserve(v_sz); // error

        st_time = GetTimeMs64();
      }

      for (int NL2 = 0; NL2 < L_sz; ++NL2) {
        // set n1'l1';n2'l2'
        e12p = L_idx[NL2];
        bool eqvp = e12p.n1 == e12p.n2 && e12p.l1 == e12p.l2;

#pragma omp barrier
#pragma omp for private(e12, Y_norm, sum_k, min_dir, min_exc)
        for (int NL1 = NL2; NL1 < L_sz; ++NL1) {
          // set n1l1;n2l2
          e12 = L_idx[NL1];

          // sqrt([la][lc][lb][ld])
          Y_norm = sqrt((2 * e12.l1 + 1) * (2 * e12p.l1 + 1) *
                        (2 * e12.l2 + 1) * (2 * e12p.l2 + 1));
          sum_k = 0.0;
          min_dir =
              ((L + e12.l2 + e12p.l1) >> 0) & 1; // check if L+lb+lc is odd
          min_exc = ((L + e12.l1 + e12p.l1) >> 0) & 1;
          for (int k = 0; k <= k_max; ++k) {
            /* for direct check if:
              (-)^{L+lb+lc}=(-)^{L+la+ld},
              |la-lc| <= k <= la+lc,
              |lb-ld| <= k <= lb+ld */
            if (min_dir == (((L + e12.l1 + e12p.l2) >> 0) & 1) &&
                ((abs(e12.l1 - e12p.l1) <= k) && (k <= e12.l1 + e12p.l1)) &&
                ((abs(e12.l2 - e12p.l2) <= k) && (k <= e12.l2 + e12p.l2)) &&
                !(((k + e12.l1 + e12p.l1) >> 0) & 1) &&
                !(((k + e12.l2 + e12p.l2) >> 0) & 1) &&
                !std::isinf(1.0 / rk_in[(k + 1) * pti_sz])) {
              sum_k += pow(-1, min_dir) *
                       intfn::fsltrLob(k, qsz, pti_sz, lc_sz, lci_sz, e12.n1,
                                       e12.l1, e12.n2, e12.l2, e12p.n1, e12p.l1,
                                       e12p.n2, e12p.l2, qw_o, pq_dx, rk, rk_in,
                                       wfn_o, wfn_i) *
                       wig3jj(2 * e12.l1, 2 * k, 2 * e12p.l1, 0, 0, 0) *
                       wig3jj(2 * e12.l2, 2 * k, 2 * e12p.l2, 0, 0, 0) *
                       wig6jj(2 * e12p.l1, 2 * e12p.l2, 2 * L, 2 * e12.l2,
                              2 * e12.l1, 2 * k);
            }
            /* for exchange check if:
              (-)^{L+la+lc}=(-)^{L+lb+ld},
              |la-ld| <= k <= la+ld,
              |lc-lb| <= k <= lc+lb */
            if (min_exc == (((L + e12.l2 + e12p.l2) >> 0) & 1) &&
                ((abs(e12.l1 - e12p.l2) <= k) && (k <= e12.l1 + e12p.l2)) &&
                ((abs(e12.l2 - e12p.l1) <= k) && (k <= e12.l2 + e12p.l1)) &&
                !(((k + e12.l1 + e12p.l2) >> 0) & 1) &&
                !(((k + e12.l2 + e12p.l1) >> 0) & 1) &&
                !std::isinf(1.0 / rk_in[(k + 1) * pti_sz])) {
              sum_k += pow(-1, min_exc) *
                       intfn::fsltrLob(k, qsz, pti_sz, lc_sz, lci_sz, e12.n1,
                                       e12.l1, e12.n2, e12.l2, e12p.n2, e12p.l2,
                                       e12p.n1, e12p.l1, qw_o, pq_dx, rk, rk_in,
                                       wfn_o, wfn_i) *
                       wig3jj(2 * e12.l1, 2 * k, 2 * e12p.l2, 0, 0, 0) *
                       wig3jj(2 * e12.l2, 2 * k, 2 * e12p.l1, 0, 0, 0) *
                       wig6jj(2 * e12p.l1, 2 * e12p.l2, 2 * L, 2 * e12.l1,
                              2 * e12.l2, 2 * k);
            }
          }

          // write symmetric V_12 as upper triangular
          double v12 = pow(-1, e12.l1 + e12.l2) * Y_norm * sum_k;

          if (e12.n1 == e12.n2 && e12.l1 == e12.l2)
            v12 *= 0.7071067811865475244008444e0;

          if (eqvp)
            v12 *= 0.7071067811865475244008444e0;

          v_mat[(2 * L_sz - NL2 - 1) * NL2 / 2 + NL1] = v12;
        }
      }

#pragma omp single
      {
        std::cout << "loop time: "
                  << ((double)(GetTimeMs64() - st_time) / 1000.0) << "s\n";
        // save upper triangular V_12
        v_dim[0] = v_sz;
        outfile_name = cpot + "V12_" + std::to_string(L) + ".h5";
        outfile = std::make_unique<H5::H5File>(
            H5::H5File(outfile_name, H5F_ACC_TRUNC));
        V_set =
            std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
                "V_12", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, v_dim))));
        V_set->write(&v_mat[0], H5::PredType::NATIVE_DOUBLE);
        outfile->close();
        v_mat.clear();
      }
    }

    wig_temp_free();
  }

  wig_table_free();
  omp_destroy_lock(&copylock);

  return 0;
}
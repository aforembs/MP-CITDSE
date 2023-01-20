#include "w1e.hpp"

namespace glx {
const double GLobx3[1] = {0.0};
const double GLobx4[2] = {-0.447213595499957939281834733746,
                          0.447213595499957939281834733746};
const double GLobx5[3] = {-0.65465367070797714379829245625, 0.0,
                          0.65465367070797714379829245625};
const double GLobx6[4] = {
    -0.765055323929464692851002973959, -0.285231516480645096314150994041,
    0.285231516480645096314150994041, 0.765055323929464692851002973959};
const double GLobx7[5] = {
    -0.830223896278566929872032213968, -0.468848793470714213803771881909, 0.0,
    0.468848793470714213803771881909, 0.830223896278566929872032213968};
const double GLobx8[6] = {
    -0.871740148509606615337445761221, -0.5917001814331423021445107314,
    -0.20929921790247886876865726035,  0.20929921790247886876865726035,
    0.5917001814331423021445107314,    0.871740148509606615337445761221};
const double GLobx9[7] = {
    -0.899757995411460157312345244418, -0.677186279510737753445885427091,
    -0.363117463826178158710752068709, 0.0,
    0.363117463826178158710752068709,  0.677186279510737753445885427091,
    0.899757995411460157312345244418};
const double GLobx10[8] = {
    -0.919533908166458813828932660822, -0.73877386510550507500310617486,
    -0.477924949810444495661175092731, -0.165278957666387024626219765958,
    0.165278957666387024626219765958,  0.477924949810444495661175092731,
    0.73877386510550507500310617486,   0.919533908166458813828932660822};

const double *GLobx[9] = {0,      GLobx3, GLobx4, GLobx5, GLobx6,
                          GLobx7, GLobx8, GLobx9, GLobx10};
} // namespace glx

int w1e::readConfig(std::string file, int &qsz, int &R_max, int &l_max,
                    std::string &pot, std::string &quad_type,
                    std::string &quad_file, std::string &in_quad_layout,
                    std::string &pt_file) {
  YAML::Node settings = YAML::LoadFile(file);

  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "Core Potential:                              " << pot
            << std::endl;
  qsz = settings["Global_Settings"]["Outer_quadrature_size"].as<int>();
  std::cout << "No. of outer quadrature points:              " << qsz
            << std::endl;
  quad_type = settings["Global_Settings"]["Quadrature_type"].as<std::string>();
  std::cout << "Outer quadrature type:                       " << quad_type
            << std::endl;
  if (quad_type.compare("user-defined") == 0) {
    quad_file = settings["Global_Settings"]["quad_file"].as<std::string>();
    std::cout << "Quadrature read from file:                   " << quad_file
              << std::endl;
  }
  in_quad_layout =
      settings["Global_Settings"]["Inner_Layout"].as<std::string>();
  std::cout << "Layout of inner quadrature points:           " << in_quad_layout
            << std::endl;
  if (in_quad_layout.compare("user-defined") == 0) {
    pt_file = settings["Global_Settings"]["in_point_file"].as<std::string>();
    std::cout << "Quadrature read from file:                   " << pt_file
              << std::endl;
  }

  R_max = settings["Basis_Settings"]["R_max"].as<int>();
  std::cout << "Maximum R:                                   " << R_max
            << std::endl;
  l_max = settings["Basis_Settings"]["l_max"].as<int>();
  std::cout << "Maximum l:                                   " << l_max
            << std::endl;
  return 0;
}

// New decoupled implementation WIP
// Here gl_x holds the actual 0->R positions of the gl points
int prOuter(int n, int bo, int off, int ooff, std::vector<double> &kkn,
            std::vector<double> &q_x, std::vector<double> &Cf,
            std::vector<double> &p_out, std::vector<double> &p_der) {
  double Pl = 0;
  double Plp = 0;
  int qi = 0;
  size_t pi = 0;
  auto ksz = kkn.size();
  auto kn = gsl_vector_alloc(ksz);
  kn->size = ksz;
  kn->stride = 1;
  kn->data = kkn.data();
  kn->owner = 0;

  size_t nbreak = n + 2 - bo;
  auto dB = gsl_matrix_alloc(bo, 2);
  gsl_bspline_workspace *bw = gsl_bspline_alloc(bo, nbreak);
  bw->knots = kn;

  // Loop over knot regions
  for (auto i = bo - 1; i < n; ++i) {
    auto i1 = i + 1;
    pi = i;

    // While glq_pt is between current knots
    while (q_x[qi] > kkn[i] && q_x[qi] < kkn[i1]) {
      Pl = 0;
      Plp = 0;
      // Calculate the B-spline values on the fly
      gsl_bspline_deriv_eval_nonzero(q_x[qi], 1, dB, &pi, &pi, bw);
      for (auto j = 0; j < bo; ++j) {
        Pl += Cf[off + i1 - bo + j] * gsl_matrix_get(dB, j, 0);
        Plp += Cf[off + i1 - bo + j] * gsl_matrix_get(dB, j, 1);
      }
      p_out[ooff + qi] = Pl;
      p_der[ooff + qi] = Plp;
      ++qi;
    }
  }

  gsl_bspline_free(bw);
  gsl_matrix_free(dB);

  return 0;
}

int w1e::defaultPointLayout(int qsz, std::vector<uint8_t> &pq_dx, int &pti_sz) {
  pq_dx.resize(qsz);
  bool odd_flag = (qsz >> 0) & 1;
  int iv_sz = (qsz + 2 - odd_flag) / 2;
  std::fill(pq_dx.begin(), pq_dx.end(), 1);
  // default inner quadrature layout
  int idm = iv_sz - 1;
  int idp = idm + odd_flag;
  for (int i = 8; i > 1; --i) {
    pq_dx[idm] = i;
    pq_dx[idm - 1] = i;
    pq_dx[idp] = i;
    pq_dx[idp + 1] = i;
    idm -= 2;
    idp += 2;
  }

  pti_sz = std::reduce(pq_dx.begin(), pq_dx.end(), static_cast<int>(0));
  std::cout << "No. of inner points: " << pti_sz << "\n";

  return 0;
}

int w1e::userPointLayout(char ftype, std::string pt_file, int qsz,
                         std::vector<uint8_t> &pq_dx, int &pti_sz) {
  pq_dx.resize(qsz);

  switch (ftype) {
  case 't': {
    std::ifstream fl(pt_file);
    std::string temp;

    for (int i = 0; i < qsz; ++i) {
      std::getline(fl, temp);
      std::istringstream iss(temp);
      int st_num;
      iss >> st_num;
      pq_dx[i] = static_cast<uint8_t>(st_num);
    }
    break;
  }
  case 'b': {
    std::ifstream fl(pt_file, std::ios::in | std::ios::binary);
    fl.read(reinterpret_cast<char *>(pq_dx.data()), qsz);
    break;
  }
  }

  pti_sz = std::reduce(pq_dx.begin(), pq_dx.end(), static_cast<int>(0));
  std::cout << "No. of inner points: " << pti_sz << "\n";

  return 0;
}

int genRin(int qsz, int pti_sz, std::vector<double> &q_x,
           std::vector<uint8_t> &pq_dx, std::vector<double> &ri) {
  double dl, sl;

  ri.reserve(pti_sz);
  double rm1 = 0.0;
  for (int i = 0; i < qsz; ++i) {
    auto r1 = q_x[i];
    uint8_t p_num = pq_dx[i];
    dl = (r1 - rm1) * 0.5;
    sl = (r1 + rm1) * 0.5;
    for (uint8_t pt = 0; pt < p_num; ++pt) {
      ri.emplace_back(dl * glx::GLobx[p_num][pt] + sl);
    }
    rm1 = r1;
  }

  return 0;
}

int prInner(int n, int bo, int off, int ooff, std::vector<double> &kkn,
            std::vector<double> &Cf, std::vector<double> &ri,
            std::vector<double> &p_in) {
  double Pl = 0;
  int r_id = 0;
  size_t pi = 0;
  auto ksz = kkn.size();
  auto kn = gsl_vector_alloc(ksz);
  kn->size = ksz;
  kn->stride = 1;
  kn->data = kkn.data();
  kn->owner = 0;

  size_t nbreak = n + 2 - bo;
  auto B = gsl_vector_alloc(bo);
  gsl_bspline_workspace *bw = gsl_bspline_alloc(bo, nbreak);
  bw->knots = kn;

  // Loop over knot regions
  for (auto i = bo - 1; i < n; ++i) {
    auto i1 = i + 1;
    pi = i;

    // While glq_pt is between current knots
    while (ri[r_id] > kkn[i] && ri[r_id] < kkn[i1] &&
           r_id < static_cast<int>(ri.size())) {
      Pl = 0;
      // Calculate the B-spline values on the fly
      gsl_bspline_eval_nonzero(ri[r_id], B, &pi, &pi, bw);
      for (auto j = 0; j < bo; ++j) {
        Pl += Cf[off + i1 - bo + j] * B->data[j];
      }
      p_in[ooff + r_id] = Pl;
      ++r_id;
    }
  }

  gsl_bspline_free(bw);
  gsl_vector_free(B);

  return 0;
}

int w1e::genGaussLegendre(int qsz, int R_max, std::vector<double> &q_x,
                          std::vector<double> &q_w) {
  double R_half = R_max * 0.5;
  // generate GL nodes and weights scaled to [0,R_max]
  fastgl::QuadPair gl_i;
  for (int i = 1; i <= qsz; ++i) {
    gl_i = fastgl::GLPair(qsz, i);
    q_x[qsz - i] = R_half * (gl_i.x() + 1.0);
    q_w[qsz - i] = R_half * gl_i.weight;
  }

  return 0;
}

int w1e::readQuad(int qsz, int R_max, std::string quad_file, char type,
                  std::vector<double> &q_x, std::vector<double> &q_w) {
  double R_half = R_max * 0.5;
  double R_const{R_half};
  switch (type) {
  case 't': {
    std::ifstream fl(quad_file);
    std::string temp;

    for (int i = 0; i < qsz; ++i) {
      std::getline(fl, temp);
      std::istringstream iss(temp);
      double x_m1_p1, w_m1_p1;
      iss >> x_m1_p1 >> w_m1_p1;
      q_x[i] = R_half * (x_m1_p1 + 1.0);
      q_w[i] = R_half * w_m1_p1;
    }
    break;
  }
  case 'b': {
    std::ifstream fl(quad_file, std::ios::in | std::ios::binary);
    fl.read(reinterpret_cast<char *>(q_x.data()), qsz * sizeof(double));
    fl.read(reinterpret_cast<char *>(q_w.data()), qsz * sizeof(double));
    std::transform(q_x.begin(), q_x.end(), q_x.begin(),
                   [&R_const](auto &a) { return a * R_const + R_const; });
    std::transform(q_w.begin(), q_w.end(), q_w.begin(),
                   [&R_const](auto &a) { return a * R_const; });
    break;
  }
  }

  return 0;
}

int w1e::genWfn(std::string pot, int qsz, int pti_sz, int l_max,
                std::vector<double> &q_x, std::vector<double> &q_w,
                std::vector<uint8_t> &pq_dx) {
  int n, k, nen;
  std::string outfile_name, filename;
  std::unique_ptr<H5::H5File> outfile = nullptr, file = nullptr;
  std::unique_ptr<H5::DataSet> Po = nullptr, Pi = nullptr, Pidx = nullptr,
                               Cf_set = nullptr, glr = nullptr, glw = nullptr,
                               glri = nullptr;

  // read knots
  filename = pot + std::to_string(0) + ".h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  file->openAttribute("N").read(H5::PredType::NATIVE_INT32, &n);
  file->openAttribute("K").read(H5::PredType::NATIVE_INT32, &k);
  auto rset =
      std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Knots")));
  auto nkn = rset->getSpace().getSimpleExtentNpoints();

  std::vector<double> kkn(nkn);
  rset->read(kkn.data(), H5::PredType::NATIVE_DOUBLE); // read knots
  rset = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
  nen = rset->getSpace().getSimpleExtentNpoints();
  file->close();

  int nnpt = nen * qsz;

  std::vector<double> wfn_o(nnpt), wfn_p(nnpt);
  std::vector<double> Cf(n * nen);

  std::vector<double> ri;
  genRin(qsz, pti_sz, q_x, pq_dx, ri);
  std::vector<double> wfn_i(nen * pti_sz);
  hsize_t dimms_o[2] = {static_cast<hsize_t>(nen), static_cast<hsize_t>(qsz)};
  hsize_t dimms_i[2] = {static_cast<hsize_t>(nen),
                        static_cast<hsize_t>(pti_sz)};
  hsize_t dimm_q[1] = {static_cast<hsize_t>(qsz)};
  hsize_t dimm_i[1] = {static_cast<hsize_t>(pti_sz)};

#pragma omp parallel
  {
    for (auto l = 0; l <= l_max; ++l) {
#pragma omp single
      {
        filename = pot + std::to_string(l) + ".h5";
        file =
            std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
        Cf_set = std::make_unique<H5::DataSet>(
            H5::DataSet(file->openDataSet("Coeff")));
        Cf_set->read(&Cf[0], H5::PredType::NATIVE_DOUBLE);
        file->close();
      }

#pragma omp parallel for
      for (auto i = 0; i < nen; ++i) {
        prOuter(n, k, i * n, i * qsz, kkn, q_x, Cf, wfn_o, wfn_p);
        prInner(n, k, i * n, i * pti_sz, kkn, Cf, ri, wfn_i);
      }

#pragma omp single
      {
        outfile_name = pot + "_w1e" + std::to_string(l) + ".h5";
        outfile = std::make_unique<H5::H5File>(
            H5::H5File(outfile_name, H5F_ACC_TRUNC));
        // Write wave functions on outer points
        Po = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
            "Pr_o", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms_o))));
        Po->write(wfn_o.data(), H5::PredType::NATIVE_DOUBLE);
        // Write wave function derivatives on outer points
        Po = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
            "Pr_p", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms_o))));
        Po->write(wfn_p.data(), H5::PredType::NATIVE_DOUBLE);
        // Write wave functions on inner points
        Pi = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
            "Pr_i", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms_i))));
        Pi->write(wfn_i.data(), H5::PredType::NATIVE_DOUBLE);
        // Write inner quadrature point distribution
        Pidx = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
            "Pr_idx", H5::PredType::NATIVE_UCHAR, H5::DataSpace(1, dimm_q))));
        Pidx->write(pq_dx.data(), H5::PredType::NATIVE_UCHAR);
        // Write outer quadrature nodes on (0,R_max)
        glr = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
            "Qr_o", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, dimm_q))));
        glr->write(q_x.data(), H5::PredType::NATIVE_DOUBLE);
        // Write outer quadrature wieghts on (0,R_max)
        glw = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
            "Qw_o", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, dimm_q))));
        glw->write(q_w.data(), H5::PredType::NATIVE_DOUBLE);
        // Write inner quadrature nodes on (0,R_max)
        glri = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
            "Qr_i", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, dimm_i))));
        glri->write(ri.data(), H5::PredType::NATIVE_DOUBLE);
        outfile->close();
      }
    }
  }

  return 0;
}
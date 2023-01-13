#include "bsp_gsl.hpp"

int bsp::genKnots(int n, int k, int R_max, double fkn, char type,
                  std::vector<double> &kkn) {
  double par = 0.0, lin = 0.0;
  static constexpr double PI = 3.141592653589793238463;
  kkn.reserve(n + k);
  for (int i = 0; i < k; ++i) {
    kkn.emplace_back(0.0);
  }

  if (type == 'l') {
    par = R_max / (double)(n - k + 1);

    for (int i = k; i < n; ++i) {
      lin += par;
      kkn.emplace_back(lin);
    }

  } else if (type == 'e') {
    par = log(R_max / fkn) / (double)(n - k + 1);

    for (int i = k; i < n; ++i) {
      kkn.emplace_back(fkn * exp((double)(i - k) * par));
    }

  } else if (type == 's') {
    par = -log((2. / PI) * asin(fkn / R_max)) / log((double)(n - k + 1));

    for (int i = k; i < n; ++i) {
      kkn.emplace_back(
          R_max *
          sin((PI / 2.) * pow((double)(i - k + 1) / (double)(n - k + 1), par)));
    }

  } else {
    std::cout << "Invalid knot type !\n";
    return 1;
  }
  for (int i = n; i < n + k; ++i) {
    kkn.emplace_back(R_max);
  }
  return 0;
}

int bsp::genKnots(int n, int k, int R_max, std::string file, char type,
                  std::vector<double> &kkn) {
  std::string pt;
  int line_num = 0;
  kkn.reserve(n + k);
  for (int i = 0; i < k; ++i) {
    kkn.emplace_back(0.0);
  }

  if (type == 't') {
    std::ifstream fl(file);

    while (std::getline(fl, pt)) {
      ++line_num;
      kkn.emplace_back(std::stod(pt));
    }
    // assert(line_num==n);
  } else if (type == 'b') {
    std::ifstream fl(file, std::ios::in | std::ios::binary);
    fl.unsetf(std::ios::skipws);

    kkn.insert(std::end(kkn), std::istream_iterator<double>(fl),
               std::istream_iterator<double>());
  } else {
    std::cout << "Invalid knot file type !\n";
    return -1;
  }

  for (int i = n; i < n + k; ++i) {
    kkn.emplace_back(R_max);
  }
  return 0;
}

int bsp::wrKnotsH5(int n, int k, int R_max, double fkn, char type,
                   std::string file, std::vector<double> &kkn) {
  std::string tp_string;

  if (type == 'l') {
    tp_string = "linear";
  } else if (type == 'e') {
    tp_string = "exponential";
  } else if (type == 's') {
    tp_string = "sine";
  } else if (type == 'c') {
    tp_string = "user-defied";
  }

  auto fl = std::make_unique<H5::H5File>(H5::H5File(file, H5F_ACC_TRUNC));
  // Check if the file was opened
  if (!fl) {
    std::cerr << "# H5::H5File:: file couldn't opened: " << file << "\n";
    exit(-1);
  }

  hsize_t nKnots_d[1] = {kkn.size()};
  hsize_t att_space[1] = {1};
  hsize_t str_space[1] = {sizeof(tp_string) / sizeof(char)};
  H5::StrType str_ty(H5::PredType::C_S1, H5T_VARIABLE);

  H5::Attribute N = fl->createAttribute("N", H5::PredType::NATIVE_INT32,
                                        H5::DataSpace(1, att_space));
  H5::Attribute K = fl->createAttribute("K", H5::PredType::NATIVE_INT32,
                                        H5::DataSpace(1, att_space));
  H5::Attribute R = fl->createAttribute("R", H5::PredType::NATIVE_DOUBLE,
                                        H5::DataSpace(1, att_space));
  H5::Attribute TP =
      fl->createAttribute("Type", str_ty, H5::DataSpace(1, str_space));
  H5::Attribute F = fl->createAttribute("FKN", H5::PredType::NATIVE_DOUBLE,
                                        H5::DataSpace(1, att_space));

  H5::DataSet Knots = fl->createDataSet("Knots", H5::PredType::NATIVE_DOUBLE,
                                        H5::DataSpace(1, nKnots_d));

  N.write(H5::PredType::NATIVE_INT32, &n);
  K.write(H5::PredType::NATIVE_INT32, &k);
  R.write(H5::PredType::NATIVE_INT32, &R_max);
  TP.write(str_ty, &tp_string);
  F.write(H5::PredType::NATIVE_DOUBLE, &fkn);
  Knots.write(&kkn[0], H5::PredType::NATIVE_DOUBLE);

  return 0;
}

int bsp::splines(int n, int k, int glq_pt, std::vector<double> &gl_x,
                 std::vector<double> &kkn, std::vector<double> &splines,
                 std::vector<double> &splinesp) {
  int i1, B_sz = n * glq_pt * k;
  double dl, sl, x;
  auto kn = gsl_vector_alloc(kkn.size());
  kn->size = kkn.size();
  kn->stride = 1;
  kn->data = kkn.data();
  kn->owner = 0;

  size_t nbreak = n + 2 - k;
  auto dB = gsl_matrix_alloc(k, 2);
  gsl_bspline_workspace *bw = gsl_bspline_alloc(k, nbreak);
  bw->knots = kn;

  splines.reserve(B_sz);
  splinesp.reserve(B_sz);

  for (int i = 0; i < (k - 1) * glq_pt * k; ++i) {
    splines.emplace_back(0.0);
    splinesp.emplace_back(0.0);
  }

  size_t pi = 0;
  for (auto i = k - 1; i < n; ++i) {
    i1 = i + 1;
    dl = kkn[i1] - kkn[i];
    sl = kkn[i1] + kkn[i];

    for (int p = 0; p < glq_pt; ++p) {
      x = dl * 0.5 * gl_x[p] + sl * 0.5; // x-transformation

      pi = i;
      gsl_bspline_deriv_eval_nonzero(x, 1, dB, &pi, &pi, bw);
      for (int j = 0; j < k; ++j) {
        splines.emplace_back(gsl_matrix_get(dB, j, 0));
        splinesp.emplace_back(gsl_matrix_get(dB, j, 1));
      }
    }
  }
  gsl_bspline_free(bw);
  gsl_matrix_free(dB);

  return 0;
}

int bsp::splineInt(int n, int k, int glq_pt, std::vector<double> &gl_w,
                   std::vector<double> &gl_x, std::vector<double> &ov,
                   std::vector<double> &spl, std::vector<double> &kkn,
                   std::unique_ptr<ModelV> &Vptr) {
  int j_max, t_min, t_max, tm1;
  double ovlp, bsum, dl, sl, x;

  for (int i = 0; i < n; ++i) {
    j_max = std::min(i + k, n);
    for (int j = i; j < j_max; ++j) {
      t_min = std::max(k, j + 2);
      t_max = std::min(i + k + 1, n + 2);
      ovlp = 0.0;

      for (int t = t_min; t <= t_max; ++t) {
        tm1 = t - 1;
        dl = (kkn[t] - kkn[tm1]) * 0.5;
        sl = (kkn[t] + kkn[tm1]) * 0.5;
        bsum = 0.0;
        for (int p = 0; p < glq_pt; ++p) {
          x = dl * gl_x[p] + sl;
          bsum += gl_w[p] * spl[i + 1 - t + k + k * (p + (tm1)*glq_pt)] *
                  spl[j + 1 - t + k + k * (p + (tm1)*glq_pt)] * Vptr->V(x);
        }
        ovlp += dl * bsum;
      }
      // index into ov here needs Fortran (col major)
      ov[j * k + k - 1 + i - j] = ovlp;
    }
  }
  return 0;
}
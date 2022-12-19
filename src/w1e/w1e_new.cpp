#include "w1e_new.hpp"

namespace xw {

const double GLobx4[1] = {0.447213595499957939281834733746};
const double GLobx5[1] = {0.65465367070797714379829245625};
const double GLobx6[2] = {0.765055323929464692851002973959,
                          0.285231516480645096314150994041};
const double GLobx7[2] = {0.830223896278566929872032213968,
                          0.468848793470714213803771881909};
const double GLobx8[3] = {0.871740148509606615337445761221,
                          0.5917001814331423021445107314,
                          0.20929921790247886876865726035};
const double GLobx9[3] = {0.899757995411460157312345244418,
                          0.677186279510737753445885427091,
                          0.363117463826178158710752068709};
const double GLobx10[4] = {
    0.919533908166458813828932660822, 0.73877386510550507500310617486,
    0.477924949810444495661175092731, 0.165278957666387024626219765958};

const double *GLobx[11] = {0,      0,      0,      0,      GLobx4, GLobx5,
                           GLobx6, GLobx7, GLobx8, GLobx9, GLobx10};

const double GLobw3[2] = {0.333333333333333333333333333333,
                          1.333333333333333333333333333333};
const double GLobw4[2] = {0.166666666666666666666666666667,
                          0.83333333333333333333333333333};
const double GLobw5[3] = {0.1, 0.544444444444444444444444444444,
                          0.71111111111111111111111111111};
const double GLobw6[3] = {0.0666666666666666666666666666667,
                          0.378474956297846980316612808212,
                          0.554858377035486353016720525121};
const double GLobw7[4] = {
    0.047619047619047619047619047619, 0.27682604736156594801070040629,
    0.431745381209862623417871022281, 0.487619047619047619047619047619};
const double GLobw8[4] = {
    0.0357142857142857142857142857143, 0.21070422714350603938299206578,
    0.34112269248350436476424067711, 0.4124587946587038815670529714};
const double GLobw9[5] = {
    0.0277777777777777777777777777778, 0.16549536156080552504633972003,
    0.274538712500161735280705618579, 0.34642851097304634511513153214,
    0.371519274376417233560090702948};
const double GLobw10[5] = {
    0.0222222222222222222222222222222, 0.133305990851070111126227170755,
    0.224889342063126452119457821731,  0.292042683679683757875582257374,
    0.327539761183897456656510527917,
};

const double *Globw[11] = {0,      0,      0,      GLobw3, GLobw4, GLobw5,
                           GLobw6, GLobw7, GLobw8, GLobw9, GLobw10};
} // namespace xw

// New decoupled implementation WIP
// Here gl_x holds the actual 0->R positions of the gl points
int prIndp(int n, int bo, int off, int ooff, std::vector<double> &kkn,
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

int prInner(int n, int bo, int off, int ooff, std::vector<double> &kkn,
            std::vector<double> &q_x, std::vector<double> &Cf,
            std::vector<double> &p_out, std::vector<double> &p_der) {
  double Pl = 0;
  double Plp = 0;
  int qi = 0;
  size_t pi = 0;
  auto ksz = kkn.size();
  auto qsz = q_x.size();
  auto kn = gsl_vector_alloc(ksz);
  kn->size = ksz;
  kn->stride = 1;
  kn->data = kkn.data();
  kn->owner = 0;

  size_t nbreak = n + 2 - bo;
  auto B = gsl_vector_alloc(bo);
  gsl_bspline_workspace *bw = gsl_bspline_alloc(bo, nbreak);
  bw->knots = kn;

  double max_dx = 0.0;
  bool odd_flag = (qsz >> 0) & 1;

  int it_idx = (qsz - odd_flag) / 2;
  max_dx = q_x[it_idx + 1] - q_x[it_idx];

  int iv_sz = (qsz + 2 - odd_flag) / 2;
  std::vector<int> pt_q_dx(iv_sz, 3);

  int id = iv_sz - 1;
  for (int i = 10; i > 3; --i) {
    pt_q_dx[id] = i;
    pt_q_dx[id - 1] = i;
    id -= 2;
  }

  // Loop over knot regions
  for (auto i = bo - 1; i < n; ++i) {
    auto i1 = i + 1;
    pi = i;

    // While glq_pt is between current knots
    while (q_x[qi] > kkn[i] && q_x[qi] < kkn[i1]) {
      Pl = 0;
      Plp = 0;
      // Calculate the B-spline values on the fly
      gsl_bspline_eval_nonzero(q_x[qi], B, &pi, &pi, bw);
      for (auto j = 0; j < bo; ++j) {
        Pl += Cf[off + i1 - bo + j] * B->data[j];
      }
      p_out[ooff + qi] = Pl;
      ++qi;
    }
  }

  gsl_bspline_free(bw);
  gsl_vector_free(B);

  return 0;
}

int w1e::GenWfn(std::string pot, int glq_pt, int R_max, int l_max,
                std::string integrator) {
  int n, k, nen;
  std::string outfile_name, filename;
  std::unique_ptr<H5::H5File> outfile = nullptr, file = nullptr;
  std::unique_ptr<H5::DataSet> Po = nullptr, Pi = nullptr, Cf_set = nullptr;

  // read knots
  filename = pot + std::to_string(0) + ".h5";
  file = std::make_unique<H5::H5File>(H5::H5File(filename, H5F_ACC_RDONLY));
  file->openAttribute("N").read(H5::PredType::NATIVE_INT32, &n);
  file->openAttribute("K").read(H5::PredType::NATIVE_INT32, &k);
  auto rset =
      std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("Knots")));
  auto nkn = rset->getSpace().getSimpleExtentNpoints();

  std::vector<double> kkn(nkn);
  rset->read(&kkn[0], H5::PredType::NATIVE_DOUBLE); // read knots
  rset = std::make_unique<H5::DataSet>(H5::DataSet(file->openDataSet("En")));
  nen = rset->getSpace().getSimpleExtentNpoints();
  file->close();

  int npt = n * glq_pt;
  int nnpt = nen * npt;

  std::vector<double> wfn_o(nnpt), wfnp(nnpt);
  std::vector<double> Cf(n * nen);
  hsize_t dimms_o[2] = {(hsize_t)nen, (hsize_t)npt};

  // generate GL nodes and weights over B-splines support
  std::vector<double> gl_x(glq_pt);
  std::vector<double> gl_w(glq_pt);
  fastgl::QuadPair gl_i;
  for (int i = 1; i <= glq_pt; ++i) {
    gl_i = fastgl::GLPair(glq_pt, i);
    gl_x[glq_pt - i] = 1.0 + gl_i.x() * R_max;
    gl_w[glq_pt - i] = gl_i.weight;
  }

  if (integrator.compare("trapezoid") == 0) {
#pragma omp parallel
    {
      for (auto l = 0; l <= l_max; ++l) {
#pragma omp single
        {
          filename = pot + std::to_string(l) + ".h5";
          file = std::make_unique<H5::H5File>(
              H5::H5File(filename, H5F_ACC_RDONLY));
          Cf_set = std::make_unique<H5::DataSet>(
              H5::DataSet(file->openDataSet("Coeff")));
          Cf_set->read(&Cf[0], H5::PredType::NATIVE_DOUBLE);
          file->close();
        }

#pragma omp parallel for
        for (auto i = 0; i < nen; ++i) {
          prIndp(n, k, i * n, i * npt, kkn, gl_x, Cf, wfn_o, wfn_p);
        }
// HDF5 is not thread safe hence the omp singles
#pragma omp single
        {
          outfile_name = pot + "_w1e" + std::to_string(l) + ".h5";
          outfile = std::make_unique<H5::H5File>(
              H5::H5File(outfile_name, H5F_ACC_TRUNC));
          Po = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
              "Pr_o", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms_o))));
          Po->write(&wfn_o[0], H5::PredType::NATIVE_DOUBLE);
          outfile->close();

          outfile_name = pot + "_w1ep" + std::to_string(l) + ".h5";
          outfile = std::make_unique<H5::H5File>(
              H5::H5File(outfile_name, H5F_ACC_TRUNC));
          Po = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
              "Pr_p", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms_o))));
          Po->write(&wfnp[0], H5::PredType::NATIVE_DOUBLE);
          outfile->close();
        }
      }
    }
  } else if (integrator.compare("mixed") == 0) {

    std::vector<double> wfn_i(nnpt);

    std::vector<double> Ssp;
    bsp::Lob3Splines(n, k, glq_pt, gl_x, kkn, Ssp);

#pragma omp parallel
    {
      for (auto l = 0; l <= l_max; ++l) {
#pragma omp single
        {
          filename = pot + std::to_string(l) + ".h5";
          file = std::make_unique<H5::H5File>(
              H5::H5File(filename, H5F_ACC_RDONLY));
          Cf_set = std::make_unique<H5::DataSet>(
              H5::DataSet(file->openDataSet("Coeff")));
          Cf_set->read(&Cf[0], H5::PredType::NATIVE_DOUBLE);
          file->close();
        }

#pragma omp parallel for
        for (auto i = 0; i < nen; ++i) {
          prIndp(n, k, i * n, i * npt, kkn, gl_x, Cf, wfn_o, wfn_p);
          prInner(n, k, i * n, i * npt, kkn, gl_x, gl_ix, Cf, );
        }

#pragma omp single
        {
          outfile_name = pot + "_w1e" + std::to_string(l) + ".h5";
          outfile = std::make_unique<H5::H5File>(
              H5::H5File(outfile_name, H5F_ACC_TRUNC));
          Po = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
              "Pr_o", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms_o))));
          Po->write(&wfn_o[0], H5::PredType::NATIVE_DOUBLE);
          Pi = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
              "Pr_i", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms_o))));
          Pi->write(&wfn_i[0], H5::PredType::NATIVE_DOUBLE);
          outfile->close();

          outfile_name = pot + "_w1ep" + std::to_string(l) + ".h5";
          outfile = std::make_unique<H5::H5File>(
              H5::H5File(outfile_name, H5F_ACC_TRUNC));
          Po = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
              "Pr_p", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms_o))));
          Po->write(&wfnp[0], H5::PredType::NATIVE_DOUBLE);
          outfile->close();
        }
      }
    }
  } else {
    std::cout << "Integrator not supported, try 'trapezoid' or 'mixed'\n";
    return -1;
  }

  return 0;
}
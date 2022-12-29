#include "w1e_new.hpp"

namespace glx {
const double GLobx3[1] = {0.0};
const double GLobx4[2] = {-0.447213595499957939281834733746, 0.447213595499957939281834733746};
const double GLobx5[3] = {-0.65465367070797714379829245625, 0.0, 0.65465367070797714379829245625};
const double GLobx6[4] = {-0.765055323929464692851002973959, -0.285231516480645096314150994041,
                          0.285231516480645096314150994041, 0.765055323929464692851002973959};
const double GLobx7[5] = {-0.830223896278566929872032213968, -0.468848793470714213803771881909, 0.0,
                          0.468848793470714213803771881909, 0.830223896278566929872032213968};
const double GLobx8[6] = {-0.871740148509606615337445761221, -0.5917001814331423021445107314,
                          -0.20929921790247886876865726035, 0.20929921790247886876865726035,
                          0.5917001814331423021445107314, 0.871740148509606615337445761221};
const double GLobx9[7] = {-0.899757995411460157312345244418, -0.677186279510737753445885427091,
                          -0.363117463826178158710752068709, 0.0, 0.363117463826178158710752068709,
                          0.677186279510737753445885427091, 0.899757995411460157312345244418};
const double GLobx10[8] = {-0.919533908166458813828932660822, -0.73877386510550507500310617486,
                          -0.477924949810444495661175092731, -0.165278957666387024626219765958,
                          0.165278957666387024626219765958, 0.477924949810444495661175092731,
                          0.73877386510550507500310617486, 0.919533908166458813828932660822};

const double *GLobx[9] = {0.0, GLobx3, GLobx4, GLobx5,
                           GLobx6, GLobx7, GLobx8, GLobx9, GLobx10};
} // namespace glx

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

int genPtaR(int R_max, int qsz, std::vector<double> &q_x, std::vector<uint8_t> &pq_dx, 
            int &pti_sz, std::vector<double> &ri) {
  double dl, sl;
  bool odd_flag = (qsz >> 0) & 1;
  int iv_sz = (qsz + 2 - odd_flag) / 2;
  pq_dx.resize(qsz + 1);
  std::fill(pq_dx.begin(), pq_dx.end(), 1);

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

  pti_sz = std::accumulate(pq_dx.begin(), pq_dx.end(), static_cast<int>(0));
  std::cout << "pti_sz: " << pti_sz << "\n";

  ri.reserve(pti_sz);
  ri.emplace_back(q_x[0] * 0.5);
  for (int i=0; i<qsz-1; ++i) {
    auto i1 = i+1;
    auto p_num = pq_dx[i1];
    dl = (q_x[i1] - q_x[i]) * 0.5;
    sl = (q_x[i1] + q_x[i]) * 0.5;
    for (int pt=0; pt<p_num; ++pt) {
      ri.emplace_back(dl * glx::Globx[p_num][pt] + sl);
    }
  }
  dl = (R_max - q_x[i]) * 0.5;
  sl = (R_max + q_x[i]) * 0.5;
  for (int pt=0; pt<pq_dx[qsz]; ++pt) {
    ri.emplace_back(dl * glx::Globx[pq_dx[qsz]][pt] + sl);
  }

  return 0;
}

int prInner(int n, int bo, int nen, int off, int ooff, std::vector<double> &kkn,
            std::vector<double> &Cf, std::vector<uint8_t> &pq_dx, 
            std::vector<double> &ri, std::vector<double> &p_in) {
  double Pl = 0;
  double Plp = 0;
  int r_id = 0;
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

  // Loop over knot regions
  for (auto i = bo - 1; i < n; ++i) {
    auto i1 = i + 1;

    // While glq_pt is between current knots
    while (ri[r_id] > kkn[i] && ri[r_id] < kkn[i1]) {
      // Calculate the B-spline values on the fly
      pi=i;
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

int w1e::GenWfn(std::string pot, int glq_pt, int R_max, int l_max,
                std::string integrator) {
  int n, k, nen;
  std::string outfile_name, filename;
  std::unique_ptr<H5::H5File> outfile = nullptr, file = nullptr;
  std::unique_ptr<H5::DataSet> Po = nullptr, Pi = nullptr, 
      Pidx = nullptr, Cf_set = nullptr;

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
  fastgl::QuadPair gl_i;
  for (int i = 1; i <= glq_pt; ++i) {
    gl_i = fastgl::GLPair(glq_pt, i);
    gl_x[glq_pt - i] = 1.0 + gl_i.x() * R_max;
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
          Po->write(&wfn_p[0], H5::PredType::NATIVE_DOUBLE);
          outfile->close();
        }
      }
    }
  } else if (integrator.compare("mixed") == 0) {

    std::vector<uint8_t> pq_dx;
    std::vector<double> ri;
    int pti_sz;
    genPtaR(R_max, glq_pt, q_x, pq_dx, pti_sz, ri);
    std::vector<double> wfn_i(nen*pti_sz);

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
          prInner(n, k, i * n, i * npt, kkn, gl_x, Cf, pq_dx, wfn_i);
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
          Pidx = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
              "Pr_idx", H5::PredType::NATIVE_UCHAR, H5::DataSpace(1, static_cast<hsize_t>(pti_sz)))));
          Pidx->write(&pq_dx[0], H5::PredType::NATIVE_UCHAR);
          glr = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
              "R_o", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, static_cast<hsize_t>(glq_pt)))));
          glri = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
              "R_i", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, static_cast<hsize_t>(pti_sz)))));
          outfile->close();

          outfile_name = pot + "_w1ep" + std::to_string(l) + ".h5";
          outfile = std::make_unique<H5::H5File>(
              H5::H5File(outfile_name, H5F_ACC_TRUNC));
          Po = std::make_unique<H5::DataSet>(H5::DataSet(outfile->createDataSet(
              "Pr_p", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms_o))));
          Po->write(&wfn_p[0], H5::PredType::NATIVE_DOUBLE);
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
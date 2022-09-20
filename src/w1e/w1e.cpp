#include "w1e.hpp"

int w1e::ReadConfig(std::string file, int &glq_pt, int &l_max, std::string &pot,
                    std::string &integrator) {
  YAML::Node settings = YAML::LoadFile(file);

  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "Core Potential:                              " << pot
            << std::endl;
  glq_pt = settings["Global_Settings"]["GL_quad_points"].as<int>();
  std::cout << "No. of GL-quadrature points between knots:   " << glq_pt
            << std::endl;

  l_max = settings["Basis_Settings"]["l_max"].as<int>();
  std::cout << "Maximum l:                                   " << l_max
            << std::endl;

  integrator = settings["R12_Settings"]["integrator"].as<std::string>();
  std::cout << "Integration scheme (inner integral):         " << integrator
            << std::endl;
  return 0;
}

int Pr(int n, int bo, int glq_pt, int off, int ooff, std::vector<double> &kkn,
       std::vector<double> &gl_x, std::vector<double> &Bsp,
       std::vector<double> &Cf, std::vector<double> &p_out) {
  double Pl = 0;

  for (auto i = bo - 1; i < n; ++i) {
    auto i1 = i + 1;
    auto i1bo = (i1 - bo) * glq_pt;

    for (auto p = 0; p < glq_pt; ++p) {
      Pl = 0;
      for (auto j = 0; j < bo; ++j) {
        Pl += Cf[off + i1 - bo + j] * Bsp[j + bo * (p + i * glq_pt)];
      }
      p_out[ooff + i1bo + p] = Pl;
    }
  }

  return 0;
}

int Prlob3(int n, int bo, int glq_pt, int off, int ooff,
           std::vector<double> &kkn, std::vector<double> &gl_x,
           std::vector<double> &Ssp, std::vector<double> &Cf,
           std::vector<double> &p_in) {
  double rm1 = 0.0, r = 0.0, rlob = 0.0;
  double Pl = 0;
  double dl, sl;

  for (auto i = bo - 1; i < n; ++i) {
    auto i1 = i + 1;
    auto i1bo = (i1 - bo) * glq_pt;
    dl = (kkn[i1] - kkn[i]) * 0.5;
    sl = (kkn[i1] + kkn[i]) * 0.5;

    for (auto p = 0; p < glq_pt; ++p) {
      r = dl * gl_x[p] + sl;
      rlob = (rm1 + r) * 0.5;
      Pl = 0;

      for (auto j = 0; j < bo; ++j) {
        Pl += Cf[off + i1 - bo + j - (kkn[i] > rlob)] *
              Ssp[j + bo * (p + i * glq_pt)];
      }

      p_in[ooff + i1bo + p] = Pl;
      rm1 = r;
    }
  }

  return 0;
}

int Prlob4(int n, int bo, int glq_pt, int off, int ooff,
           std::vector<double> &kkn, std::vector<double> &gl_x,
           std::vector<double> &Ssp, std::vector<double> &Cf,
           std::vector<double> &p_in) {
  double rm1 = 0.0, r = 0.0, rloba = 0.0, rlobb = 0.0;
  double Pla = 0, Plb = 0;
  double dl, sl, dlob, slob;
  int i1, i1bo, i2bo, sp, ai, bi, ibo1j;

  constexpr double Lob4p = 0.4472135954999579392818347e0;

  for (auto i = bo - 1; i < n; ++i) {
    i1 = i + 1;
    i1bo = (i1 - bo) * glq_pt;
    i2bo = 2 * i1bo;
    dl = (kkn[i1] - kkn[i]) * 0.5;
    sl = (kkn[i1] + kkn[i]) * 0.5;

    sp = 0;
    for (auto p = 0; p < glq_pt; ++p, sp += 2) {
      r = dl * gl_x[p] + sl;
      dlob = (r - rm1) * 0.5 * Lob4p;
      slob = (r + rm1) * 0.5;
      rloba = -dlob + slob;
      rlobb = dlob + slob;
      Pla = 0, Plb = 0;

      for (auto j = 0; j < bo; ++j) {
        ibo1j = i1 - bo + j;

        ai = ibo1j - (kkn[i] > rloba);
        bi = ibo1j - (kkn[i] > rlobb);
        Pla += Cf[off + ai] * Ssp[j + bo * (sp + i * 2 * glq_pt)];
        Plb += Cf[off + bi] * Ssp[j + bo * (sp + 1 + i * 2 * glq_pt)];
      }
      p_in[ooff + i2bo + sp] = Pla;
      p_in[ooff + i2bo + sp + 1] = Plb;
      rm1 = r;
    }
  }

  return 0;
}

int w1e::GenWfn(std::string pot, int glq_pt, int l_max,
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
    gl_x[glq_pt - i] = gl_i.x();
    gl_w[glq_pt - i] = gl_i.weight;
  }

  // generate B-splines
  std::vector<double> Bsp, Psp;
  bsp::Splines(n, k, glq_pt, gl_x, kkn, Bsp);
  bsp::SplinesP(n, k, glq_pt, gl_x, kkn, Psp);

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
          Pr(n, k, glq_pt, i * n, i * npt, kkn, gl_x, Bsp, Cf, wfn_o);
          Pr(n, k, glq_pt, i * n, i * npt, kkn, gl_x, Psp, Cf, wfnp);
        }

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
  } else if (integrator.compare("mixed") == 0 ||
             integrator.compare("glob3") == 0) {

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
          Pr(n, k, glq_pt, i * n, i * npt, kkn, gl_x, Bsp, Cf, wfn_o);
          Pr(n, k, glq_pt, i * n, i * npt, kkn, gl_x, Psp, Cf, wfnp);
          Prlob3(n, k, glq_pt, i * n, i * npt, kkn, gl_x, Ssp, Cf, wfn_i);
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
  } else if (integrator.compare("glob4") == 0) {

    hsize_t dimms_i[2] = {(hsize_t)nen, (hsize_t)2 * npt};
    std::vector<double> wfn_i(2 * nnpt);

    std::vector<double> Ssp;
    bsp::Lob4Splines(n, k, glq_pt, gl_x, kkn, Ssp);

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
          Pr(n, k, glq_pt, i * n, i * npt, kkn, gl_x, Bsp, Cf, wfn_o);
          Pr(n, k, glq_pt, i * n, i * npt, kkn, gl_x, Psp, Cf, wfnp);
          Prlob4(n, k, glq_pt, i * n, i * 2 * npt, kkn, gl_x, Ssp, Cf, wfn_i);
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
              "Pr_i", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dimms_i))));
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
  }

  return 0;
}
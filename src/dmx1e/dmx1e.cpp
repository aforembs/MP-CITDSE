#include "dmx1e.h"

int dmx1e::ReadConfig(std::string file, int &glq_pt,
                      std::string &pot, int &l_max,
                      char &gauge) {
  YAML::Node settings = YAML::LoadFile(file);

  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "Core Potential:                              "
            << pot << std::endl;
  glq_pt = settings["Global_Settings"]["GL_quad_points"].as<int>();
  std::cout << "No. of GL-quadrature points between knots:   "
            << glq_pt << std::endl;
  gauge = settings["Global_Settings"]["gauge"].as<char>();
  std::cout << "Gauge type ('l' length/'v' velocity):        "
            << gauge << std::endl;

  l_max = settings["Basis_Settings"]["l_max"].as<int>();
  std::cout << "Maximum one electron angular momentum:       "
            << l_max << std::endl;            
  return 0;
}

int dmx1e::GenDipole(std::string cpot, int l_max, int glq_pt, char gauge) {
  int n=0;   // no. of points
  int bo=0;   // max B-spline order
  int nkn=0; // no. of knots
  int lp1=0;
  double kl, t_ab;
  std::vector<double> kkn;
  std::string filename;

  H5::DataSpace cspace;
  H5::DataSpace memspace;
  auto D_set = std::unique_ptr<H5::DataSet>();
  hsize_t d_dim[2];

  // read knots
  filename = cpot + std::to_string(0) + ".h5";
  auto file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  file->openAttribute("N").read(H5::PredType::NATIVE_INT32, &n);
  file->openAttribute("K").read(H5::PredType::NATIVE_INT32, &bo);
  auto rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Knots")));
  nkn = rset->getSpace().getSimpleExtentNpoints();

  kkn.reserve(nkn);
  rset->read(&kkn[0], H5::PredType::NATIVE_DOUBLE); //read knots
  file->close();

  auto lc_sz = n*n*glq_pt;
  std::vector<double> wfn(lc_sz*(l_max+1));
  std::vector<double> D(n*n);

  d_dim[0]=n;
  d_dim[1]=n;
  // read wavefunctions for all l
  for(int l=0; l<=l_max; ++l) {
    filename = cpot + "_w1e" + std::to_string(l) + ".h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Pr_o")));
    rset->read(&wfn[l*lc_sz], H5::PredType::NATIVE_DOUBLE);
    file->close();
  }

  // generate GL nodes and weights over B-splines support
  std::vector<double> gl_x(glq_pt);
  std::vector<double> gl_w(glq_pt);
  fastgl::QuadPair gl_i;
  for(int i=1; i<=glq_pt; ++i) {
    gl_i = fastgl::GLPair(glq_pt, i);
    gl_x[glq_pt-i] = gl_i.x(); 
    gl_w[glq_pt-i] = gl_i.weight;
  }

  std::vector<double> wfnp;

  switch (gauge) {
    case 'v':
      wfnp.reserve(lc_sz);
      #pragma omp parallel
      {
        for(int l=0; l<l_max; ++l) {
          #pragma omp single
          {
            filename = cpot + "_w1ep" + std::to_string(l+1) + ".h5";
            file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
            rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Pr_p")));
            rset->read(&wfnp[0], H5::PredType::NATIVE_DOUBLE);
            file->close();
          }
          lp1=l+1;
          kl=sqrt((double)(lp1*lp1)/(4.0*lp1*lp1-1));
          for(int n2=0; n2<n; ++n2) {
            #pragma omp barrier
            #pragma omp for private(t_ab)
            for(int n1=0; n1<n; ++n1) {
              t_ab = kl*tvelGL(n,glq_pt,bo,lc_sz,n1,l,n2,lp1,gl_w,gl_x,kkn,wfn,wfnp);
              D[n1+n*n2]=t_ab;
            }
          }
          #pragma omp single
          {
          filename = cpot + std::to_string(l) + std::to_string(l+1) + gauge + ".h5";
          file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_TRUNC));
          D_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet(
                    "d_if", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, d_dim))));
          D_set->write(&D[0], H5::PredType::NATIVE_DOUBLE);
          file->close();
          }
        }
      }
      wfnp.clear();
      break;
    case 'l':
      #pragma omp parallel
      {
        for(int l=0; l<l_max; ++l) {
          lp1=l+1;
          kl=sqrt((double)(lp1*lp1)/(4.0*lp1*lp1-1));
          for(int n2=0; n2<n; ++n2) {
            #pragma omp barrier
            #pragma omp for private(t_ab)
            for(int n1=0; n1<n; ++n1) {
              t_ab = kl*tlenGL(n,glq_pt,bo,lc_sz,n1,l,n2,lp1,gl_w,gl_x,kkn,wfn);
              D[n1+n*n2]=t_ab;
            }
          }
          #pragma omp single
          {
          filename = cpot + std::to_string(l) + std::to_string(l+1) + gauge + ".h5";
          file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_TRUNC));
          D_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet(
                    "d_if", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, d_dim))));
          D_set->write(&D[0], H5::PredType::NATIVE_DOUBLE);
          file->close();
          }
        }
      }
      break;
  }

  return 0;
}
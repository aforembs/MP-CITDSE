#include "r12.h"
#include "time_tst.h"

int r_12::ReadConfig(std::string file, int &glq_pt,
                    std::string &pot, int &L_max,
                    char &gauge, std::string &integrator) {
  YAML::Node settings = YAML::LoadFile(file);

  pot = settings["Global_Settings"]["potential"].as<std::string>();
  std::cout << "Core Potential:                              "
            << pot << std::endl;
  glq_pt = settings["Global_Settings"]["GL_quad_points"].as<int>();
  std::cout << "No. of GL-quadrature points between knots:   "
            << glq_pt << std::endl;
  L_max = settings["Global_Settings"]["L_max"].as<int>();
  std::cout << "Maximum total two electron angular momentum: "
            << L_max << std::endl;
  gauge = settings["Global_Settings"]["gauge"].as<char>();
  std::cout << "Gauge type ('l' length/'v' velocity):        "
            << gauge << std::endl;

  integrator = settings["R12_Settings"]["integrator"].as<std::string>();
  std::cout << "Integration scheme (inner integral):         "
            << integrator << std::endl;
  return 0;
}

int Rpowk(std::vector<double> &r_out,int n, int bo, 
          int glq_pt, std::vector<double> &gl_x, 
          std::vector<double> &kkn, int k_max) {
  auto nqpt=n*glq_pt;
  r_out.reserve((k_max+1)*nqpt);
  auto km1=0;
  auto i1=0;
  
  double dl, sl, ri;
  auto ibo1bo=0, iqpt=0;
  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    dl = (kkn[i1] - kkn[i])*0.5;
    sl = (kkn[i1] + kkn[i])*0.5;
    ibo1bo = (i1-bo)*glq_pt;

    for(auto p=0; p<glq_pt; ++p){
      ri = dl*gl_x[p] + sl;
      r_out[p+ibo1bo+nqpt] = ri;
      r_out[p+ibo1bo] = 1.0;
    }
  }

  for(auto k=2; k<=k_max; ++k) {
    km1=k-1;
    for(auto i=0; i<n+1-bo; ++i) {
      iqpt = i*glq_pt;
      for(auto p=0; p<glq_pt; ++p){
        r_out[p+iqpt+k*nqpt] = r_out[p+iqpt+nqpt]*r_out[p+iqpt+km1*nqpt];
      }
    }
  }

  return 0;
}

int Rpowk(std::vector<double> &r_out, std::vector<double> &r_in, 
          int n, int bo, int glq_pt, std::vector<double> &gl_x, 
          std::vector<double> &kkn, int k_max) {
  auto nqpt=n*glq_pt;
  r_out.reserve((k_max+1)*nqpt);
  r_in.reserve((k_max+1)*nqpt);
  auto km1=0;
  auto i1=0;
  double rm1 = 0.0;
  
  double rlm, dl, sl, ri;
  auto ibo1bo=0, iqpt=0;
  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    dl = (kkn[i1] - kkn[i])*0.5;
    sl = (kkn[i1] + kkn[i])*0.5;
    ibo1bo = (i1-bo)*glq_pt;

    for(auto p=0; p<glq_pt; ++p){
      ri = dl*gl_x[p] + sl;
      rlm = (rm1+ri)*0.5;
      r_out[p+ibo1bo+nqpt] = ri;
      r_out[p+ibo1bo] = 1.0;
      r_in[p+ibo1bo+nqpt] = rlm;
      r_in[p+ibo1bo] = 1.0;
      rm1=ri;
    }
  }

  for(auto k=2; k<=k_max; ++k) {
    km1=k-1;
    for(auto i=0; i<n+1-bo; ++i) {
      iqpt = i*glq_pt;
      for(auto p=0; p<glq_pt; ++p){
        r_out[p+iqpt+k*nqpt] = r_out[p+iqpt+nqpt]*r_out[p+iqpt+km1*nqpt];
        r_in[p+iqpt+k*nqpt] = r_in[p+iqpt+nqpt]*r_in[p+iqpt+km1*nqpt];
      }
    }
  }

  return 0;
}

int Rpowklob4(std::vector<double> &r_out, std::vector<double> &r_in, 
          int n, int bo, int glq_pt, std::vector<double> &gl_x, 
          std::vector<double> &kkn, int k_max) {
  auto nqpt=n*glq_pt;
  auto nqpt2=2*nqpt;
  r_out.reserve((k_max+1)*nqpt);
  r_in.reserve((k_max+1)*nqpt2);
  auto km1=0;
  auto i1=0;
  double rm1 = 0.0;
  constexpr double Lob4p = 0.4472135954999579392818347e0;


  double dl, sl, ri, dlob, slob;
  double rloba, rlobb;
  auto ibo1gl=0, ibo1gl2=0, iqpt=0, iqpt2=0;
  auto sp=0;
  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    dl = (kkn[i1] - kkn[i])*0.5;
    sl = (kkn[i1] + kkn[i])*0.5;
    ibo1gl = (i1-bo)*glq_pt;
    ibo1gl2=2*ibo1gl;

    sp=0;
    for(auto p=0; p<glq_pt; ++p, sp+=2){
      ri = dl*gl_x[p] + sl;
      dlob = (ri-rm1)*0.5*Lob4p;
      slob = (ri+rm1)*0.5;
      rloba = -dlob + slob;
      rlobb = dlob + slob;
      r_out[p+ibo1gl+nqpt] = ri;
      r_out[p+ibo1gl] = 1.0;
      r_in[sp+ibo1gl2+nqpt2] = rloba;
      r_in[sp+1+ibo1gl2+nqpt2] = rlobb;
      r_in[sp+ibo1gl2] = 1.0;
      r_in[sp+1+ibo1gl2] = 1.0;
      rm1=ri;
    }
  }

  for(auto k=2; k<=k_max; ++k) {
    km1=k-1;
    for(auto i=0; i<n+1-bo; ++i) {
      iqpt = i*glq_pt;
      iqpt2= 2*iqpt;
      sp=0;
      for(auto p=0; p<glq_pt; ++p, sp+=2){
        r_out[p+iqpt+k*nqpt] = r_out[p+iqpt+nqpt]*r_out[p+iqpt+km1*nqpt];
        r_in[sp+iqpt2+k*nqpt2] = r_in[sp+iqpt2+nqpt2]*r_in[sp+iqpt2+km1*nqpt2];
        r_in[sp+1+iqpt2+k*nqpt2] = r_in[sp+1+iqpt2+nqpt2]*r_in[sp+1+iqpt2+km1*nqpt2];
      }
    }
  }

  return 0;
}

int Rpowk(std::vector<double> &r_out, std::vector<double> &r_in, 
          int n, int bo, int glq_pt, std::vector<double> &gl_xo, 
          std::vector<double> &gl_xi, std::vector<double> &kkn, int k_max) {
  r_out.reserve((k_max+1)*n*bo);
  r_in.reserve((k_max+1)*n*bo);
  auto km1=0;

  double rm1 = 0.0;
  auto nbo=n*bo;
  double rlm, dl, sl, ri;
  auto ibo1bo=0;
  for(auto i=bo-1; i<n; ++i) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    ibo1bo = (i-bo+1)*bo;

    for(auto p=0; p<bo; ++p){
      ri = dl*gl_xo[p] + sl;
      rlm = (rm1+ri)*0.5;
      r_out[p+ibo1bo+nbo] = ri;
      r_out[p+ibo1bo] = 1.0;
      r_in[p+ibo1bo+nbo] = rlm;
      r_in[p+ibo1bo] = 1.0;
      rm1=ri;
    }
  }

  for(auto k=2; k<=k_max; ++k) {
    km1=k-1;
    rm1 = 0.0;
    for(auto i=bo-1; i<n; ++i) {
      dl = (kkn[i+1] - kkn[i])*0.5;
      sl = (kkn[i+1] + kkn[i])*0.5;
      ibo1bo = (i-bo+1)*bo;

      for(auto p=0; p<bo; ++p){
        ri = dl*gl_xo[p] + sl;
        rlm = (rm1+ri)*0.5;
        r_out[p+ibo1bo+k*nbo] = ri*r_out[p+ibo1bo+km1*nbo];
        r_in[p+ibo1bo+k*nbo] = rlm*r_in[p+ibo1bo+km1*nbo];
        rm1=ri;
      }
    }
  }

  return 0;
}

int r_12::R12MM(std::string cpot, int L_max, int glq_pt, std::string dir) {
  bool min_dir=0, min_exc=0;
  double Y_norm=0.0;
  std::vector<double> v_mat;
  double sum_k=0.0;
  idx4 e12, e12p;

  int n=0;   // no. of points
  int bo=0;   // max B-spline order
  int nkn=0; // no. of knots
  std::vector<double> kkn;

  int L_sz=0, v_sz=0;
  std::string filename;
  std::string outfile_name;

  // HDF5 defines
  std::unique_ptr<H5::H5File> outfile=nullptr;
  std::unique_ptr<H5::DataSet> V_set=nullptr;
  std::unique_ptr<H5::DataSet> L_set=nullptr;
  std::vector<idx4> L_idx;
  H5::DataSpace cspace;
  hsize_t offset[2]={0,0}, count[2], stride[2]={1,1}, block[2]={1,1};
  hsize_t dimms[2], v_dim[1];
  H5::DataSpace memspace;

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

  int ncf, sym, max_N=0, l1_m=0, l2_m=0;
  std::vector<cfg::line> cfgs;
  cfg::line max_line;

  for(auto li=0; li<=L_max; ++li) {
    cfg::ReadCfg(dir, li, sym, ncf, cfgs);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
    max_N = std::max(max_N,max_line.n2max);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.l2 < b.l2; });
    l2_m = std::max(l2_m,max_line.l2);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });
    l1_m = std::max(l1_m,max_line.l1);
  }
  auto k_max = l1_m+l2_m;
  auto lc_sz = max_N*n*glq_pt;

  // reserve space for wfns
  auto e1_lm = std::max(l1_m,l2_m);
  std::vector<double> wfn_o(lc_sz*(e1_lm+1));
  std::vector<double> wfn_i(lc_sz*(e1_lm+1));

  count[0] = max_N;
  count[1] = n*glq_pt;
  dimms[0] = count[0]; 
  dimms[1] = count[1];
  memspace.setExtentSimple(2, dimms, NULL);
  // read wavefunctions for all l
  for(int l=0; l<=e1_lm; ++l) {
    filename = cpot + "_w1e" + std::to_string(l) + ".h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Pr_o")));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(&wfn_o[l*lc_sz], H5::PredType::NATIVE_DOUBLE, memspace, cspace);
    rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Pr_i")));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(&wfn_i[l*lc_sz], H5::PredType::NATIVE_DOUBLE, memspace, cspace);
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

  // Precalculate r^k
  std::vector<double> rk, rk_mid;
  Rpowk(rk, rk_mid, n, bo, glq_pt, gl_x, kkn, k_max);

  omp_lock_t copylock;
  omp_init_lock(&copylock);
  uint64_t st_time;

  wig_table_init(2*(k_max+2),6);
  
  #pragma omp parallel private(e12p)
  {
  wig_thread_temp_init(2*(k_max+2));

  for(int L=0; L<=L_max; ++L) {
    
    // Read n1l1;n2l2 indices for NL states
    #pragma omp single
    {
    filename = cpot + std::to_string(L) + "idx.h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
    L_sz = L_set->getSpace().getSimpleExtentNpoints()/4;
    L_idx.reserve(L_sz);
    L_set->read(&L_idx[0], H5::PredType::NATIVE_INT32);
    file->close();

    std::cout << "L: " << L <<" L_sz: "<< L_sz << "\n";

    v_sz = L_sz*(L_sz+1)/2;
    v_mat.reserve(v_sz); // error

    st_time = GetTimeMs64();
    }

    for(int NL2=0; NL2<L_sz; ++NL2) {
      //set n1'l1';n2'l2'
      e12p =L_idx[NL2];
      bool eqvp = e12p.n1==e12p.n2 && e12p.l1==e12p.l2;

      #pragma omp barrier
      #pragma omp for private(e12,Y_norm,sum_k,min_dir,min_exc)
      for(int NL1=NL2; NL1<L_sz; ++NL1) {
        //set n1l1;n2l2
        e12 = L_idx[NL1];

        // sqrt([la][lc][lb][ld])
        Y_norm = sqrt((2*e12.l1+1)*(2*e12p.l1+1)*(2*e12.l2+1)*(2*e12p.l2+1));
        sum_k=0.0;
        min_dir = ((L+e12.l2+e12p.l1) >> 0) & 1; // check if L+lb+lc is odd
        min_exc = ((L+e12.l1+e12p.l1) >> 0) & 1;
        for(int k=0; k<=k_max; ++k) {
          /* for direct check if:
            (-)^{L+lb+lc}=(-)^{L+la+ld},
            |la-lc| <= k <= la+lc,
            |lb-ld| <= k <= lb+ld */
          if(min_dir ==(((L+e12.l1+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l1)<=k) && (k<=e12.l1+e12p.l1)) &&
             ((abs(e12.l2-e12p.l2)<=k) && (k<=e12.l2+e12p.l2)) &&
             !(((k+e12.l1+e12p.l1) >> 0) & 1) &&
             !(((k+e12.l2+e12p.l2) >> 0) & 1) &&
             !std::isinf(1.0/rk_mid[(k+1)*2*n*glq_pt])) {
            sum_k += pow(-1,min_dir)*FsltrMM(k, n, bo, glq_pt, lc_sz,
                      e12.n1, e12.l1, e12.n2, e12.l2, e12p.n1, e12p.l1,
                      e12p.n2, e12p.l2, gl_w, gl_x, kkn, rk, rk_mid, wfn_o, wfn_i)
                      *wig3jj(2*e12.l1,2*k,2*e12p.l1,0,0,0)
                      *wig3jj(2*e12.l2,2*k,2*e12p.l2,0,0,0)
                      *wig6jj(2*e12p.l1,2*e12p.l2,2*L,2*e12.l2,2*e12.l1,2*k);
          }
          /* for exchange check if:
            (-)^{L+la+lc}=(-)^{L+lb+ld},
            |la-ld| <= k <= la+ld,
            |lc-lb| <= k <= lc+lb */
          if(min_exc ==(((L+e12.l2+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l2)<=k) && (k<=e12.l1+e12p.l2)) &&
             ((abs(e12.l2-e12p.l1)<=k) && (k<=e12.l2+e12p.l1)) &&
             !(((k+e12.l1+e12p.l2) >> 0) & 1) &&
             !(((k+e12.l2+e12p.l1) >> 0) & 1) &&
             !std::isinf(1.0/rk_mid[(k+1)*2*n*glq_pt])) {
            sum_k += pow(-1,min_exc)*FsltrMM(k, n, bo, glq_pt, lc_sz,
                      e12.n1, e12.l1, e12.n2, e12.l2, e12p.n2, e12p.l2,
                      e12p.n1, e12p.l1, gl_w, gl_x, kkn, rk, rk_mid, wfn_o, wfn_i)
                      *wig3jj(2*e12.l1,2*k,2*e12p.l2,0,0,0)
                      *wig3jj(2*e12.l2,2*k,2*e12p.l1,0,0,0)
                      *wig6jj(2*e12p.l1,2*e12p.l2,2*L,2*e12.l1,2*e12.l2,2*k);
          }
        }

        // write symmetric V_12 as upper triangular
        double v12 = pow(-1,e12.l1+e12.l2)*Y_norm*sum_k;

        if(e12.n1==e12.n2 && e12.l1==e12.l2) 
          v12 *= 0.7071067811865475244008444e0;

        if(eqvp) 
          v12 *= 0.7071067811865475244008444e0;

        v_mat[(2*L_sz-NL2-1)*NL2/2 + NL1] = v12;
      }
    }

    #pragma omp single
    {
    std::cout << "loop time: " << ((double)(GetTimeMs64()-st_time)/1000.0) << "s\n";
    // save upper triangular V_12
    v_dim[0] = v_sz;
    outfile_name = cpot + "V12_" + std::to_string(L) + ".h5";
    outfile = std::unique_ptr<H5::H5File>(new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    V_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(outfile->createDataSet(
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

int r_12::R12Glob4(std::string cpot, int L_max, int glq_pt, std::string dir) {
  bool min_dir=0, min_exc=0;
  double Y_norm=0.0;
  std::vector<double> v_mat;
  double sum_k=0.0;
  idx4 e12, e12p;

  int n=0;   // no. of points
  int bo=0;   // max B-spline order
  int nkn=0; // no. of knots
  std::vector<double> kkn;

  int L_sz=0, v_sz=0;
  std::string filename;
  std::string outfile_name;

  // HDF5 defines
  std::unique_ptr<H5::H5File> outfile=nullptr;
  std::unique_ptr<H5::DataSet> V_set=nullptr;
  std::unique_ptr<H5::DataSet> L_set=nullptr;
  std::vector<idx4> L_idx;
  H5::DataSpace cspace;
  hsize_t offset[2]={0,0},  stride[2]={1,1}, block[2]={1,1};
  hsize_t count_o[2], count_i[2], dimms_o[2], dimms_i[2], v_dim[1];
  H5::DataSpace memspace_o, memspace_i;

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

  int ncf, sym, max_N=0, l1_m=0, l2_m=0;
  std::vector<cfg::line> cfgs;
  cfg::line max_line;

  for(auto li=0; li<=L_max; ++li) {
    cfg::ReadCfg(dir, li, sym, ncf, cfgs);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
    max_N = std::max(max_N,max_line.n2max);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.l2 < b.l2; });
    l2_m = std::max(l2_m,max_line.l2);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });
    l1_m = std::max(l1_m,max_line.l1);
  }
  auto k_max = l1_m+l2_m;
  auto lc_sz = max_N*n*glq_pt;

  // reserve space for wfns
  auto e1_lm = std::max(l1_m,l2_m);
  std::vector<double> wfn_o(lc_sz*(e1_lm+1));
  std::vector<double> wfn_i(2*lc_sz*(e1_lm+1));

  count_o[0] = max_N;
  count_o[1] = n*glq_pt;
  dimms_o[0] = count_o[0]; 
  dimms_o[1] = count_o[1];
  memspace_o.setExtentSimple(2, dimms_o, NULL);
  count_i[0] = max_N;
  count_i[1] = 2*n*glq_pt;
  dimms_i[0] = count_i[0]; 
  dimms_i[1] = count_i[1];
  memspace_i.setExtentSimple(2, dimms_i, NULL);
  // read wavefunctions for all l
  for(int l=0; l<=e1_lm; ++l) {
    filename = cpot + "_w1e" + std::to_string(l) + ".h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Pr_o")));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count_o, offset, stride, block);
    rset->read(&wfn_o[l*lc_sz], H5::PredType::NATIVE_DOUBLE, memspace_o, cspace);
    rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Pr_i")));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count_i, offset, stride, block);
    rset->read(&wfn_i[l*lc_sz], H5::PredType::NATIVE_DOUBLE, memspace_i, cspace);
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


  // Precalculate r^k
  std::vector<double> rk, rk_mid;
  Rpowklob4(rk, rk_mid, n, bo, glq_pt, gl_x, kkn, k_max);

  omp_lock_t copylock;
  omp_init_lock(&copylock);
  uint64_t st_time;

  wig_table_init(2*(k_max+2),6);
  

  #pragma omp parallel private(e12p)
  {
  wig_thread_temp_init(2*(k_max+2));

  for(int L=0; L<=L_max; ++L) {
    
    // Read n1l1;n2l2 indices for NL states
    #pragma omp single
    {
    filename = cpot + std::to_string(L) + "idx.h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
    L_sz = L_set->getSpace().getSimpleExtentNpoints()/4;
    L_idx.reserve(L_sz);
    L_set->read(&L_idx[0], H5::PredType::NATIVE_INT32);
    file->close();

    std::cout << "L: " << L <<" L_sz: "<< L_sz << "\n";

    v_sz = L_sz*(L_sz+1)/2;
    v_mat.reserve(v_sz); // error

    st_time = GetTimeMs64();
    }

    for(int NL2=0; NL2<L_sz; ++NL2) {
      //set n1'l1';n2'l2'
      e12p =L_idx[NL2];
      bool eqvp = e12p.n1==e12p.n2 && e12p.l1==e12p.l2;

      #pragma omp barrier
      #pragma omp for private(e12,Y_norm,sum_k,min_dir,min_exc)
      for(int NL1=NL2; NL1<L_sz; ++NL1) {
        //set n1l1;n2l2
        e12 = L_idx[NL1];

        // sqrt([la][lc][lb][ld])
        Y_norm = sqrt((2*e12.l1+1)*(2*e12p.l1+1)*(2*e12.l2+1)*(2*e12p.l2+1));
        sum_k=0.0;
        min_dir = ((L+e12.l2+e12p.l1) >> 0) & 1; // check if L+lb+lc is odd
        min_exc = ((L+e12.l1+e12p.l1) >> 0) & 1;
        for(int k=0; k<=k_max; ++k) {
          /* for direct check if:
            (-)^{L+lb+lc}=(-)^{L+la+ld},
            |la-lc| <= k <= la+lc,
            |lb-ld| <= k <= lb+ld */
          if(min_dir ==(((L+e12.l1+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l1)<=k) && (k<=e12.l1+e12p.l1)) &&
             ((abs(e12.l2-e12p.l2)<=k) && (k<=e12.l2+e12p.l2)) &&
             !(((k+e12.l1+e12p.l1) >> 0) & 1) &&
             !(((k+e12.l2+e12p.l2) >> 0) & 1) &&
             !std::isinf(1.0/rk_mid[(k+1)*2*n*glq_pt])) {
            sum_k += pow(-1,min_dir)*FsltrLob4GL(k, n, bo, glq_pt, lc_sz,
                      e12.n1, e12.l1, e12.n2, e12.l2, e12p.n1, e12p.l1,
                      e12p.n2, e12p.l2, gl_w, gl_x, kkn, rk, rk_mid, wfn_o, wfn_i)
                      *wig3jj(2*e12.l1,2*k,2*e12p.l1,0,0,0)
                      *wig3jj(2*e12.l2,2*k,2*e12p.l2,0,0,0)
                      *wig6jj(2*e12p.l1,2*e12p.l2,2*L,2*e12.l2,2*e12.l1,2*k);
          }
          /* for exchange check if:
            (-)^{L+la+lc}=(-)^{L+lb+ld},
            |la-ld| <= k <= la+ld,
            |lc-lb| <= k <= lc+lb */
          if(min_exc ==(((L+e12.l2+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l2)<=k) && (k<=e12.l1+e12p.l2)) &&
             ((abs(e12.l2-e12p.l1)<=k) && (k<=e12.l2+e12p.l1)) &&
             !(((k+e12.l1+e12p.l2) >> 0) & 1) &&
             !(((k+e12.l2+e12p.l1) >> 0) & 1) &&
             !std::isinf(1.0/rk_mid[(k+1)*2*n*glq_pt])) {
            sum_k += pow(-1,min_exc)*FsltrLob4GL(k, n, bo, glq_pt, lc_sz,
                      e12.n1, e12.l1, e12.n2, e12.l2, e12p.n2, e12p.l2,
                      e12p.n1, e12p.l1, gl_w, gl_x, kkn, rk, rk_mid, wfn_o, wfn_i)
                      *wig3jj(2*e12.l1,2*k,2*e12p.l2,0,0,0)
                      *wig3jj(2*e12.l2,2*k,2*e12p.l1,0,0,0)
                      *wig6jj(2*e12p.l1,2*e12p.l2,2*L,2*e12.l1,2*e12.l2,2*k);
          }
        }

        // write symmetric V_12 as upper triangular
        double v12 = pow(-1,e12.l1+e12.l2)*Y_norm*sum_k;

        if(e12.n1==e12.n2 && e12.l1==e12.l2) 
          v12 *= 0.7071067811865475244008444e0;

        if(eqvp) 
          v12 *= 0.7071067811865475244008444e0;

        v_mat[(2*L_sz-NL2-1)*NL2/2 + NL1] = v12;
      }
    }

    #pragma omp single
    {
    std::cout << "loop time: " << ((double)(GetTimeMs64()-st_time)/1000.0) << "s\n";
    // save upper triangular V_12
    v_dim[0] = v_sz;
    outfile_name = cpot + "V12_" + std::to_string(L) + ".h5";
    outfile = std::unique_ptr<H5::H5File>(new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    V_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(outfile->createDataSet(
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

int r_12::R12Glob3(std::string cpot, int L_max, int glq_pt, std::string dir) {
  bool min_dir=0, min_exc=0;
  double Y_norm=0.0;
  std::vector<double> v_mat;
  double sum_k=0.0;
  idx4 e12, e12p;

  int n=0;   // no. of points
  int bo=0;   // max B-spline order
  int nkn=0; // no. of knots
  std::vector<double> kkn;

  int L_sz=0, v_sz=0;
  std::string filename;
  std::string outfile_name;

  // HDF5 defines
  std::unique_ptr<H5::H5File> outfile=nullptr;
  std::unique_ptr<H5::DataSet> V_set=nullptr;
  std::unique_ptr<H5::DataSet> L_set=nullptr;
  std::vector<idx4> L_idx;
  H5::DataSpace cspace;
  hsize_t offset[2]={0,0}, count[2], stride[2]={1,1}, block[2]={1,1};
  hsize_t dimms[2], v_dim[1];
  H5::DataSpace memspace;

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

  int ncf, sym, max_N=0, l1_m=0, l2_m=0;
  std::vector<cfg::line> cfgs;
  cfg::line max_line;

  for(auto li=0; li<=L_max; ++li) {
    cfg::ReadCfg(dir, li, sym, ncf, cfgs);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
    max_N = std::max(max_N,max_line.n2max);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.l2 < b.l2; });
    l2_m = std::max(l2_m,max_line.l2);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });
    l1_m = std::max(l1_m,max_line.l1);
  }
  auto k_max = l1_m+l2_m;
  auto lc_sz = max_N*n*glq_pt;

  // reserve space for coefficients
  auto e1_lm = std::max(l1_m,l2_m);
  std::vector<double> wfn_o(lc_sz*(e1_lm+1));
  std::vector<double> wfn_i(lc_sz*(e1_lm+1));

  count[0] = max_N;
  count[1] = n*glq_pt;
  dimms[0] = count[0]; 
  dimms[1] = count[1];
  memspace.setExtentSimple(2, dimms, NULL);
  // read wavefunctions for all l
  for(int l=0; l<=e1_lm; ++l) {
    filename = cpot + "_w1e" + std::to_string(l) + ".h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Pr_o")));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(&wfn_o[l*lc_sz], H5::PredType::NATIVE_DOUBLE, memspace, cspace);
    rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Pr_i")));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(&wfn_i[l*lc_sz], H5::PredType::NATIVE_DOUBLE, memspace, cspace);
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

  // Precalculate r^k
  std::vector<double> rk, rk_mid;
  Rpowk(rk, rk_mid, n, bo, glq_pt, gl_x, kkn, k_max);

  omp_lock_t copylock;
  omp_init_lock(&copylock);
  uint64_t st_time;

  wig_table_init(2*(k_max+2),6);
  

  #pragma omp parallel private(e12p)
  {
  wig_thread_temp_init(2*(k_max+2));

  for(int L=0; L<=L_max; ++L) {
    
    // Read n1l1;n2l2 indices for NL states
    #pragma omp single
    {
    filename = cpot + std::to_string(L) + "idx.h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
    L_sz = L_set->getSpace().getSimpleExtentNpoints()/4;
    L_idx.reserve(L_sz);
    L_set->read(&L_idx[0], H5::PredType::NATIVE_INT32);
    file->close();

    std::cout << "L: " << L <<" L_sz: "<< L_sz << "\n";

    v_sz = L_sz*(L_sz+1)/2;
    v_mat.reserve(v_sz); // error

    st_time = GetTimeMs64();
    }

    for(int NL2=0; NL2<L_sz; ++NL2) {
      //set n1'l1';n2'l2'
      e12p =L_idx[NL2];
      bool eqvp = e12p.n1==e12p.n2 && e12p.l1==e12p.l2;

      #pragma omp barrier
      #pragma omp for private(e12,Y_norm,sum_k,min_dir,min_exc)
      for(int NL1=NL2; NL1<L_sz; ++NL1) {
        //set n1l1;n2l2
        e12 = L_idx[NL1];

        // sqrt([la][lc][lb][ld])
        Y_norm = sqrt((2*e12.l1+1)*(2*e12p.l1+1)*(2*e12.l2+1)*(2*e12p.l2+1));
        sum_k=0.0;
        min_dir = ((L+e12.l2+e12p.l1) >> 0) & 1; // check if L+lb+lc is odd
        min_exc = ((L+e12.l1+e12p.l1) >> 0) & 1;
        for(int k=0; k<=k_max; ++k) {
          /* for direct check if:
            (-)^{L+lb+lc}=(-)^{L+la+ld},
            |la-lc| <= k <= la+lc,
            |lb-ld| <= k <= lb+ld */
          if(min_dir ==(((L+e12.l1+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l1)<=k) && (k<=e12.l1+e12p.l1)) &&
             ((abs(e12.l2-e12p.l2)<=k) && (k<=e12.l2+e12p.l2)) &&
             !(((k+e12.l1+e12p.l1) >> 0) & 1) &&
             !(((k+e12.l2+e12p.l2) >> 0) & 1) &&
             !std::isinf(1.0/rk_mid[(k+1)*n*glq_pt])) {
            sum_k += pow(-1,min_dir)*FsltrLob3GL(k, n, bo, glq_pt, lc_sz,
                      e12.n1, e12.l1, e12.n2, e12.l2, e12p.n1, e12p.l1,
                      e12p.n2, e12p.l2, gl_w, gl_x, kkn, rk, rk_mid, wfn_o, wfn_i)
                      *wig3jj(2*e12.l1,2*k,2*e12p.l1,0,0,0)
                      *wig3jj(2*e12.l2,2*k,2*e12p.l2,0,0,0)
                      *wig6jj(2*e12p.l1,2*e12p.l2,2*L,2*e12.l2,2*e12.l1,2*k);
          }
          /* for exchange check if:
            (-)^{L+la+lc}=(-)^{L+lb+ld},
            |la-ld| <= k <= la+ld,
            |lc-lb| <= k <= lc+lb */
          if(min_exc ==(((L+e12.l2+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l2)<=k) && (k<=e12.l1+e12p.l2)) &&
             ((abs(e12.l2-e12p.l1)<=k) && (k<=e12.l2+e12p.l1)) &&
             !(((k+e12.l1+e12p.l2) >> 0) & 1) &&
             !(((k+e12.l2+e12p.l1) >> 0) & 1) &&
             !std::isinf(1.0/rk_mid[(k+1)*n*glq_pt])) {
            sum_k += pow(-1,min_exc)*FsltrLob3GL(k, n, bo, glq_pt, lc_sz,
                      e12.n1, e12.l1, e12.n2, e12.l2, e12p.n2, e12p.l2,
                      e12p.n1, e12p.l1, gl_w, gl_x, kkn, rk, rk_mid, wfn_o, wfn_i)
                      *wig3jj(2*e12.l1,2*k,2*e12p.l2,0,0,0)
                      *wig3jj(2*e12.l2,2*k,2*e12p.l1,0,0,0)
                      *wig6jj(2*e12p.l1,2*e12p.l2,2*L,2*e12.l1,2*e12.l2,2*k);
          }
        }

        // write symmetric V_12 as upper triangular
        double v12 = pow(-1,e12.l1+e12.l2)*Y_norm*sum_k;

        if(e12.n1==e12.n2 && e12.l1==e12.l2) 
          v12 *= 0.7071067811865475244008444e0;

        if(eqvp) 
          v12 *= 0.7071067811865475244008444e0;

        v_mat[(2*L_sz-NL2-1)*NL2/2 + NL1] = v12;
      }
    }

    #pragma omp single
    {
    std::cout << "loop time: " << ((double)(GetTimeMs64()-st_time)/1000.0) << "s\n";
    // save upper triangular V_12
    v_dim[0] = v_sz;
    outfile_name = cpot + "V12_" + std::to_string(L) + ".h5";
    outfile = std::unique_ptr<H5::H5File>(new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    V_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(outfile->createDataSet(
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

int r_12::R12Trap(std::string cpot, int L_max, int glq_pt, std::string dir) {
  bool min_dir=0, min_exc=0;
  double Y_norm=0.0;
  std::vector<double> v_mat;
  double sum_k=0.0;
  idx4 e12, e12p;

  int n=0;   // no. of points
  int bo=0;   // max B-spline order
  int nkn=0; // no. of knots
  std::vector<double> kkn;

  int L_sz=0, v_sz=0;
  std::string filename;
  std::string outfile_name;

  // HDF5 defines
  std::unique_ptr<H5::H5File> outfile=nullptr;
  std::unique_ptr<H5::DataSet> V_set=nullptr;
  std::unique_ptr<H5::DataSet> L_set=nullptr;
  std::vector<idx4> L_idx;
  H5::DataSpace cspace;
  hsize_t offset[2]={0,0}, count[2], stride[2]={1,1}, block[2]={1,1};
  hsize_t dimms[2], v_dim[1];
  H5::DataSpace memspace;

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

  int ncf, sym, max_N=0, l1_m=0, l2_m=0;
  std::vector<cfg::line> cfgs;
  cfg::line max_line;

  for(auto li=0; li<=L_max; ++li) {
    cfg::ReadCfg(dir, li, sym, ncf, cfgs);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
    max_N = std::max(max_N,max_line.n2max);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.l2 < b.l2; });
    l2_m = std::max(l2_m,max_line.l2);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });
    l1_m = std::max(l1_m,max_line.l1);
  }
  auto k_max = l1_m+l2_m;
  auto lc_sz = max_N*n*glq_pt;

  // reserve space for coefficients
  auto e1_lm = std::max(l1_m,l2_m);
  std::vector<double> wfn_o(lc_sz*(e1_lm+1));

  count[0] = max_N;
  count[1] = n*glq_pt;
  dimms[0] = count[0]; 
  dimms[1] = count[1];
  memspace.setExtentSimple(2, dimms, NULL);
  // read wavefunctions for all l
  for(int l=0; l<=e1_lm; ++l) {
    filename = cpot + "_w1e" + std::to_string(l) + ".h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Pr_o")));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(&wfn_o[l*lc_sz], H5::PredType::NATIVE_DOUBLE, memspace, cspace);
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

  // Precalculate r^k
  std::vector<double> rk;//, rk_mid;
  Rpowk(rk, n, bo, glq_pt, gl_x, kkn, k_max);

  omp_lock_t copylock;
  omp_init_lock(&copylock);
  uint64_t st_time;

  wig_table_init(2*(k_max+2),6);
  

  #pragma omp parallel private(e12p)
  {
  wig_thread_temp_init(2*(k_max+2));

  for(int L=0; L<=L_max; ++L) {
    
    // Read n1l1;n2l2 indices for NL states
    #pragma omp single
    {
    filename = cpot + std::to_string(L) + "idx.h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
    L_sz = L_set->getSpace().getSimpleExtentNpoints()/4;
    L_idx.reserve(L_sz);
    L_set->read(&L_idx[0], H5::PredType::NATIVE_INT32);
    file->close();

    std::cout << "L: " << L <<" L_sz: "<< L_sz << "\n";

    v_sz = L_sz*(L_sz+1)/2;
    v_mat.reserve(v_sz); // error

    st_time = GetTimeMs64();
    }

    for(int NL2=0; NL2<L_sz; ++NL2) {
      //set n1'l1';n2'l2'
      e12p =L_idx[NL2];
      bool eqvp = e12p.n1==e12p.n2 && e12p.l1==e12p.l2;

      #pragma omp barrier
      #pragma omp for private(e12,Y_norm,sum_k,min_dir,min_exc)
      for(int NL1=NL2; NL1<L_sz; ++NL1) {
        //set n1l1;n2l2
        e12 = L_idx[NL1];

        // sqrt([la][lc][lb][ld])
        Y_norm = sqrt((2*e12.l1+1)*(2*e12p.l1+1)*(2*e12.l2+1)*(2*e12p.l2+1));
        sum_k=0.0;
        min_dir = ((L+e12.l2+e12p.l1) >> 0) & 1; // check if L+lb+lc is odd
        min_exc = ((L+e12.l1+e12p.l1) >> 0) & 1;
        for(int k=0; k<=k_max; ++k) {
          /* for direct check if:
            (-)^{L+lb+lc}=(-)^{L+la+ld},
            |la-lc| <= k <= la+lc,
            |lb-ld| <= k <= lb+ld */
          if(min_dir ==(((L+e12.l1+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l1)<=k) && (k<=e12.l1+e12p.l1)) &&
             ((abs(e12.l2-e12p.l2)<=k) && (k<=e12.l2+e12p.l2)) &&
             !(((k+e12.l1+e12p.l1) >> 0) & 1) &&
             !(((k+e12.l2+e12p.l2) >> 0) & 1) &&
             !std::isinf(1.0/rk[(k+1)*n*glq_pt])) {
            sum_k += pow(-1,min_dir)*FsltrTrapGL(k, n, bo, glq_pt, lc_sz,
                      e12.n1, e12.l1, e12.n2, e12.l2, e12p.n1, e12p.l1,
                      e12p.n2, e12p.l2, gl_w, gl_x, kkn, rk, wfn_o)
                      *wig3jj(2*e12.l1,2*k,2*e12p.l1,0,0,0)
                      *wig3jj(2*e12.l2,2*k,2*e12p.l2,0,0,0)
                      *wig6jj(2*e12p.l1,2*e12p.l2,2*L,2*e12.l2,2*e12.l1,2*k);
          }
          /* for exchange check if:
            (-)^{L+la+lc}=(-)^{L+lb+ld},
            |la-ld| <= k <= la+ld,
            |lc-lb| <= k <= lc+lb */
          if(min_exc ==(((L+e12.l2+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l2)<=k) && (k<=e12.l1+e12p.l2)) &&
             ((abs(e12.l2-e12p.l1)<=k) && (k<=e12.l2+e12p.l1)) &&
             !(((k+e12.l1+e12p.l2) >> 0) & 1) &&
             !(((k+e12.l2+e12p.l1) >> 0) & 1) &&
             !std::isinf(1.0/rk[(k+1)*n*glq_pt])) {
            sum_k += pow(-1,min_exc)*FsltrTrapGL(k, n, bo, glq_pt, lc_sz,
                      e12.n1, e12.l1, e12.n2, e12.l2, e12p.n2, e12p.l2,
                      e12p.n1, e12p.l1, gl_w, gl_x, kkn, rk, wfn_o)
                      *wig3jj(2*e12.l1,2*k,2*e12p.l2,0,0,0)
                      *wig3jj(2*e12.l2,2*k,2*e12p.l1,0,0,0)
                      *wig6jj(2*e12p.l1,2*e12p.l2,2*L,2*e12.l1,2*e12.l2,2*k);
          }
        }

        // write symmetric V_12 as upper triangular
        double v12 = pow(-1,e12.l1+e12.l2)*Y_norm*sum_k;

        if(e12.n1==e12.n2 && e12.l1==e12.l2) 
          v12 *= 0.7071067811865475244008444e0;

        if(eqvp) 
          v12 *= 0.7071067811865475244008444e0;

        v_mat[(2*L_sz-NL2-1)*NL2/2 + NL1] = v12;
      }
    }

    #pragma omp single
    {
    std::cout << "loop time: " << ((double)(GetTimeMs64()-st_time)/1000.0) << "s\n";
    // save upper triangular V_12
    v_dim[0] = v_sz;
    outfile_name = cpot + "V12_" + std::to_string(L) + ".h5";
    outfile = std::unique_ptr<H5::H5File>(new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    V_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(outfile->createDataSet(
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
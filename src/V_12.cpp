#include "V_12.h"

double Fsltr(int k, int n, int bo,
            int na, int la, int nb, int lb,
            int nc, int lc, int nd, int ld,
            std::vector<uint> &N_max,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<int> &offset,
            std::vector<double> &Cf) {
  int kp1=k+1;
  double *Cl1i_pt=&Cf[offset[la]];
  double *Cl1p_pt=&Cf[offset[lc]];
  double *Cl2i_pt=&Cf[offset[lb]];
  double *Cl2p_pt=&Cf[offset[ld]];
  double Pl1i=0, Pl1p=0, Pl2i=0, Pl2p=0;
  double dl, sl, loc_GL, r2, r1, pr2, chi;

  // first calculate Qk for all of 0->R
  double Qk=0;
  double Jk=0;
  double Fk=0;
  int bidx=0;
  for(auto i=bo-1; i<n; ++i, ++bidx) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(int p=0; p<bo; ++p){
      r2 = dl*gl_x[p] + sl;
      Pl2i = 0; Pl2p = 0;

      for(int j=0; j<bo; ++j) {
        Pl2i += Cl2i_pt[nb*n+i-bo+1+j]*Bsp[j+bo*(p+bidx*bo)];
        Pl2p += Cl2p_pt[nd*n+i-bo+1+j]*Bsp[j+bo*(p+bidx*bo)];
      }
      loc_GL+=gl_w[p]*pow(r2,-kp1)*Pl2i*Pl2p;
    }
    Qk+=dl*loc_GL;
  }

  bidx=0;
  for(auto i=bo-1; i<n; ++i, ++bidx) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(int p=0; p<bo; ++p){ // need to work around this index for chi calculation
      r1 = dl*gl_x[p] + sl;
      Pl1i = 0; Pl1p = 0;
      Pl2i = 0; Pl2p = 0;


      for(int j=0; j<bo; ++j) {
        Pl1i += Cl1i_pt[na*n+i-bo+1+j]*Bsp[j+bo*(p+bidx*bo)];
        Pl1p += Cl1p_pt[nc*n+i-bo+1+j]*Bsp[j+bo*(p+bidx*bo)];
        Pl2i += Cl2i_pt[nb*n+i-bo+1+j]*Bsp[j+bo*(p+bidx*bo)];
        Pl2p += Cl2p_pt[nd*n+i-bo+1+j]*Bsp[j+bo*(p+bidx*bo)];
      }
      pr2 = Pl2i*Pl2p;
      // chi(r1)
      Jk+=dl*gl_w[p]*pow(r1,k)*pr2;
      Qk-=dl*gl_w[p]*pow(r1,-kp1)*pr2;
      chi=pow(r1,-kp1)*Jk+pow(r1,k)*Qk;

      // Fk12;1'2'
      loc_GL += gl_w[p]*Pl1i*Pl1p*chi;
      // if (i==bo-1) {std::cout << Jk << " " << r1 << " " << gl_w[p]*dl
      // << " " << Qk << " " << dl*loc_GL << "\n";}
      // if(nc==la&&la==lb&&lc==ld&&la==0&&k==0) {
      //   std::cout << std::setiosflags(std::ios::scientific)
      //         << std::setprecision(12)<< r1 << " " << Pl1i << "\n";
      // }
    }
    Fk+=dl*loc_GL;
  }
  std::cout << Fk << '\n';
  return Fk;
}

double Fsltr_alt(int k, int n, int bo,
            int na, int la, int nb, int lb,
            int nc, int lc, int nd, int ld,
            std::vector<uint> &N_max,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<int> &offset,
            std::vector<double> &Cf) {
  int kp1=k+1;
  double *Cl1i_pt=&Cf[offset[la]];
  double *Cl1p_pt=&Cf[offset[lc]];
  double *Cl2i_pt=&Cf[offset[lb]];
  double *Cl2p_pt=&Cf[offset[ld]];
  double Pl1i=0, Pl1p=0, Pl2i=0, Pl2p=0;
  double dl, sl, dl2, sl2, loc_GL, r2, r1, pr2, chi;

  // first calculate Qk for all of 0->R
  double Qk=0;
  double Jk=0;
  double Fk=0;
  int bidx=0; int bidx2=0;
  for(auto i=bo-1; i<n; ++i, ++bidx) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(int p=0; p<bo; ++p){ // need to work around this index for chi calculation
      r1 = dl*gl_x[p] + sl;
      Pl1i = 0; Pl1p = 0;

      for(int j=0; j<bo; ++j) {
        Pl1i += Cl1i_pt[na*n+i-bo+1+j]*Bsp[j+bo*(p+bidx*bo)];
        Pl1p += Cl1p_pt[nc*n+i-bo+1+j]*Bsp[j+bo*(p+bidx*bo)];
      }

      bidx2=0;
      Jk=0; Qk=0;
      for(auto i2=bo-1; i2<n; ++i2, ++bidx2) {
      dl2 = (kkn[i2+1] - kkn[i2])*0.5;
      sl2 = (kkn[i2+1] + kkn[i2])*0.5;
        for(int p2=0; p2<bo; ++p2){
          r2 = dl*gl_x[p2] + sl;
          Pl2i = 0; Pl2p = 0;

          for(int j2=0; j2<bo; ++j2) {
            Pl2i += Cl2i_pt[nb*n+i2-bo+1+j2]*Bsp[j2+bo*(p2+bidx2*bo)];
            Pl2p += Cl2p_pt[nd*n+i2-bo+1+j2]*Bsp[j2+bo*(p2+bidx2*bo)];
          }
          pr2 = Pl2i*Pl2p;
          if(r2<=r1) {
            Jk+=dl2*gl_w[p2]*pow(r2/r1,k)*pr2;
          } else if (r2>r1) {
            Qk+=dl2*gl_w[p2]*pow(r1/r2,kp1)*pr2;
          }
        }
      }
      // chi(r1)
      chi=Jk+Qk;

      // Fk12;1'2'
      Fk += dl*gl_w[p]*Pl1i*Pl1p*chi/r1;
      if (i==bo-1) {std::cout << r1 << " " << Fk << "\n";}

    }
    //Fk+=dl*loc_GL;
  }

  return Fk;
}

int V12(std::string cpot, uint L_max, std::vector<uint> &N_sz) {
  uint min_dir=0, min_exc=0;
  double Y_norm=0.0;
  std::vector<double> v_mat; // standin for V_12 matrix
  double sum_k=0.0;
  idx4 e12, e12p;

  int n=0;   // no. of points
  int bo=0;   // max B-spline order
  int nkn=0; // no. of knots
  int nSt=0; // no. of states
  std::vector<double> kkn;
  std::vector<double> Cf;
  std::vector<double*> C(L_max+1);

  uint L_sz=0, v_sz=0;
  std::string filename;
  std::string outfile_name;
  H5::H5File *file=nullptr;
  H5::H5File *outfile=nullptr;
  std::vector<idx4> L_idx;
  H5::DataSet *L_set=nullptr;
  H5::DataSet *V_set=nullptr;
  H5::DataSet *rset=nullptr;
  H5::DataSpace cspace;
  hsize_t offset[2], count[2], stride[2], block[2];
  hsize_t dimms[2], v_dim[1];
  offset[0]=0; offset[1]=0;
  count[0] =0;
  stride[0]=1; stride[1]=1;
  block[0] =1; block[1]=1;
  dimms[0] =0;
  H5::DataSpace memspace;

  int tot_states = 0;
  for(auto &b : N_sz) {
    tot_states += b;
  }

  filename = cpot + std::to_string(0) + ".h5";
  file = new H5::H5File(filename, H5F_ACC_RDONLY);
  file->openAttribute("N").read(H5::PredType::NATIVE_INT32, &n);
  file->openAttribute("K").read(H5::PredType::NATIVE_INT32, &bo);
  rset = new H5::DataSet(file->openDataSet("Knots"));
  nkn = rset->getSpace().getSimpleExtentNpoints();

  kkn.reserve(nkn);
  rset->read(&kkn[0], H5::PredType::NATIVE_DOUBLE); //read knots
  delete rset;

  rset = new H5::DataSet(file->openDataSet("En"));
  nSt = rset->getSpace().getSimpleExtentNpoints();
  delete rset;
  delete file;

  Cf.reserve(tot_states*n);
  C[0]=&Cf[1];
  Cf[0]=0.0; Cf[nSt+1] =0.0;
  int nt=0;
  std::vector<int> nst_prev(L_max+1);
  nst_prev[0] = nt;
  for(int i=1; i<=L_max; ++i) {
    nt += N_sz[i-1];
    nst_prev[i] = nt*n;
    C[i] = &Cf[nt*n+1];
    Cf[nt*n]=0.0;
    Cf[i*n-1]=0.0;
  }

  // read coefficients for all l
  for(int l=0; l<=L_max; ++l) {
    count[0] = N_sz[l]; count[1] = nSt;
    dimms[0] = count[0]; dimms[1] = count[1];
    memspace.setExtentSimple(2, dimms, NULL);

    filename = cpot + std::to_string(l) + ".h5";
    file = new H5::H5File(filename, H5F_ACC_RDONLY);
    rset = new H5::DataSet(file->openDataSet("Coeff"));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(C[l], H5::PredType::NATIVE_DOUBLE, memspace, cspace);
    delete rset;
    delete file;
  }

  std::vector<double> gl_x(bo);
  std::vector<double> gl_w(bo);
  fastgl::QuadPair gl_i;
  for(int i=1; i<=bo; ++i) {
    gl_i = fastgl::GLPair(bo, i); // generate GL nodes and weights over B-splines support
    gl_x[bo-i] = gl_i.x(); 
    gl_w[bo-i] = gl_i.weight;
  }

  // generate B-splines
  std::vector<double> Bsplines;
  Bsplines.reserve(nSt*bo*bo);
  bsp::Splines(n, bo, gl_x, kkn, Bsplines);

  int L_real_size=0;
  filename = cpot + std::to_string(0) + "idx.h5";
  file = new H5::H5File(filename, H5F_ACC_RDONLY);
  L_set = new H5::DataSet(file->openDataSet("idx"));
  L_sz = 1;
  L_real_size = L_set->getSpace().getSimpleExtentNpoints()/4;
  L_idx.resize(L_real_size);
  delete L_set;
  delete file;

  for(uint L=0; L<=L_max; ++L) {
    std::cout << "L: " << L << "\n";
    // Read indices n1l1;n2l2
    filename = cpot + std::to_string(L) + "idx.h5";
    file = new H5::H5File(filename, H5F_ACC_RDONLY);
    L_set = new H5::DataSet(file->openDataSet("idx"));
    L_real_size = L_set->getSpace().getSimpleExtentNpoints()/4;
    L_idx.resize(L_real_size);
    L_set->read(&L_idx[0], H5::PredType::NATIVE_UINT32);
    delete L_set;
    delete file;

    v_sz = L_sz*(L_sz+1)/2;
    v_mat.reserve(v_sz);

    for(uint NL2=0; NL2<L_sz; ++NL2) {
      //set n1'l1';n2'l2'
      e12p = L_idx[NL2];
      for(uint NL1=NL2; NL1<L_sz; ++NL1){
        //set n1l1;n2l2
        e12 = L_idx[NL1];

        Y_norm = sqrt((2*e12.l1+1)*(2*e12p.l1+1)*(2*e12.l2+1)*(2*e12p.l2+1));
        sum_k=0.0;
        for(uint k=0; k<=L_max; ++k) { // should include l_max!=L_max ?
          min_dir = ((L+e12.l2+e12p.l1) >> 0) & 1;
          if(min_dir ==(((L+e12.l1+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l1)<=k) && (k<=e12.l1+e12p.l1)) &&
             ((abs(e12.l2-e12p.l2)<=k) && (k<=e12.l2+e12p.l2))) {
            sum_k += pow(-1,min_dir)*Fsltr(k, n, bo, e12.n1, e12.l1, e12.n2, e12.l2,
                                          e12p.n1, e12p.l1, e12p.n2, e12p.l2,
                                          N_sz, gl_w, gl_x, kkn, Bsplines, nst_prev, Cf)
                  *wigner_3j0(e12.l1,k,e12p.l1)*wigner_3j0(e12.l2,k,e12p.l2)
                  *wigner_6j(e12p.l1,e12p.l2,L,e12.l2,e12.l1,k);
          }
          min_exc = ((L+e12.l1+e12p.l1) >> 0) & 1;
          if(min_exc ==(((L+e12.l2+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l2)<=k) && (k<=e12.l1+e12p.l2)) &&
             ((abs(e12.l2-e12p.l1)<=k) && (k<=e12.l2+e12p.l1))) {
            sum_k += pow(-1,min_exc)*Fsltr(k, n, bo, e12.n1, e12.l1, e12.n2, e12.l2,
                                          e12p.n2, e12p.l2, e12p.n1, e12p.l1,
                                          N_sz, gl_w, gl_x, kkn, Bsplines, nst_prev, Cf)
                  *wigner_3j0(e12.l1,k,e12p.l2)*wigner_3j0(e12.l2,k,e12p.l1)
                  *wigner_6j(e12p.l1,e12p.l2,L,e12.l1,e12.l2,k);
          }
        }
        // write symmetric V_12 as upper triangular
        v_mat[(2*L_sz-NL2-1)*NL2/2 + NL1] = pow(-1,(e12.l1+e12.l2))*Y_norm*sum_k;
      } 
    }
    // save upper triangular V_12
    v_dim[0] = v_sz;
    outfile_name = cpot + "V12_" + std::to_string(L) + ".h5";
    outfile = new H5::H5File(outfile_name, H5F_ACC_TRUNC);
    V_set = new H5::DataSet(outfile->createDataSet("V_12", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, v_dim)));
    V_set->write(&v_mat[0], H5::PredType::NATIVE_DOUBLE);
    delete V_set;
    delete outfile;

    v_mat.clear();
  }

  return 0;
}
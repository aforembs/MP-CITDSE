#include "V_12.h"

// Fast slater integral code
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
  double Qk=0.0;
  double Jk=0;
  double Fk=0;

  std::ofstream outQk("dat/Qk"+std::to_string(nb)+std::to_string(nd)+".dat");

  for(auto i=bo-1; i<n; ++i) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<bo; ++p){
      r2 = dl*gl_x[p] + sl;
      Pl2i = 0; Pl2p = 0;

      for(auto j=0; j<bo; ++j) {
        Pl2i += Cl2i_pt[nb*n+i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
        Pl2p += Cl2p_pt[nd*n+i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
      }
      outQk << r2 << " " << Pl2i << " " << Pl2p << " "<< Pl2i*Pl2p <<" "<< gl_w[p]<<" "<<pow(r2,-kp1)<< "\n";
      loc_GL+=gl_w[p]*pow(r2,-kp1)*Pl2i*Pl2p;
    }
    // if (nb==1 && nd==0 && lb==0 && ld==0) {
    // std::cout <<"dl: "<< dl<< " loc_gl: "<< loc_GL << " Qk: " << Qk <<"\n";}
    Qk+=dl*loc_GL;
  }

  if (na==0&&nb==0&&nc==0&&nd==0&&lb==0&&ld==0) std::cout << std::setiosflags(std::ios::scientific)
        << std::setprecision(15) << "Qk: " << Qk <<"\n";

  std::ofstream outFile("dat/slt_test"+std::to_string(na)+std::to_string(la)+std::to_string(nb)+std::to_string(lb)
  +std::to_string(nc)+std::to_string(lc)+std::to_string(nd)+std::to_string(ld)+".dat", std::ofstream::out);

  // GL-quadrature won't work here!!!
  // might need to generate different points for Simpson's rule
  for(auto i=bo-1; i<n; ++i) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<bo; ++p){ // need to work around this index for chi calculation
      r1 = dl*gl_x[p] + sl;
      Pl1i = 0; Pl1p = 0;
      Pl2i = 0; Pl2p = 0;

      for(auto j=0; j<bo; ++j) {
        Pl1i += Cl1i_pt[na*n+i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
        Pl1p += Cl1p_pt[nc*n+i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
        Pl2i += Cl2i_pt[nb*n+i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
        Pl2p += Cl2p_pt[nd*n+i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
      }
      pr2 = Pl2i*Pl2p;
      // chi(r1)
      Jk+=dl*gl_w[p]*pow(r1,k)*pr2;
      Qk-=dl*gl_w[p]*pow(r1,-kp1)*pr2;
      chi=pow(r1,-kp1)*Jk+pow(r1,k)*Qk;

        outFile <<r1<<" "<<Jk<<" "<<Qk<<" "<<pr2<<" "<<chi<< "\n";

      // Fk12;1'2'
      loc_GL += gl_w[p]*Pl1i*Pl1p*chi;
    }
    Fk+=dl*loc_GL;
  }
  outFile.close();

  return Fk;
}

// Fast slater integral code with simpson's rule inner int
double FsltrSimp(int k, int n, int bo,
            int na, int la, int nb, int lb,
            int nc, int lc, int nd, int ld,
            std::vector<uint> &N_max,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Ssp,
            std::vector<int> &offset,
            std::vector<double> &Cf) {
  int kp1=k+1;
  double *Cl1i_pt=&Cf[offset[la]];
  double *Cl1p_pt=&Cf[offset[lc]];
  double *Cl2i_pt=&Cf[offset[lb]];
  double *Cl2p_pt=&Cf[offset[ld]];
  double Pl1i=0, Pl1p=0, Pl2i=0, Pl2p=0;
  double Pl2ira=0, Pl2pra=0, Pl2irb=0, Pl2prb=0;
  double dl, sl, loc_GL, r2, r1, pr2, pr2a, pr2b, chi;

  // first calculate Qk for all of 0->R
  double Qk=0.0;
  double Jk=0;
  double Fk=0;

  for(auto i=bo-1; i<n; ++i) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<bo; ++p){
      r2 = dl*gl_x[p] + sl;
      Pl2i = 0; Pl2p = 0;

      for(auto j=0; j<bo; ++j) {
        Pl2i += Cl2i_pt[nb*n+i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
        Pl2p += Cl2p_pt[nd*n+i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
      }
      loc_GL+=gl_w[p]*pow(r2,-kp1)*Pl2i*Pl2p;
    }
    Qk+=dl*loc_GL;
  }

  if (na==0&&nb==0&&nc==0&&nd==0&&lb==0&&ld==0) std::cout << std::setiosflags(std::ios::scientific)
        << std::setprecision(15) << "Qk: " << Qk <<"\n";

  std::ofstream outFile("dat/slt_test"+std::to_string(na)+std::to_string(la)+std::to_string(nb)+std::to_string(lb)
  +std::to_string(nc)+std::to_string(lc)+std::to_string(nd)+std::to_string(ld)+".dat", std::ofstream::out);

  double rm1=0.0;
  double rs38a, rs38b;
  double fm1j=0.0, fm1q=0.0;
  auto sp=0;
  double fact, pr1k, pr1km;
  auto ibo1j=0;
  auto ai=0;
  auto bi=0;
  // GL-quadrature won't work here!!!
  // might need to generate different points for Simpson's rule
  for(auto i=bo-1; i<n; ++i) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    sp=0;
    for(auto p=0; p<bo; ++p, sp+=2){ // need to work around this index for chi calculation
      r1 = dl*gl_x[p] + sl;
      rs38a = (2*rm1+r1)/3;
      rs38b = (rm1+2*r1)/3;
      Pl1i = 0; Pl1p = 0;
      Pl2i = 0; Pl2p = 0;
      Pl2ira=0; Pl2pra=0; 
      Pl2irb=0; Pl2prb=0;

      for(auto j=0; j<bo; ++j) {
        ibo1j = i-bo+1+j;
        Pl1i += Cl1i_pt[na*n+ibo1j]*Bsp[j+bo*(p+i*bo)];
        Pl1p += Cl1p_pt[nc*n+ibo1j]*Bsp[j+bo*(p+i*bo)];
        Pl2i += Cl2i_pt[nb*n+ibo1j]*Bsp[j+bo*(p+i*bo)];
        Pl2p += Cl2p_pt[nd*n+ibo1j]*Bsp[j+bo*(p+i*bo)];

        ai=ibo1j-(kkn[i]>rs38a);
        bi=ibo1j-(kkn[i]>rs38b);
        Pl2ira += Cl2i_pt[nb*n+ai]*Ssp[j+bo*(sp+i*2*bo)];
        Pl2pra += Cl2p_pt[nd*n+ai]*Ssp[j+bo*(sp+i*2*bo)];

        Pl2irb += Cl2i_pt[nb*n+bi]*Ssp[j+bo*(sp+1+i*2*bo)];
        Pl2prb += Cl2p_pt[nd*n+bi]*Ssp[j+bo*(sp+1+i*2*bo)];
      }
      pr2 = Pl2i*Pl2p;
      pr2a= Pl2ira*Pl2pra;
      pr2b= Pl2irb*Pl2prb;
      // chi(r1)
      fact = (r1-rm1)*0.125;
      pr1k = pow(r1,k);
      pr1km= pow(r1,-kp1);
      Jk+=fact*(fm1j+3*pow(rs38a,k)*pr2a+3*pow(rs38b,k)*pr2b+pr1k*pr2);
      Qk-=fact*(fm1q+3*pow(rs38a,-kp1)*pr2a+3*pow(rs38b,-kp1)*pr2b+pr1km*pr2);
      chi=pr1km*Jk+pr1k*Qk;

        outFile <<r1<<" "<<rs38a<<" "<<rs38b<<" "<<Pl2i<<" "<<Pl2ira<<" "<<Pl2irb<<" "<<chi<< "\n";

      // Fk12;1'2'
      loc_GL += gl_w[p]*Pl1i*Pl1p*chi;

      rm1=r1;
      fm1j=pr1k*pr2;
      fm1q=pr1km*pr2;
    }
    Fk+=dl*loc_GL;
  }
  outFile.close();

  return Fk;
}

// alternative slater integral code
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
  double dl, sl, sl2, dl2, r2, r1, pr2, chi;

  // first calculate Qk for all of 0->R
  double Qk=0;
  double Jk=0;
  double Fk=0;
  for(auto i=bo-1; i<n; ++i) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;

    for(int p=0; p<bo; ++p){ // need to work around this index for chi calculation
      r1 = dl*gl_x[p] + sl;
      Pl1i = 0; Pl1p = 0;

      for(int j=0; j<bo; ++j) {
        Pl1i += Cl1i_pt[na*n+i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
        Pl1p += Cl1p_pt[nc*n+i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
      }

      Jk=0; Qk=0;
      for(auto i2=bo-1; i2<n; ++i2) {
      dl2 = (kkn[i2+1] - kkn[i2])*0.5;
      sl2 = (kkn[i2+1] + kkn[i2])*0.5;
        for(int p2=0; p2<bo; ++p2){
          r2 = dl2*gl_x[p2] + sl2;
          Pl2i = 0; Pl2p = 0;

          for(int j2=0; j2<bo; ++j2) {
            Pl2i += Cl2i_pt[nb*n+i2-bo+1+j2]*Bsp[j2+bo*(p2+i2*bo)];
            Pl2p += Cl2p_pt[nd*n+i2-bo+1+j2]*Bsp[j2+bo*(p2+i2*bo)];
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
      //if (i==bo-1) {std::cout << r1 << " " << Fk << "\n";}

    }
    //Fk+=dl*loc_GL;
  }

  return Fk;
}

int V12(std::string cpot, uint L_max, std::vector<uint> &N_sz) {
  uint min_dir=0, min_exc=0;
  double Y_norm=0.0;
  std::vector<double> v_mat;
  double sum_k=0.0;
  idx4 e12, e12p;

  int n=0;   // no. of points
  int bo=0;   // max B-spline order
  int nkn=0; // no. of knots
  std::vector<double> kkn;
  std::vector<double> Cf;
  std::vector<double*> C(L_max+1);

  uint L_sz=0, v_sz=0;
  std::string filename;
  std::string outfile_name;

  // HDF5 defines
  std::unique_ptr<H5::H5File> outfile=nullptr;
  std::unique_ptr<H5::DataSet> V_set=nullptr;
  std::unique_ptr<H5::DataSet> L_set=nullptr;
  std::vector<idx4> L_idx;
  H5::DataSpace cspace;
  hsize_t offset[2], count[2], stride[2], block[2];
  hsize_t dimms[2], v_dim[1];
  offset[0]=0; 
  offset[1]=0;
  count[0] =0;
  stride[0]=1; 
  stride[1]=1;
  block[0] =1; 
  block[1]=1;
  H5::DataSpace memspace;

  int tot_states = 0;
  for(auto &b : N_sz) {
    tot_states += b;
  }

  // read knots
  filename = cpot + std::to_string(0) + ".h5";
  auto file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  file->openAttribute("N").read(H5::PredType::NATIVE_INT32, &n);
  file->openAttribute("K").read(H5::PredType::NATIVE_INT32, &bo);
  auto rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Knots")));
  nkn = rset->getSpace().getSimpleExtentNpoints();

  kkn.reserve(nkn);
  rset->read(&kkn[0], H5::PredType::NATIVE_DOUBLE); //read knots

  // reserve space for coefficients
  Cf.reserve(tot_states*n);
  C[0]=&Cf[0];
  int nt=0;
  std::vector<int> nst_prev(L_max+1);
  nst_prev[0] = nt;
  for(uint i=1; i<=L_max; ++i) {
    nt += N_sz[i-1];
    nst_prev[i] = nt*n;
    C[i] = &Cf[nt*n];
  }

  // read coefficients for all l
  for(uint l=0; l<=L_max; ++l) {
    count[0] = N_sz[l]; 
    count[1] = n;
    dimms[0] = count[0]; 
    dimms[1] = count[1];
    memspace.setExtentSimple(2, dimms, NULL);

    filename = cpot + std::to_string(l) + ".h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Coeff")));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(C[l], H5::PredType::NATIVE_DOUBLE, memspace, cspace);
  }

  // generate GL nodes and weights over B-splines support
  std::vector<double> gl_x(bo);
  std::vector<double> gl_w(bo);
  fastgl::QuadPair gl_i;
  for(int i=1; i<=bo; ++i) {
    gl_i = fastgl::GLPair(bo, i);
    gl_x[bo-i] = gl_i.x(); 
    gl_w[bo-i] = gl_i.weight;
  }

  // generate B-splines
  std::vector<double> Bsplines;
  bsp::Splines(n, bo, gl_x, kkn, Bsplines);

  std::vector<double> Ssp;
  bsp::SimpSplines(n, bo, gl_x, kkn, Ssp);

  // read size of n1l1;n2l2 index vector
  int L_real_size=0;
  // filename = cpot + std::to_string(0) + "idx.h5";
  // file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  // auto L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
  
  // L_real_size = L_set->getSpace().getSimpleExtentNpoints()/4;
  // L_idx.resize(L_real_size);
  L_sz = 2; // read only N=1 and N=2 for each L

  for(uint L=0; L<=L_max; ++L) {
    std::cout << "L: " << L << "\n";
    // Read n1l1;n2l2 indices for NL states
    filename = cpot + std::to_string(L) + "idx.h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
    L_real_size = L_set->getSpace().getSimpleExtentNpoints()/4;
    L_idx.resize(L_real_size);
    L_set->read(&L_idx[0], H5::PredType::NATIVE_UINT32);

    v_sz = L_sz*(L_sz+1)/2;
    v_mat.reserve(v_sz);

    for(uint NL2=0; NL2<L_sz; ++NL2) {
      //set n1'l1';n2'l2'
      e12p = L_idx[NL2];
      std::cout << e12p.n1<<" "<<e12p.l1<<" "<<e12p.n2<<" "<<e12p.l2 <<"\n";
      for(uint NL1=NL2; NL1<L_sz; ++NL1){
        //set n1l1;n2l2
        e12 = L_idx[NL1];
        std::cout <<e12.n1<<" "<<e12.l1<<" "<<e12.n2<<" "<<e12.l2 <<"\n";
        // sqrt([la][lc][lb][ld])
        Y_norm = sqrt((2*e12.l1+1)*(2*e12p.l1+1)*(2*e12.l2+1)*(2*e12p.l2+1));
        sum_k=0.0;
        for(uint k=0; k<=L_max+1; ++k) { 
          /* for direct check if:
            (-)^{L+lb+lc}=(-)^{L+la+ld},
            |la-lc| <= k <= la+lc,
            |lb-ld| <= k <= lb+ld */
          min_dir = ((L+e12.l2+e12p.l1) >> 0) & 1; // check if L+lb+lc is even
          if(min_dir ==(((L+e12.l1+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l1)<=k) && (k<=e12.l1+e12p.l1)) &&
             ((abs(e12.l2-e12p.l2)<=k) && (k<=e12.l2+e12p.l2))) {
            sum_k += pow(-1,min_dir)*FsltrSimp(k, n, bo, e12.n1, e12.l1, e12.n2, e12.l2,
                                          e12p.n1, e12p.l1, e12p.n2, e12p.l2,
                                          N_sz, gl_w, gl_x, kkn, Bsplines, Ssp, nst_prev, Cf)
                  *wigner_3j0(e12.l1,k,e12p.l1)*wigner_3j0(e12.l2,k,e12p.l2)
                  *wigner_6j(e12p.l1,e12p.l2,L,e12.l2,e12.l1,k);
          }
          /* for exchange check if:
            (-)^{L+la+lc}=(-)^{L+lb+ld},
            |la-ld| <= k <= la+ld,
            |lc-lb| <= k <= lc+lb */
          min_exc = ((L+e12.l1+e12p.l1) >> 0) & 1;
          if(min_exc ==(((L+e12.l2+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l2)<=k) && (k<=e12.l1+e12p.l2)) &&
             ((abs(e12.l2-e12p.l1)<=k) && (k<=e12.l2+e12p.l1))) {
            sum_k += pow(-1,min_exc)*FsltrSimp(k, n, bo, e12.n1, e12.l1, e12.n2, e12.l2,
                                          e12p.n2, e12p.l2, e12p.n1, e12p.l1,
                                          N_sz, gl_w, gl_x, kkn, Bsplines, Ssp, nst_prev, Cf)
                  *wigner_3j0(e12.l1,k,e12p.l2)*wigner_3j0(e12.l2,k,e12p.l1)
                  *wigner_6j(e12p.l1,e12p.l2,L,e12.l1,e12.l2,k);
          }
        }
        std::cout << sum_k<<" "<<pow(-1,(e12.l1+e12.l2))*Y_norm*sum_k<< "\n";
        // write symmetric V_12 as upper triangular
        v_mat[(2*L_sz-NL2-1)*NL2/2 + NL1] = pow(-1,(e12.l1+e12.l2))*Y_norm*sum_k;
      } 
    }
    // save upper triangular V_12
    v_dim[0] = v_sz;
    outfile_name = cpot + "V12_" + std::to_string(L) + ".h5";
    outfile = std::unique_ptr<H5::H5File>(new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    V_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(outfile->createDataSet(
                    "V_12", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, v_dim))));
    V_set->write(&v_mat[0], H5::PredType::NATIVE_DOUBLE);

    v_mat.clear();
  }

  return 0;
}
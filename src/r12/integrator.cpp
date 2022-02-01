#include "integrator.h"

double FsltrGL2(int k, int n, int bo, int gl2,
            // int na, int la, int nb, int lb,
            // int nc, int lc, int nd, int ld,
            std::vector<double> &gl_ow, 
            std::vector<double> &gl_ox,
            std::vector<double> &gl_iw, 
            std::vector<double> &gl_ix, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Gsp,
            // std::vector<int> &offset,
            std::vector<double> &Cl1i_pt,
            std::vector<double> &Cl1p_pt,
            std::vector<double> &Cl2i_pt,
            std::vector<double> &Cl2p_pt) {
  int kp1=k+1;
  // double *Cl1i_pt=&Cf[offset[la]];
  // double *Cl1p_pt=&Cf[offset[lc]];
  // double *Cl2i_pt=&Cf[offset[lb]];
  // double *Cl2p_pt=&Cf[offset[ld]];
  double Pl1i=0, Pl1p=0, Pl2i=0, Pl2p=0;
  double Pl2igli=0, Pl2pgli=0;
  double dl, sl, loc_GL, r2, r1, gl2r, pr2, dli, chi;

  // int nna=na*n, nnb=nb*n, nnc=nc*n, nnd=nd*n;

  // first calculate Qk for all of 0->R
  double Qk=0.0;
  double Jk=0;
  double Fk=0;

  double rm1=0.0;
  double jki=0.0, qki=0.0;
  auto sp=0;
  auto ibo1j=0, jbopi=0;
  auto ai=0;

  for(auto i=bo-1; i<n; ++i) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<bo; ++p){
      r2 = dl*gl_ox[p] + sl;
      Pl2i = 0; Pl2p = 0;

      for(auto j=0; j<bo; ++j) {
        jbopi = j+bo*(p+i*bo);
        ibo1j = i-bo+1+j;
        Pl2i += Cl2i_pt[ibo1j]*Bsp[jbopi];
        Pl2p += Cl2p_pt[ibo1j]*Bsp[jbopi];
      }
      loc_GL+=gl_ow[p]*pow(r2,-kp1)*Pl2i*Pl2p;
    }
    Qk+=dl*loc_GL;
  }

  std::ofstream outFile("dat/slt_test.dat", std::ofstream::out);


  for(auto i=bo-1; i<n; ++i) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    sp=0;
    for(auto p=0; p<bo; ++p, sp+=gl2){
      r1 = dl*gl_ox[p] + sl;
      Pl1i = 0; Pl1p = 0;

      for(auto j=0; j<bo; ++j) {
        ibo1j = i-bo+1+j;
        jbopi = j+bo*(p+i*bo);
        Pl1i += Cl1i_pt[ibo1j]*Bsp[jbopi];
        Pl1p += Cl1p_pt[ibo1j]*Bsp[jbopi];
      }

      // second, inner gaussian quadrature
      jki=0; qki=0;
      dli = (r1-rm1)*0.5;
      for(auto p2=0; p2<gl2; p2+=4) {
        Pl2igli=0; Pl2pgli=0;
        gl2r = dli*gl_ix[p2] + (r1+rm1)*0.5;
        for(auto j=0; j<bo; ++j) {
          ibo1j = i-bo+1+j;
          jbopi = j+bo*(sp+p2+i*gl2*bo);
          ai=ibo1j-(kkn[i]>gl2r);
          Pl2igli += Cl2i_pt[ai]*Gsp[jbopi];
          Pl2pgli += Cl2p_pt[ai]*Gsp[jbopi];
        }
        pr2 = Pl2igli*Pl2pgli;
        jki += gl_iw[p2]*pow(gl2r,k)*pr2;
        qki += gl_iw[p2]*pow(gl2r,-kp1)*pr2;

        Pl2igli=0; Pl2pgli=0;
        gl2r = dli*gl_ix[p2+1] + (r1+rm1)*0.5;
        for(auto j=0; j<bo; ++j) {
          ibo1j = i-bo+1+j;
          jbopi = j+bo*(sp+p2+1+i*gl2*bo);
          ai=ibo1j-(kkn[i]>gl2r);
          Pl2igli += Cl2i_pt[ai]*Gsp[jbopi];
          Pl2pgli += Cl2p_pt[ai]*Gsp[jbopi];
        }
        pr2 = Pl2igli*Pl2pgli;
        jki += gl_iw[p2+1]*pow(gl2r,k)*pr2;
        qki += gl_iw[p2+1]*pow(gl2r,-kp1)*pr2;

        Pl2igli=0; Pl2pgli=0;
        gl2r = dli*gl_ix[p2+2] + (r1+rm1)*0.5;
        for(auto j=0; j<bo; ++j) {
          ibo1j = i-bo+1+j;
          jbopi = j+bo*(sp+p2+2+i*gl2*bo);
          ai=ibo1j-(kkn[i]>gl2r);
          Pl2igli += Cl2i_pt[ai]*Gsp[jbopi];
          Pl2pgli += Cl2p_pt[ai]*Gsp[jbopi];
        }
        pr2 = Pl2igli*Pl2pgli;
        jki += gl_iw[p2+2]*pow(gl2r,k)*pr2;
        qki += gl_iw[p2+2]*pow(gl2r,-kp1)*pr2;

        Pl2igli=0; Pl2pgli=0;
        gl2r = dli*gl_ix[p2+3] + (r1+rm1)*0.5;
        for(auto j=0; j<bo; ++j) {
          ibo1j = i-bo+1+j;
          jbopi = j+bo*(sp+p2+3+i*gl2*bo);
          ai=ibo1j-(kkn[i]>gl2r);
          Pl2igli += Cl2i_pt[ai]*Gsp[jbopi];
          Pl2pgli += Cl2p_pt[ai]*Gsp[jbopi];
        }
        pr2 = Pl2igli*Pl2pgli;
        jki += gl_iw[p2+3]*pow(gl2r,k)*pr2;
        qki += gl_iw[p2+3]*pow(gl2r,-kp1)*pr2;
      }

      // chi(r1)
      Jk+=dli*jki;
      Qk-=dli*qki;
      chi=pow(r1,-kp1)*Jk+pow(r1,k)*Qk;

      outFile <<r1<<" "<<Jk<<" "<<Qk<<" "<<chi<< "\n";

      // Glq outer
      loc_GL += gl_ow[p]*Pl1i*Pl1p*chi;

      rm1=r1;
    }
    Fk+=dl*loc_GL;
  }

  return Fk;
}

// Fast slater integral code with Gauss Lobatto inner int
double FsltrLob4GL(int k, int n, int bo,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Ssp,
            std::vector<double> &Cl1i_pt,
            std::vector<double> &Cl1p_pt,
            std::vector<double> &Cl2i_pt,
            std::vector<double> &Cl2p_pt) {
  int kp1=k+1;
  double Pl1i=0, Pl1p=0, Pl2i=0, Pl2p=0;
  double Pl2ira=0, Pl2pra=0, Pl2irb=0, Pl2prb=0;
  double dl, sl, loc_GL, r2, r1, pr2, pr2a, pr2b, chi;

  // first calculate Qk for all of 0->R
  double Qk=0.0;
  double Jk=0;
  double Fk=0;

  double Lob4o = 0.1666666666666666666666667e0;
  double Lob4i = 0.8333333333333333333333333e0;
  double Lob4p = 0.4472135954999579392818347e0;

  for(auto i=bo-1; i<n; ++i) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<bo; ++p){
      r2 = dl*gl_x[p] + sl;
      Pl2i = 0; Pl2p = 0;

      for(auto j=0; j<bo; ++j) {
        Pl2i += Cl2i_pt[i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
        Pl2p += Cl2p_pt[i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
      }
      loc_GL+=gl_w[p]*pow(r2,-kp1)*Pl2i*Pl2p;
    }
    Qk+=dl*loc_GL;
  }

  // std::ofstream outFile("dat/slt_test"+std::to_string(na)+std::to_string(la)+std::to_string(nb)+std::to_string(lb)
  // +std::to_string(nc)+std::to_string(lc)+std::to_string(nd)+std::to_string(ld)+".dat", std::ofstream::out);

  double rm1=0.0;
  double dlob, slob;
  double rloba, rlobb;
  double fm1j=0.0, fm1q=0.0;
  auto sp=0;
  double pr1k, pr1km;
  auto ibo1j=0, jbopi=0;
  auto ai=0;
  auto bi=0;

  for(auto i=bo-1; i<n; ++i) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    sp=0;
    for(auto p=0; p<bo; ++p, sp+=2){
      r1 = dl*gl_x[p] + sl;
      dlob = (r1-rm1)*0.5;
      slob = (r1+rm1)*0.5;
      rloba = -dlob*Lob4p + slob;
      rlobb = dlob*Lob4p  + slob;
      Pl1i = 0; Pl1p = 0;
      Pl2i = 0; Pl2p = 0;
      Pl2ira=0; Pl2pra=0; 
      Pl2irb=0; Pl2prb=0;

      for(auto j=0; j<bo; ++j) {
        ibo1j = i-bo+1+j;
        jbopi = j+bo*(p+i*bo);
        Pl1i += Cl1i_pt[ibo1j]*Bsp[jbopi];
        Pl1p += Cl1p_pt[ibo1j]*Bsp[jbopi];
        Pl2i += Cl2i_pt[ibo1j]*Bsp[jbopi];
        Pl2p += Cl2p_pt[ibo1j]*Bsp[jbopi];

        // G-Lobatto inner
        ai=ibo1j-(kkn[i]>rloba);
        bi=ibo1j-(kkn[i]>rlobb);
        Pl2ira += Cl2i_pt[ai]*Ssp[j+bo*(sp+i*2*bo)];
        Pl2pra += Cl2p_pt[ai]*Ssp[j+bo*(sp+i*2*bo)];

        Pl2irb += Cl2i_pt[bi]*Ssp[j+bo*(sp+1+i*2*bo)];
        Pl2prb += Cl2p_pt[bi]*Ssp[j+bo*(sp+1+i*2*bo)];
      }
      pr2 = Pl2i*Pl2p;
      pr2a= Pl2ira*Pl2pra;
      pr2b= Pl2irb*Pl2prb;

      // chi(r1)
      pr1k = pow(r1,k);
      pr1km= pow(r1,-kp1);
      Jk+=dlob*(Lob4o*fm1j+Lob4i*pow(rloba,k)*pr2a+Lob4i*pow(rlobb,k)*pr2b+Lob4o*pr1k*pr2);
      Qk-=dlob*(Lob4o*fm1q+Lob4i*pow(rloba,-kp1)*pr2a+Lob4i*pow(rlobb,-kp1)*pr2b+Lob4o*pr1km*pr2);
      chi=pr1km*Jk+pr1k*Qk;

      // outFile <<r1<<" "<<Jk<<" "<<Qk<<" "<<chi<< "\n";

      // Glq outer
      loc_GL += gl_w[p]*Pl1i*Pl1p*chi;

      rm1=r1;
      fm1j=pr1k*pr2;
      fm1q=pr1km*pr2;
    }
    Fk+=dl*loc_GL;
  }

  return Fk;
}

double FsltrLob4GL(int k, int n, int bo, int glq_pt,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Ssp,
            std::vector<double> &rk,
            std::vector<double> &rk_mid,
            std::vector<double> &Cl1i_pt,
            std::vector<double> &Cl2i_pt,
            std::vector<double> &p1p_buff,
            std::vector<double> &p2p_buff,
            std::vector<double> &p2p_mid,
            std::vector<double> &p2is) {
  int kp1=k+1;
  int nqpt=n*glq_pt;
  int nqpt2=2*nqpt;
  int i1, i1bo, i2bo;
  double Pl1i=0, Pl2i=0;
  double Pl2ira=0, Pl2irb=0;
  double dl, sl, loc_GL, r2, r1, pr2, pr2a, pr2b, chi;
  double dlob, slob;
  double rloba, rlobb;
  auto ibo1j=0;

  // first calculate Qk for all of 0->R
  double Qk=0.0;
  double Jk=0;
  double Fk=0;
  double rm1=0.0;
  auto ai=0;
  auto bi=0;
  double fm1q=0.0;
  auto sp=0;

  double Lob4o = 0.1666666666666666666666667e0;
  double Lob4i = 0.8333333333333333333333333e0;
  double Lob4p = 0.4472135954999579392818347e0;

  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    i1bo=(i1-bo)*glq_pt;
    i2bo=i1bo*2;
    dl = (kkn[i1] - kkn[i])*0.5;
    sl = (kkn[i1] + kkn[i])*0.5;

    sp=0;
    for(auto p=0; p<glq_pt; ++p, sp+=2){
      r2 = dl*gl_x[p] + sl;
      dlob = (r2-rm1)*0.5;
      slob = (r2+rm1)*0.5;
      rloba = -dlob*Lob4p + slob;
      rlobb = dlob*Lob4p  + slob;
      Pl2i = 0; 
      Pl2ira=0; Pl2irb=0;

      for(auto j=0; j<bo; ++j) {
        ibo1j=i1-bo+j;
        Pl2i += Cl2i_pt[ibo1j]*Bsp[j+bo*(p+i*glq_pt)];

        // G-Lobatto inner
        ai=ibo1j-(kkn[i]>rloba);
        bi=ibo1j-(kkn[i]>rlobb);
        Pl2ira += Cl2i_pt[ai]*Ssp[j+bo*(sp+i*2*glq_pt)]; // save these also
        Pl2irb += Cl2i_pt[bi]*Ssp[j+bo*(sp+1+i*2*glq_pt)];
      }
      p2is[i1bo+p]=Pl2i;

      Qk+=dlob*(fm1q+Lob4i*(1.0/rk_mid[sp+i2bo+kp1*nqpt2])*Pl2ira*p2p_mid[i2bo+sp]+
                Lob4i*(1.0/rk_mid[sp+1+i2bo+kp1*nqpt2])*Pl2irb*p2p_mid[i2bo+sp+1]+
                Lob4o*(1.0/rk[p+i1bo+kp1*nqpt])*Pl2i*p2p_buff[i1bo+p]);

      fm1q=Lob4o*(1.0/rk[p+i1bo+kp1*nqpt])*Pl2i*p2p_buff[i1bo+p];
      rm1=r2;
    }
  }

  rm1=0.0;
  fm1q=0.0;
  double fm1j=0.0;
  double pr1k, pr1km;

  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    i1bo=(i1-bo)*glq_pt;
    i2bo=i1bo*2;
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    sp=0;
    for(auto p=0; p<glq_pt; ++p, sp+=2){
      r1 = dl*gl_x[p] + sl;
      dlob = (r1-rm1)*0.5;
      slob = (r1+rm1)*0.5;
      rloba = -dlob*Lob4p + slob;
      rlobb = dlob*Lob4p  + slob;
      Pl1i = 0;
      Pl2ira=0;
      Pl2irb=0;

      for(auto j=0; j<bo; ++j) {
        ibo1j = i1-bo+j;
        Pl1i += Cl1i_pt[ibo1j]*Bsp[j+bo*(p+i*glq_pt)];

        // G-Lobatto inner
        ai=ibo1j-(kkn[i]>rloba);
        bi=ibo1j-(kkn[i]>rlobb);
        Pl2ira += Cl2i_pt[ai]*Ssp[j+bo*(sp+i*2*glq_pt)];
        Pl2irb += Cl2i_pt[bi]*Ssp[j+bo*(sp+1+i*2*glq_pt)];
      }
      pr2 = p2is[i1bo+p]*p2p_buff[i1bo+p];
      pr2a= Pl2ira*p2p_mid[i2bo+sp];
      pr2b= Pl2irb*p2p_mid[i2bo+sp+1];

      // chi(r1)
      pr1k = rk[p+i1bo+k*nqpt];
      pr1km= 1.0/rk[p+i1bo+kp1*nqpt];
      Jk+=dlob*(fm1j+Lob4i*rk_mid[sp+i2bo+k*nqpt2]*pr2a+
          Lob4i*rk_mid[sp+1+i2bo+k*nqpt2]*pr2b+Lob4o*pr1k*pr2);
      Qk-=dlob*(fm1q+Lob4i*(1.0/rk_mid[sp+i2bo+kp1*nqpt2])*pr2a+
          Lob4i*(1.0/rk_mid[sp+1+i2bo+kp1*nqpt2])*pr2b+Lob4o*pr1km*pr2);
      chi=pr1km*Jk+pr1k*Qk;

      // Glq outer
      loc_GL += gl_w[p]*Pl1i*p1p_buff[i1bo+p]*chi;

      rm1=r1;
      fm1j=Lob4o*pr1k*pr2;
      fm1q=Lob4o*pr1km*pr2;
    }
    Fk+=dl*loc_GL;
  }

  return Fk;
}

double FsltrLob3GL(int k, int n, int bo, int glq_pt,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Ssp,
            std::vector<double> &rk,
            std::vector<double> &rk_mid,
            std::vector<double> &Cl1i_pt,
            std::vector<double> &Cl2i_pt,
            std::vector<double> &p1p_buff,
            std::vector<double> &p2p_buff,
            std::vector<double> &p2p_mid,
            std::vector<double> &p2is) {
  int kp1=k+1;
  int nqpt=n*glq_pt;
  int i1, i1bo;
  double Pl1i=0, Pl2i=0;
  double Pl2ira=0;
  double dl, sl, loc_GL, r1, pr2, pr2a, chi;
  double rm1=0.0;
  double rlob, dlob;
  double fm1j=0.0, fm1q=0.0;
  auto ai=0;
  double pr1k, pr1km;
  auto ibo1j=0, jbopi=0;

  // first calculate Qk for all of 0->R
  double Qk=0.0;

  constexpr double Lobo = 0.3333333333333333333333333e0;
  constexpr double Lobi = 0.1333333333333333333333333e1;

  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    i1bo=(i1-bo)*glq_pt;
    dl = (kkn[i1] - kkn[i])*0.5;
    sl = (kkn[i1] + kkn[i])*0.5;

    // need to add inner lobatto points
    for(auto p=0; p<glq_pt; ++p){
      r1 = dl*gl_x[p] + sl;
      rlob = (r1+rm1)*0.5;
      dlob = (r1-rm1)*0.5;
      Pl2i = 0;
      Pl2ira=0;
      for(auto j=0; j<bo; ++j) {
        jbopi = j+bo*(p+i*glq_pt);
        Pl2i += Cl2i_pt[i1-bo+j]*Bsp[jbopi];
        
        ai=i1-bo+j-(kkn[i]>rlob);
        Pl2ira += Cl2i_pt[ai]*Ssp[jbopi]; //save
      }
      p2is[i1bo+p]=Pl2i;

      Qk+=dlob*(fm1q+Lobi*(1.0/rk_mid[p+i1bo+kp1*nqpt])*Pl2ira*p2p_mid[i1bo+p]
                +Lobo*(1.0/rk[p+i1bo+kp1*nqpt])*Pl2i*p2p_buff[i1bo+p]);

      fm1q=Lobo*(1.0/rk[p+i1bo+kp1*nqpt])*Pl2i*p2p_buff[i1bo+p];
      rm1=r1;
    }
  }

  double Jk=0.0;
  double Fk=0.0;
  rm1=0.0;
  fm1q=0.0;

  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    i1bo=(i1-bo)*glq_pt;
    dl = (kkn[i1] - kkn[i])*0.5;
    sl = (kkn[i1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<glq_pt; ++p){
      r1 = dl*gl_x[p] + sl;
      rlob = (r1+rm1)*0.5;
      dlob = (r1-rm1)*0.5;
      Pl1i = 0;
      Pl2ira=0;

      for(auto j=0; j<bo; ++j) {
        ibo1j = i1-bo+j;
        jbopi = j+bo*(p+i*glq_pt);
        Pl1i += Cl1i_pt[ibo1j]*Bsp[jbopi];

        // G-Lobatto inner
        ai=ibo1j-(kkn[i]>rlob);
        Pl2ira += Cl2i_pt[ai]*Ssp[jbopi];
      }
      pr2 = p2is[i1bo+p]*p2p_buff[i1bo+p];//p2ps[(i1-bo)*bo+p];
      pr2a= Pl2ira*p2p_mid[i1bo+p]; //Pl2pra;

      // chi(r1)
      pr1k = rk[p+i1bo+k*nqpt];
      pr1km= 1.0/rk[p+i1bo+kp1*nqpt];
      Jk+=dlob*(fm1j+Lobi*rk_mid[p+i1bo+k*nqpt]*pr2a+Lobo*pr1k*pr2);
      Qk-=dlob*(fm1q+Lobi*(1.0/rk_mid[p+i1bo+kp1*nqpt])*pr2a+Lobo*pr1km*pr2);
      chi=pr1km*Jk+pr1k*Qk;

      // Glq outer
      loc_GL+=gl_w[p]*Pl1i*p1p_buff[i1bo+p]*chi;

      rm1=r1;
      fm1j=Lobo*pr1k*pr2;
      fm1q=Lobo*pr1km*pr2;
    }
    Fk+=dl*loc_GL;
  }
  return Fk;
}

double FsltrTrapGL(int k, int n, int bo, int glq_pt,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &rk,
            std::vector<double> &Cl1i_pt,
            std::vector<double> &Cl2i_pt,
            std::vector<double> &p1p_buff,
            std::vector<double> &p2p_buff,
            std::vector<double> &p2is) {
  int kp1=k+1;
  int nqpt=n*glq_pt;
  double Pl1i=0, Pl2i=0;
  double dl, sl, loc_GL, r2, r1, pr2, chi;

  // first calculate Qk for all of 0->R
  double Qk=0.0;
  double Jk=0;
  double Fk=0;
  double rm1=0.0;
  double dltr;
  double fm1j=0.0, fm1q=0.0;
  double pr1k, pr1km;
  auto ibo1j=0, jbopi=0, i1=0, i1bo=0;

  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    dl = (kkn[i1] - kkn[i])*0.5;
    sl = (kkn[i1] + kkn[i])*0.5;
    i1bo=(i1-bo)*glq_pt;

    for(auto p=0; p<glq_pt; ++p){
      r2 = dl*gl_x[p] + sl;
      dltr=(r2-rm1)*0.5;
      Pl2i = 0; Pl2i = 0;

      for(auto j=0; j<bo; ++j) {
        Pl2i += Cl2i_pt[i1-bo+j]*Bsp[j+bo*(p+i*glq_pt)];
      }
      p2is[i1bo+p]=Pl2i;
      Qk+=dltr*(fm1q+(1.0/rk[p+i1bo+kp1*nqpt])*Pl2i*p2p_buff[i1bo+p]);

      fm1q=(1.0/rk[p+i1bo+kp1*nqpt])*Pl2i*p2p_buff[i1bo+p];
      rm1=r2;
    }
  }

  rm1=0.0;
  fm1q=0.0;
  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    i1bo=(i1-bo)*glq_pt;
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<glq_pt; ++p){
      r1 = dl*gl_x[p] + sl;
      dltr=(r1-rm1)*0.5;
      Pl1i = 0;

      for(auto j=0; j<bo; ++j) {
        ibo1j = i1-bo+j;
        jbopi = j+bo*(p+i*glq_pt);
        Pl1i += Cl1i_pt[ibo1j]*Bsp[jbopi];
      }
      pr2 = p2is[i1bo+p]*p2p_buff[i1bo+p];

      // chi(r1)
      pr1k = rk[p+i1bo+k*nqpt];
      pr1km= 1.0/rk[p+i1bo+kp1*nqpt];
      Jk+=dltr*(fm1j+pr1k*pr2);
      Qk-=dltr*(fm1q+pr1km*pr2);
      chi=pr1km*Jk+pr1k*Qk;

      // Glq outer
      loc_GL += gl_w[p]*Pl1i*p1p_buff[i1bo+p]*chi;

      rm1=r1;
      fm1j=pr1k*pr2;
      fm1q=pr1km*pr2;
    }
    Fk+=dl*loc_GL;
  }
  return Fk;
}
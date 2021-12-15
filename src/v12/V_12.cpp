#include "V_12.h"
#include "time_tst.h"
#include <omp.h>

// Fast slater integral code
double Fsltr(int k, int n, int bo,
            int na, int la, int nb, int lb,
            int nc, int lc, int nd, int ld,
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
double FsltrSimpGL(int k, int n, int bo,
            int na, int la, int nb, int lb,
            int nc, int lc, int nd, int ld,
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

  std::ofstream outFile("dat/slt_test"+std::to_string(na)+std::to_string(la)+std::to_string(nb)+std::to_string(lb)
  +std::to_string(nc)+std::to_string(lc)+std::to_string(nd)+std::to_string(ld)+".dat", std::ofstream::out);

  double rm1=0.0;
  double rs38a, rs38b;
  double fm1j=0.0, fm1q=0.0;
  auto sp=0;
  double fact, pr1k, pr1km;
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
      rs38a = (2*rm1+r1)/3;
      rs38b = (rm1+2*r1)/3;
      Pl1i = 0; Pl1p = 0;
      Pl2i = 0; Pl2p = 0;
      Pl2ira=0; Pl2pra=0; 
      Pl2irb=0; Pl2prb=0;

      for(auto j=0; j<bo; ++j) {
        ibo1j = i-bo+1+j;
        jbopi = j+bo*(p+i*bo);
        Pl1i += Cl1i_pt[na*n+ibo1j]*Bsp[jbopi];
        Pl1p += Cl1p_pt[nc*n+ibo1j]*Bsp[jbopi];
        Pl2i += Cl2i_pt[nb*n+ibo1j]*Bsp[jbopi];
        Pl2p += Cl2p_pt[nd*n+ibo1j]*Bsp[jbopi];

        // Simpson's inner
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

      outFile <<r1<<" "<<Jk<<" "<<Qk<<" "<<chi<< "\n";

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

// Fast slater integral code with Gauss Lobatto inner int
double FsltrLob3GL(int k, int n, int bo,
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
  int i1;
  double Pl1i=0, Pl1p=0, Pl2i=0, Pl2p=0;
  double Pl2ira=0, Pl2pra=0;
  double dl, sl, loc_GL, r2, r1, pr2, pr2a, chi;

  // first calculate Qk for all of 0->R
  double Qk=0.0;
  double Jk=0;
  double Fk=0;

  constexpr double Lobo = 0.3333333333333333333333333e0;
  constexpr double Lobi = 0.1333333333333333333333333e1;

  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    dl = (kkn[i1] - kkn[i])*0.5;
    sl = (kkn[i1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<bo; ++p){
      r2 = dl*gl_x[p] + sl;
      Pl2i = 0; Pl2p = 0;

      for(auto j=0; j<bo; ++j) {
        Pl2i += Cl2i_pt[i1-bo+j]*Bsp[j+bo*(p+i*bo)];
        Pl2p += Cl2p_pt[i1-bo+j]*Bsp[j+bo*(p+i*bo)];
      }
      loc_GL+=gl_w[p]*pow(r2,-kp1)*Pl2i*Pl2p;
    }
    Qk+=dl*loc_GL;
  }

  // std::ofstream outFile("dat/slt_test"+std::to_string(na)+std::to_string(la)+std::to_string(nb)+std::to_string(lb)
  // +std::to_string(nc)+std::to_string(lc)+std::to_string(nd)+std::to_string(ld)+".dat", std::ofstream::out);

  double rm1=0.0;
  double rlob, dlob;
  double fm1j=0.0, fm1q=0.0;
  double pr1k, pr1km;
  auto ibo1j=0, jbopi=0;
  auto ai=0;

  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    dl = (kkn[i1] - kkn[i])*0.5;
    sl = (kkn[i1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<bo; ++p){
      r1 = dl*gl_x[p] + sl;
      rlob = (r1+rm1)*0.5;
      dlob = (r1-rm1)*0.5;
      Pl1i = 0; Pl1p = 0;
      Pl2i = 0; Pl2p = 0;
      Pl2ira=0; Pl2pra=0; 

      for(auto j=0; j<bo; ++j) {
        ibo1j = i1-bo+j;
        jbopi = j+bo*(p+i*bo);
        Pl1i += Cl1i_pt[ibo1j]*Bsp[jbopi];
        Pl1p += Cl1p_pt[ibo1j]*Bsp[jbopi];
        Pl2i += Cl2i_pt[ibo1j]*Bsp[jbopi];
        Pl2p += Cl2p_pt[ibo1j]*Bsp[jbopi];

        // G-Lobatto inner
        ai=ibo1j-(kkn[i]>rlob);
        Pl2ira += Cl2i_pt[ai]*Ssp[j+bo*(p+i*bo)];
        Pl2pra += Cl2p_pt[ai]*Ssp[j+bo*(p+i*bo)];
      }
      pr2 = Pl2i*Pl2p;
      pr2a= Pl2ira*Pl2pra;

      // chi(r1)
      pr1k = pow(r1,k);
      pr1km= pow(r1,-kp1);
      Jk+=dlob*(Lobo*fm1j+Lobi*pow(rlob,k)*pr2a+Lobo*pr1k*pr2);
      Qk-=dlob*(Lobo*fm1q+Lobi*pow(rlob,-kp1)*pr2a+Lobo*pr1km*pr2);
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

double FsltrLob3GL(int k, int n, int bo,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Ssp,
            std::vector<double> &rk,
            std::vector<double> &rk_mid,
            std::vector<double> &Cl1i_pt,
            // std::vector<double> &Cl1p_pt,
            std::vector<double> &Cl2i_pt,
            // std::vector<double> &Cl2p_pt,
            std::vector<double> &p1p_buff,
            std::vector<double> &p2p_buff,
            std::vector<double> &p2p_mid,
            std::vector<double> &p2is) {
  int kp1=k+1;
  int i1;
  double Pl1i=0, Pl2i=0;
  double Pl2ira=0;
  double dl, sl, loc_GL, r1, pr2, pr2a, chi;

  // first calculate Qk for all of 0->R
  double Qk=0.0;

  constexpr double Lobo = 0.3333333333333333333333333e0;
  constexpr double Lobi = 0.1333333333333333333333333e1;

  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    dl = (kkn[i1] - kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<bo; ++p){
      Pl2i = 0;
      for(auto j=0; j<bo; ++j) {
        Pl2i += Cl2i_pt[i1-bo+j]*Bsp[j+bo*(p+i*bo)];
        // Pl2p += Cl2p_pt[i1-bo+j]*Bsp[j+bo*(p+i*bo)];
      }
      p2is[(i1-bo)*bo+p]=Pl2i;
      // p2ps[(i1-bo)*bo+p]=Pl2p;
      loc_GL+=gl_w[p]*(1/rk[p+(i1-bo)*bo+kp1*n*bo])*Pl2i*p2p_buff[(i1-bo)*bo+p];
    }
    Qk+=dl*loc_GL;
  }

  std::ofstream outFile;

  outFile.open("dat/slt_test_wrong.dat", std::ofstream::out);


  double rm1=0.0;
  double rlob, dlob;
  double fm1j=0.0, fm1q=0.0;
  double pr1k, pr1km;
  auto ibo1j=0, jbopi=0;
  auto ai=0;

  double Jk=0.0;
  double Fk=0.0;

  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    dl = (kkn[i1] - kkn[i])*0.5;
    sl = (kkn[i1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<bo; ++p){
      r1 = dl*gl_x[p] + sl;
      rlob = (r1+rm1)*0.5;
      dlob = (r1-rm1)*0.5;
      Pl1i = 0;
      // Pl2i = 0; Pl2p = 0;
      Pl2ira=0;

      for(auto j=0; j<bo; ++j) {
        ibo1j = i1-bo+j;
        jbopi = j+bo*(p+i*bo);
        Pl1i += Cl1i_pt[ibo1j]*Bsp[jbopi];
        // Pl1p += Cl1p_pt[ibo1j]*Bsp[jbopi];
        // Pl2i += Cl2i_pt[ibo1j]*Bsp[jbopi];
        // Pl2p += Cl2p_pt[ibo1j]*Bsp[jbopi];

        // G-Lobatto inner
        ai=ibo1j-(kkn[i]>rlob);
        Pl2ira += Cl2i_pt[ai]*Ssp[j+bo*(p+i*bo)];
        // Pl2pra += Cl2p_pt[ai]*Ssp[j+bo*(p+i*bo)];
      }
      pr2 = p2is[(i1-bo)*bo+p]*p2p_buff[(i1-bo)*bo+p];//p2ps[(i1-bo)*bo+p];
      pr2a= Pl2ira*p2p_mid[(i1-bo)*bo+p]; //Pl2pra;

      // chi(r1)
      pr1k = rk[p+(i1-bo)*bo+k*n*bo];
      pr1km= 1.0/rk[p+(i1-bo)*bo+kp1*n*bo];
      Jk+=dlob*(Lobo*fm1j+Lobi*rk_mid[p+(i1-bo)*bo+k*n*bo]*pr2a+Lobo*pr1k*pr2);
      Qk-=dlob*(Lobo*fm1q+Lobi*(1.0/rk_mid[p+(i1-bo)*bo+kp1*n*bo])*pr2a+Lobo*pr1km*pr2);
      chi=pr1km*Jk+pr1k*Qk;

      // if(flag==1) {
      //   outFile <<r1<<" "<<chi<<" "<<k<<" "<<Fk<<" "
      //           <<p2is[(i1-bo)*bo+p]<<" "<<p2p_buff[(i1-bo)*bo+p]<<" "
      //           <<Pl1i<<" "<<p1p_buff[(i1-bo)*bo+p]<<"\n";
      // }

      // Glq outer
      loc_GL += gl_w[p]*Pl1i*p1p_buff[(i1-bo)*bo+p]*chi;

      rm1=r1;
      fm1j=pr1k*pr2;
      fm1q=pr1km*pr2;
    }
    Fk+=dl*loc_GL;
    outFile <<r1<<" "<<chi<<" "<<k<<" "<<Fk<<"\n";

  }
  //if(flag==1) flag=0;

  return Fk;
}

double FsltrTrapGL(int k, int n, int bo,
            std::vector<double> &gl_w, 
            std::vector<double> &gl_x, 
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Cl1i_pt,
            std::vector<double> &Cl1p_pt,
            std::vector<double> &Cl2i_pt,
            std::vector<double> &Cl2p_pt) {
  int kp1=k+1;
  double Pl1i=0, Pl1p=0, Pl2i=0, Pl2p=0;
  double dl, sl, loc_GL, r2, r1, pr2, chi;

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
  double dltr;
  double fm1j=0.0, fm1q=0.0;
  double pr1k, pr1km;
  auto ibo1j=0, jbopi=0;

  for(auto i=bo-1; i<n; ++i) {
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<bo; ++p){
      r1 = dl*gl_x[p] + sl;
      dltr=(r1-rm1)*0.5;
      Pl1i = 0; Pl1p = 0;
      Pl2i = 0; Pl2p = 0;

      for(auto j=0; j<bo; ++j) {
        ibo1j = i-bo+1+j;
        jbopi = j+bo*(p+i*bo);
        Pl1i += Cl1i_pt[ibo1j]*Bsp[jbopi];
        Pl1p += Cl1p_pt[ibo1j]*Bsp[jbopi];
        Pl2i += Cl2i_pt[ibo1j]*Bsp[jbopi];
        Pl2p += Cl2p_pt[ibo1j]*Bsp[jbopi];
      }
      pr2 = Pl2i*Pl2p;

      // chi(r1)
      pr1k = pow(r1,k);
      pr1km= pow(r1,-kp1);
      Jk+=dltr*(fm1j+pr1k*pr2);
      Qk-=dltr*(fm1q+pr1km*pr2);
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

double FsltrTrap(int k, int n, int bo, int pt,
            std::vector<double> &kkn,
            std::vector<double> &Bsp,
            std::vector<double> &Cl1i_pt,
            std::vector<double> &Cl1p_pt,
            std::vector<double> &Cl2i_pt,
            std::vector<double> &Cl2p_pt) {
  int kp1=k+1;
  double Pl1i=0, Pl1p=0, Pl2i=0, Pl2p=0;
  double a, b, step, loc_GL, r2, r1, pr2, chi;
  double dltr;

  // first calculate Qk for all of 0->R
  double Qk=0.0;
  double Jk=0;
  double Fk=0;
  double rm1=0.0;
  double fm1j=0.0, fm1q=0.0;
  double pr1k, pr1km;

  for(auto i=bo-1; i<n; ++i) {
    b = kkn[i+1];
    a = kkn[i];
    step = (b-a)/(double)pt;
    loc_GL=0.0;

    for(auto p=1; p<=pt; ++p){
      r2 = a + p*step;
      dltr = (r2-rm1)*0.5;
      Pl2i = 0; Pl2p = 0;

      for(auto j=0; j<bo; ++j) {
        Pl2i += Cl2i_pt[i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
        Pl2p += Cl2p_pt[i-bo+1+j]*Bsp[j+bo*(p+i*bo)];
      }
      pr1km=pow(r2,-kp1);
      loc_GL+=dltr*(fm1q + pr1km*Pl2i*Pl2p);
      rm1=r2;
      fm1q=pr1km*Pl2i*Pl2p;
    }
    Qk+=loc_GL;
  }

  // std::ofstream outFile("dat/slt_test"+std::to_string(na)+std::to_string(la)+std::to_string(nb)+std::to_string(lb)
  // +std::to_string(nc)+std::to_string(lc)+std::to_string(nd)+std::to_string(ld)+".dat", std::ofstream::out);

  rm1=0.0;
  
  fm1j=0.0, fm1q=0.0;
  auto ibo1j=0, jbopi=0;

  for(auto i=bo-1; i<n; ++i) {
    b = kkn[i+1];
    a = kkn[i];
    step = (b-a)/(double)pt;
    loc_GL=0.0;

    for(auto p=1; p<=pt; ++p){
      r1 = a + p*step;
      dltr=(r1-rm1)*0.5;
      Pl1i = 0; Pl1p = 0;
      Pl2i = 0; Pl2p = 0;

      for(auto j=0; j<bo; ++j) {
        ibo1j = i-bo+1+j;
        jbopi = j+bo*(p+i*bo);
        Pl1i += Cl1i_pt[ibo1j]*Bsp[jbopi];
        Pl1p += Cl1p_pt[ibo1j]*Bsp[jbopi];
        Pl2i += Cl2i_pt[ibo1j]*Bsp[jbopi];
        Pl2p += Cl2p_pt[ibo1j]*Bsp[jbopi];
      }
      pr2 = Pl2i*Pl2p;

      // chi(r1)
      pr1k = pow(r1,k);
      pr1km= pow(r1,-kp1);
      Jk+=dltr*(fm1j+pr1k*pr2);
      Qk-=dltr*(fm1q+pr1km*pr2);
      chi=pr1km*Jk+pr1k*Qk;

      // outFile <<r1<<" "<<Jk<<" "<<Qk<<" "<<chi<< "\n";

      // Glq outer
      loc_GL += Pl1i*Pl1p*chi;

      rm1=r1;
      fm1j=pr1k*pr2;
      fm1q=pr1km*pr2;
    }
    Fk+=loc_GL;
  }

  return Fk;
}

// alternative slater integral code
double Fsltr_alt(int k, int n, int bo,
            int na, int la, int nb, int lb,
            int nc, int lc, int nd, int ld,
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

    for(int p=0; p<bo; ++p){
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

int Rpowk(std::vector<double> &r_out, std::vector<double> &r_in, int n, int bo,
       std::vector<double> &gl_x, std::vector<double> &kkn, int k_max) {
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
      ri = dl*gl_x[p] + sl;
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
        ri = dl*gl_x[p] + sl;
        rlm = (rm1+ri)*0.5;
        r_out[p+ibo1bo+k*nbo] = ri*r_out[p+ibo1bo+km1*nbo];
        r_in[p+ibo1bo+k*nbo] = rlm*r_in[p+ibo1bo+km1*nbo];
        rm1=ri;
      }
    }
  }

  return 0;
}

int Rpowk(std::vector<double> &r_out, std::vector<double> &r_in, int n, int bo,
       std::vector<double> &gl_xo, std::vector<double> &gl_xi,
       std::vector<double> &kkn, int k_max) {
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

int Prprim(int n, int bo, std::vector<double> &kkn,
           std::vector<double> &gl_x, int per_l, 
           int n1, int l1, int n2, int l2,
           std::vector<double> &Bsp, std::vector<double> &Ssp, 
           std::vector<double> &Cf, 
           std::vector<double> &p1p_out, std::vector<double> &p2p_out,
           std::vector<double> &p1p_in, std::vector<double> &p2p_in) {
      double rm1 = 0.0, r=0.0;
      double dl, sl;
      for(auto i=bo-1; i<n; ++i) {
        auto i1=i+1;
        dl = (kkn[i1] - kkn[i])*0.5;
        sl = (kkn[i1] + kkn[i])*0.5;

        for(auto p=0; p<bo; ++p){
          double Pl1p = 0, Pl2p = 0, Pl2pm=0, Pl1pm=0;
          r=dl*gl_x[p] + sl;
          double rlob = (rm1+r)*0.5;

          for(auto j=0; j<bo; ++j) {
            auto ibo1j=i1-bo+j;
            Pl1p += Cf[per_l*l1+n1*n+ibo1j]*Bsp[j+bo*(p+i*bo)];
            Pl2p += Cf[per_l*l2+n2*n+ibo1j]*Bsp[j+bo*(p+i*bo)];

            auto ai=ibo1j-(kkn[i]>rlob);
            Pl1pm += Cf[per_l*l1+n1*n+ai]*Ssp[j+bo*(p+i*bo)];
            Pl2pm += Cf[per_l*l2+n2*n+ai]*Ssp[j+bo*(p+i*bo)];
          }
          p1p_out[(i1-bo)*bo+p]=Pl1p;
          p2p_out[(i1-bo)*bo+p]=Pl2p;
          p1p_in[(i1-bo)*bo+p] =Pl1pm;
          p2p_in[(i1-bo)*bo+p] =Pl2pm;
          rm1=r;
        }
      }

  return 0;
}

int V12(std::string cpot, int L_max, std::vector<uint> &N_sz) {
  int min_dir=0, min_exc=0;
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

  int L_sz=0, v_sz=0;
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
  for(int i=1; i<=L_max; ++i) {
    nt += N_sz[i-1];
    nst_prev[i] = nt*n;
    C[i] = &Cf[nt*n];
  }

  // read coefficients for all l
  for(int l=0; l<=L_max; ++l) {
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

  int bo2=3; // 4 inner 13 outer gives epsilon for linear
  std::vector<double> gl_ix(bo2);
  std::vector<double> gl_iw(bo2);
  for(int i=1; i<=bo2; ++i) {
    gl_i = fastgl::GLPair(bo2, i);
    gl_ix[bo2-i] = gl_i.x(); 
    gl_iw[bo2-i] = gl_i.weight;
  }

  // generate B-splines
  std::vector<double> Bsplines;
  bsp::Splines(n, bo, gl_x, kkn, Bsplines);

  std::vector<double> Ssp;
  // bsp::SimpSplines(n, bo, gl_x, kkn, Ssp);
  bsp::GL2Splines(n, bo, bo2, gl_x, gl_ix, kkn, Ssp);
  // bsp::LobSplines(n, bo, gl_x, kkn, Ssp);

  std::vector<double> l1p_loc(n), l2p_loc(n);

  // int L_real_size=0;
  // L_sz = 10000; // read only N=1 and N=2 for each L

  omp_lock_t copylock;
  omp_init_lock(&copylock);

  for(int L=0; L<=L_max; ++L) {
    
    // Read n1l1;n2l2 indices for NL states
    filename = cpot + std::to_string(L) + "idx.h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
    L_sz = L_set->getSpace().getSimpleExtentNpoints()/4;
    L_idx.resize(L_sz);
    L_set->read(&L_idx[0], H5::PredType::NATIVE_UINT32);
    file->close();

    std::cout << "L: " << L <<" L_sz: "<< L_sz << "\n";

    v_sz = L_sz*(L_sz+1)/2;
    v_mat.reserve(v_sz);

    uint64_t st_time = GetTimeMs64();    
    // #pragma omp parallel private(e12p)
    // {
    for(int NL2=0; NL2<L_sz; ++NL2) {
      //set n1'l1';n2'l2'
      e12p = L_idx[NL2];
      if (e12p.n1==29&&e12p.n2==29&&e12p.l1==0) {
      std::cout << e12p.n1<<" "<<e12p.l1<<" "<<e12p.n2<<" "<<e12p.l2 <<"\n";}
      //omp_set_lock(&copylock);
      #pragma omp critical
      {
      std::copy_n(Cf.begin()+offset[e12p.l1]+e12p.n1*n, n, l1p_loc.begin());
      std::copy_n(Cf.begin()+offset[e12p.l2]+e12p.n2*n, n, l2p_loc.begin());
      }
      //omp_unset_lock(&copylock);
      // #pragma omp for private(e12,Y_norm,sum_k,min_dir,min_exc)
      if (e12p.n1==29&&e12p.n2==29&&e12p.l1==0) {
      for(int NL1=NL2; NL1<L_sz; ++NL1){
        //set n1l1;n2l2
        e12 = L_idx[NL1];
        if (e12.n1==29&&e12.n2==29&&e12.l1==0){
        std::cout <<e12.n1<<" "<<e12.l1<<" "<<e12.n2<<" "<<e12.l2 <<"\n";
        std::vector<double> l1i_loc(n), l2i_loc(n);
        omp_set_lock(&copylock);
        std::copy_n(Cf.begin()+offset[e12.l1]+e12.n1*n, n, l1i_loc.begin());
        std::copy_n(Cf.begin()+offset[e12.l2]+e12.n2*n, n, l2i_loc.begin());
        omp_unset_lock(&copylock);
        // sqrt([la][lc][lb][ld])
        Y_norm = sqrt((2*e12.l1+1)*(2*e12p.l1+1)*(2*e12.l2+1)*(2*e12p.l2+1));
        sum_k=0.0;
        for(int k=0; k<=std::max(std::max(e12.l1+e12p.l1,e12.l2+e12p.l2),std::max(e12.l1+e12p.l2,e12.l2+e12p.l1)); ++k) { 
          /* for direct check if:
            (-)^{L+lb+lc}=(-)^{L+la+ld},
            |la-lc| <= k <= la+lc,
            |lb-ld| <= k <= lb+ld */
          min_dir = ((L+e12.l2+e12p.l1) >> 0) & 1; // check if L+lb+lc is even
          if(min_dir ==(((L+e12.l1+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l1)<=k) && (k<=e12.l1+e12p.l1)) &&
             ((abs(e12.l2-e12p.l2)<=k) && (k<=e12.l2+e12p.l2))) {
            sum_k += pow(-1,min_dir)*FsltrGL2(k, n, bo, bo2, //e12.n1, e12.l1, e12.n2, e12.l2,
                                          //e12p.n1, e12p.l1, e12p.n2, e12p.l2,
                                          gl_w, gl_x, gl_iw, gl_ix, kkn, Bsplines, Ssp, 
                                          l1i_loc, l1p_loc, l2i_loc, l2p_loc)
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
            sum_k += pow(-1,min_exc)*FsltrGL2(k, n, bo, bo2, //e12.n1, e12.l1, e12.n2, e12.l2,
                                          // e12p.n2, e12p.l2, e12p.n1, e12p.l1,
                                          gl_w, gl_x, gl_iw, gl_ix, kkn, Bsplines, Ssp,
                                          l1i_loc, l2p_loc, l2i_loc, l1p_loc)
                  *wigner_3j0(e12.l1,k,e12p.l2)*wigner_3j0(e12.l2,k,e12p.l1)
                  *wigner_6j(e12p.l1,e12p.l2,L,e12.l1,e12.l2,k);
          }
        }

        std::cout << std::setiosflags(std::ios::scientific)
                  << std::setprecision(15) << pow(-1,(e12.l1+e12.l2))*Y_norm*sum_k<< "\n";
        // write symmetric V_12 as upper triangular
        v_mat[(2*L_sz-NL2-1)*NL2/2 + NL1] = pow(-1,(e12.l1+e12.l2))*Y_norm*sum_k;}
      } }
    }
    // }
    omp_destroy_lock(&copylock);
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

  return 0;
}

int V12(std::string cpot, int L_max, std::string dir) {
  int min_dir=0, min_exc=0;
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

  int L_sz=0, v_sz=0;
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
  block[1] =1;
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

  int ncf, sym;
  std::vector<cfg::line> cfgs;
  cfg::ReadCfg(dir, 0, sym, ncf, cfgs);

  auto max_n2l = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
  auto max_N = max_n2l.n2max;

  auto max_line = *std::max_element(cfgs.begin(), cfgs.end(),
      [](cfg::line const &a, cfg::line const &b) { return a.l2 < b.l2; });
  auto l_m = max_line.l2;

  // reserve space for coefficients
  Cf.reserve(max_N*n*l_m);
  C[0]=&Cf[0];
  int nt=0;
  std::vector<int> nst_prev(L_max+1);
  nst_prev[0] = nt;
  for(int i=1; i<=L_max; ++i) {
    nt += max_N;
    nst_prev[i] = nt*n;
    C[i] = &Cf[nt*n];
  }

  // read coefficients for all l
  for(int l=0; l<=L_max; ++l) {
    count[0] = max_N; 
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

  // int bo2=9; // 4 inner 13 outer gives epsilon for linear
  // std::vector<double> gl_ix(bo2);
  // std::vector<double> gl_iw(bo2);
  // for(int i=1; i<=bo2; ++i) {
  //   gl_i = fastgl::GLPair(bo2, i);
  //   gl_ix[bo2-i] = gl_i.x(); 
  //   gl_iw[bo2-i] = gl_i.weight;
  // }

  // generate B-splines
  std::vector<double> Bsplines;
  bsp::Splines(n, bo, gl_x, kkn, Bsplines);

  std::vector<double> Ssp;
  // bsp::SimpSplines(n, bo, gl_x, kkn, Ssp);
  // bsp::GL2Splines(n, bo, bo2, gl_x, gl_ix, kkn, Ssp);
  bsp::Lob3Splines(n, bo, gl_x, kkn, Ssp);

  std::vector<double> l1p_loc(n), l2p_loc(n);

  // int L_real_size=0;
  // L_sz = 10000; // read only N=1 and N=2 for each L

  omp_lock_t copylock;
  omp_init_lock(&copylock);

  for(int L=0; L<=L_max; ++L) {
    
    // Read n1l1;n2l2 indices for NL states
    filename = cpot + std::to_string(L) + "idx.h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
    L_sz = L_set->getSpace().getSimpleExtentNpoints()/4;
    L_idx.resize(L_sz);
    L_set->read(&L_idx[0], H5::PredType::NATIVE_UINT32);

    std::cout << "L: " << L <<" L_sz: "<< L_sz << "\n";

    v_sz = L_sz*(L_sz+1)/2;
    v_mat.reserve(v_sz);

    uint64_t st_time = GetTimeMs64();    
    #pragma omp parallel private(e12p)
    {
    for(int NL2=0; NL2<L_sz; ++NL2) {
      //set n1'l1';n2'l2'
      e12p = L_idx[NL2];
      //omp_set_lock(&copylock);
      #pragma omp critical
      {
      std::copy_n(Cf.begin()+offset[e12p.l1]+e12p.n1*n, n, l1p_loc.begin());
      std::copy_n(Cf.begin()+offset[e12p.l2]+e12p.n2*n, n, l2p_loc.begin());
      }
      //omp_unset_lock(&copylock);
      #pragma omp for private(e12,Y_norm,sum_k,min_dir,min_exc)
      for(int NL1=NL2; NL1<L_sz; ++NL1) {
        //set n1l1;n2l2
        e12 = L_idx[NL1];

        std::vector<double> l1i_loc(n), l2i_loc(n);
        omp_set_lock(&copylock);
        std::copy_n(Cf.begin()+offset[e12.l1]+e12.n1*n, n, l1i_loc.begin());
        std::copy_n(Cf.begin()+offset[e12.l2]+e12.n2*n, n, l2i_loc.begin());
        omp_unset_lock(&copylock);
        // sqrt([la][lc][lb][ld])
        Y_norm = sqrt((2*e12.l1+1)*(2*e12p.l1+1)*(2*e12.l2+1)*(2*e12p.l2+1));
        sum_k=0.0;
        for(int k=0; k<=std::max(std::max(e12.l1+e12p.l1,e12.l2+e12p.l2),std::max(e12.l1+e12p.l2,e12.l2+e12p.l1)); ++k) { 
          /* for direct check if:
            (-)^{L+lb+lc}=(-)^{L+la+ld},
            |la-lc| <= k <= la+lc,
            |lb-ld| <= k <= lb+ld */
          min_dir = ((L+e12.l2+e12p.l1) >> 0) & 1; // check if L+lb+lc is even
          if(min_dir ==(((L+e12.l1+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l1)<=k) && (k<=e12.l1+e12p.l1)) &&
             ((abs(e12.l2-e12p.l2)<=k) && (k<=e12.l2+e12p.l2))) {
            sum_k += pow(-1,min_dir)*FsltrLob3GL(k, n, bo, //e12.n1, e12.l1, e12.n2, e12.l2,
                                          //e12p.n1, e12p.l1, e12p.n2, e12p.l2,
                                          gl_w, gl_x, kkn, Bsplines, Ssp, 
                                          l1i_loc, l1p_loc, l2i_loc, l2p_loc)
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
            sum_k += pow(-1,min_exc)*FsltrLob3GL(k, n, bo, //e12.n1, e12.l1, e12.n2, e12.l2,
                                          // e12p.n2, e12p.l2, e12p.n1, e12p.l1,
                                          gl_w, gl_x, kkn, Bsplines, Ssp,
                                          l1i_loc, l2p_loc, l2i_loc, l1p_loc)
                  *wigner_3j0(e12.l1,k,e12p.l2)*wigner_3j0(e12.l2,k,e12p.l1)
                  *wigner_6j(e12p.l1,e12p.l2,L,e12.l1,e12.l2,k);
          }
        }
        // if(e12p.n1==0&&e12p.l1==0&&e12p.n2==0&&e12p.l2==0&&
        //   e12.n1==0&&e12.l1==0&&e12.n2==0&&e12.l2==0) {
        // std::cout << std::setiosflags(std::ios::scientific)
        //           << std::setprecision(15) << pow(-1,(e12.l1+e12.l2))*Y_norm*sum_k<< "\n";
        // }
        // write symmetric V_12 as upper triangular
        v_mat[(2*L_sz-NL2-1)*NL2/2 + NL1] = pow(-1,(e12.l1+e12.l2))*Y_norm*sum_k;
      }
    }
    }
    omp_destroy_lock(&copylock);
    std::cout << "loop time: " << ((double)(GetTimeMs64()-st_time)/1000.0) << "s\n";
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

int V12_alt(std::string cpot, int L_max, std::string dir) {
  bool min_dir=0, min_exc=0;
  double Y_norm=0.0;
  std::vector<double> v_mat;
  double sum_k=0.0;
  idx4 e12, e12p;

  int n=0;   // no. of points
  int bo=0;   // max B-spline order
  int nkn=0; // no. of knots
  std::vector<double> kkn;
  std::vector<double> Cf;
  std::vector<double*> C(L_max+1); // error

  int L_sz=0, v_sz=0;
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
  block[1] =1;
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

  // reserve space for coefficients
  auto e1_lm = std::max(l1_m,l2_m);
  Cf.reserve(max_N*n*(e1_lm+1));
  C[0]=&Cf[0];
  int nt=0;
  std::vector<int> nst_prev(L_max+1);
  // nst_prev[0] = nt;
  // This could be wrong !!!
  for(int i=1; i<=e1_lm; ++i) {
    // nst_prev[i] = nt*n;
    C[i] = &Cf[max_N*n*i]; //invalid read of size 8
  }

  // read coefficients for all l
  for(int l=0; l<=e1_lm; ++l) {
    count[0] = max_N;
    count[1] = n;
    dimms[0] = count[0]; 
    dimms[1] = count[1];
    memspace.setExtentSimple(2, dimms, NULL);

    filename = cpot + std::to_string(l) + ".h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    rset = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("Coeff")));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(C[l], H5::PredType::NATIVE_DOUBLE, memspace, cspace); //invalid read of size 8
    file->close();
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

  // Try to see if this is better!!!

  // int bo2=4; // 4 inner 13 outer gives epsilon for linear
  // std::vector<double> gl_ix(bo2);
  // std::vector<double> gl_iw(bo2);
  // for(int i=1; i<=bo2; ++i) {
  //   gl_i = fastgl::GLPair(bo2, i);
  //   gl_ix[bo2-i] = gl_i.x(); 
  //   gl_iw[bo2-i] = gl_i.weight;
  // }

  // generate B-splines
  std::vector<double> Bsplines;
  bsp::Splines(n, bo, gl_x, kkn, Bsplines);

  std::vector<double> Ssp;
  // bsp::SimpSplines(n, bo, gl_x, kkn, Ssp);
  // bsp::GL2Splines(n, bo, bo2, gl_x, gl_ix, kkn, Ssp);
  bsp::Lob3Splines(n, bo, gl_x, kkn, Ssp);

  // Precalculate r^k
  std::vector<double> rk, rk_mid;
  Rpowk(rk, rk_mid, n, bo, gl_x, kkn, k_max);

  std::vector<double> p1p_buff(bo*n), p2p_buff(bo*n), p2p_mid(bo*n), p1p_mid(bo*n);

  // int L_real_size=0;
  // L_sz = 100; // read only N=1 and N=2 for each L

  omp_lock_t copylock;
  omp_init_lock(&copylock);
  omp_set_num_threads(1);
  uint64_t st_time;

  #pragma omp parallel private(e12p)
  {
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
    // #pragma omp parallel private(e12p)
    // {
    for(int NL2=0; NL2<1; ++NL2) {
      //set n1'l1';n2'l2'
      e12p = {0,4,7,4};//L_idx[NL2];
      //omp_set_lock(&copylock);
      #pragma omp single
      {
      // calculate P1' P2' here and pass it to sltr
      Prprim(n, bo, kkn, gl_x, max_N, e12p.n1, e12p.l1, e12p.n2, e12p.l2, 
          Bsplines, Ssp, Cf, p1p_buff, p2p_buff, p1p_mid, p2p_mid);
      }

      //omp_unset_lock(&copylock);
      #pragma omp barrier
      #pragma omp for private(e12,Y_norm,sum_k,min_dir,min_exc)
      for(int NL1=NL2; NL1<1; ++NL1) {
        //set n1l1;n2l2
        e12 = {4,0,20,0};//L_idx[NL1];

        std::vector<double> l1i_loc(n), l2i_loc(n);
        std::vector<double> pi(bo*n), pp(bo*n);
        omp_set_lock(&copylock);
        std::copy_n(Cf.begin()+e12.l1*max_N*n+e12.n1*n, n, l1i_loc.begin());
        std::copy_n(Cf.begin()+e12.l2*max_N*n+e12.n2*n, n, l2i_loc.begin());
        omp_unset_lock(&copylock);
        // sqrt([la][lc][lb][ld])
        Y_norm = sqrt((2*e12.l1+1)*(2*e12p.l1+1)*(2*e12.l2+1)*(2*e12p.l2+1));
        sum_k=0.0;
        min_dir = ((L+e12.l2+e12p.l1) >> 0) & 1; // check if L+lb+lc is odd
        min_exc = ((L+e12.l1+e12p.l1) >> 0) & 1;
        for(int k=std::max(std::max(abs(e12.l1-e12p.l1),abs(e12.l2-e12p.l2)),std::max(abs(e12.l1-e12p.l2),abs(e12.l2-e12p.l1))); 
            k<=std::min(std::min(e12.l1+e12p.l1,e12.l2+e12p.l2),std::max(e12.l1+e12p.l2,e12.l2+e12p.l1)); ++k) { 
          /* for direct check if:
            (-)^{L+lb+lc}=(-)^{L+la+ld},
            |la-lc| <= k <= la+lc,
            |lb-ld| <= k <= lb+ld */
          if(min_dir ==(((L+e12.l1+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l1)<=k) && (k<=e12.l1+e12p.l1)) &&
             ((abs(e12.l2-e12p.l2)<=k) && (k<=e12.l2+e12p.l2))) {
            sum_k += //pow(-1,min_dir)*
            FsltrLob3GL(k, n, bo, gl_w, gl_x, kkn, Bsplines, 
                          Ssp, rk, rk_mid, l1i_loc, l2i_loc, p1p_buff, 
                          p2p_buff, p2p_mid, pi);
                          //*wigner_3j0(e12.l1,k,e12p.l1)
                          //*wigner_3j0(e12.l2,k,e12p.l2)
                          //*wigner_6j(e12p.l1,e12p.l2,L,e12.l2,e12.l1,k);
            std::cout<< sum_k<<"\n";
          }
          /* for exchange check if:
            (-)^{L+la+lc}=(-)^{L+lb+ld},
            |la-ld| <= k <= la+ld,
            |lc-lb| <= k <= lc+lb */
          if(min_exc ==(((L+e12.l2+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l2)<=k) && (k<=e12.l1+e12p.l2)) &&
             ((abs(e12.l2-e12p.l1)<=k) && (k<=e12.l2+e12p.l1))) {
            sum_k += //pow(-1,min_exc)*
            FsltrLob3GL(k, n, bo, gl_w, gl_x, kkn, Bsplines, 
                          Ssp, rk, rk_mid, l1i_loc, l2i_loc, p2p_buff, 
                          p1p_buff, p1p_mid, pi);
                          //*wigner_3j0(e12.l1,k,e12p.l2)
                          //*wigner_3j0(e12.l2,k,e12p.l1)
                          //*wigner_6j(e12p.l1,e12p.l2,L,e12.l1,e12.l2,k);
            std::cout<< sum_k<<"\n";
          }
        }
        // omp_set_lock(&copylock);
        //   if(abs(sum_k)>1.0e+5) {
        //   std::cout<<e12p.n1<<" "<<e12p.l1<<" "<<e12p.n2<<" "<<e12p.l2<<
        //   " "<<e12.n1<<" "<<e12.l1<<" "<<e12.n2<<" "<<e12.l2<<"\n"; }
        // omp_unset_lock(&copylock);
        // if(e12p.n1==0&&e12p.l1==0&&e12p.n2==0&&e12p.l2==0&&
        //   e12.n1==3&&e12.l1==1&&e12.n2==105&&e12.l2==1) {
        // std::cout << std::setiosflags(std::ios::scientific)
        //           << std::setprecision(15) << pow(-1,(e12.l1+e12.l2))*Y_norm*sum_k<< "\n";
        //   }
        // write symmetric V_12 as upper triangular
        v_mat[(2*L_sz-NL2-1)*NL2/2 + NL1] = pow(-1,(e12.l1+e12.l2))*Y_norm*sum_k;
      }
    }

    #pragma omp single
    {
    // omp_destroy_lock(&copylock);
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
  }

  omp_destroy_lock(&copylock);

  return 0;
}
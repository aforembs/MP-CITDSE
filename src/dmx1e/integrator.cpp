#include "integrator.h"

double tvelGL(int n, int glq_pt,
              int n1, int l1, 
              int n2, int l2,
              std::vector<double> &gl_w,
              std::vector<double> &gl_x,
              std::vector<double> &kkn,
              std::vector<double> &wfn,
              std::vector<double> &wfnp) {
  double t_ab=0.0;
  int nqpt=n*glq_pt;
  int p1=l1*lc_sz+n1*nqpt;
  int p2=l2*lc_sz+n2*nqpt;
  auto i1=0, i1bo=0;
  double v2=(l1*(l1+1)-l2*(l2+1))*0.5;
  double dl, sl, loc_GL, r1, pra;

  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    i1bo=(i1-bo)*glq_pt;
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<glq_pt; ++p){
      r1 = dl*gl_x[p] + sl;
      pra = wfn[p1+p+i1bo];

      loc_GL+=gl_w[p]*(pra*wfnp[n2*nqpt+p+i1bo]-v2*pra*wfn[p2+p+i1bo]/r1);
    }
    t_ab+=dl*loc_GL;
  }
  return t_ab;
}

double tlenGL(int n, int glq_pt,
              int n1, int l1, 
              int n2, int l2,
              std::vector<double> &gl_w,
              std::vector<double> &gl_x,
              std::vector<double> &kkn,
              std::vector<double> &wfn) {
  double t_ab=0.0;
  int nqpt=n*glq_pt;
  int p1=l1*lc_sz+n1*nqpt;
  int p2=l2*lc_sz+n2*nqpt;
  auto i1=0, i1bo=0;
  double dl, sl, loc_GL, r1;

  for(auto i=bo-1; i<n; ++i) {
    i1=i+1;
    i1bo=(i1-bo)*glq_pt;
    dl = (kkn[i+1] - kkn[i])*0.5;
    sl = (kkn[i+1] + kkn[i])*0.5;
    loc_GL=0.0;

    for(auto p=0; p<glq_pt; ++p){
      r1 = dl*gl_x[p] + sl;

      loc_GL+=gl_w[p]*r1*wfn[p1+p+i1bo]*wfn[p2+p+i1bo];
    }
    t_ab+=dl*loc_GL;
  }

  return t_ab;
}
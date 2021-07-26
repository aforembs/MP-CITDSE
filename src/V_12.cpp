#include "V_12.h"

int V12() {
  uint min_dir=0, min_exc=0;
  double Y_norm=0.0;
  double F_dir=1.0; // dummy slater integrals
  double F_exc=F_dir;
  double v_mat; // standin for V_12 matrix
  double sum_k=0.0;

  for(uint L=0; L<=L_max; ++L) {
    // Read indices n1l1;n2l2

    for(uint NL1=0; NL1<sz_L; ++NL1) {
      //set n1l1;n2l2
      for(uint NL2=0; NL2<=NL1; ++NL2){
        //set n1'l1';n2'l2'
        Y_norm = sqrt((2*l1+1)*(2*l1p+1)*(2*l2+1)*(2*l2p+1));
        sum_k=0.0;
        for(uint k=0; k<=l1e_max; ++k) {
          min_dir = ((L+l2+l1p) >> 0) & 1;
          if(min_dir ==(((L+l1+l2p) >> 0) & 1) &&
             ((abs(l1-l1p)<=k) && (k<=l1+l1p)) &&
             ((abs(l2-l2p)<=k) && (k<=l2+l2p))) {
            sum_k += pow(-1,min_dir)*F_dir*wigner_3j0(l1,k,l1p)*wigner_3j0(l2,k,l2p)
                    *wigner_6j(l1p,l2p,L,l2,l1,k);
          }
          min_exc = ((L+l1+l1p) >> 0) & 1;
          if(min_exc ==(((L+l2+l2p) >> 0) & 1) &&
             ((abs(l1-l2p)<=k) && (k<=l1+l2p)) &&
             ((abs(l2-l1p)<=k) && (k<=l2+l1p))) {
            sum_k += pow(-1,min_exc)*F_exc*wigner_3j0(l1,k,l2p)*wigner_3j0(l2,k,l1p)
                    *wigner_6j(l1p,l2p,L,l1,l2,k);
          }
        }
        v_mat = pow(-1,(l1+l2))*Y_norm*sum_k;
      }
    }
  }

  return 0;
}
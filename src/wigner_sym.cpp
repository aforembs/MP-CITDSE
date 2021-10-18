#include "wigner_sym.h"

uint fact(uint x) {
  uint tab[] = {1,1,2,6,24,120,720};

  if (x<=6)
    return tab[x];
  else {
    uint res = 720;
    for(uint i = 7; i <= x; ++i)
      res *= i;
    return res;
  }
}

double triangle(uint a, uint b, uint c) {
  double res = sqrt(fact(a-b+c)*fact(a+b-c)*fact(-a+b+c)/(double)fact(a+c+b+1));
  return res;
}

double wigner_3j0(uint j1, uint j2, uint j3) {
  uint t_half = (j1+j2+j3)/2;

  return pow(-1,-t_half)*fact(t_half)*triangle(j1,j3,j2)/
         fact(t_half-j1)/fact(t_half-j2)/fact(t_half-j3); 
}

double wigner_6j(uint j1, uint j2, uint j3, 
                uint j4, uint j5, uint j6) {
  uint max_z = std::min(std::min(2*j4,-j2+j4+j6),j4-j5+j3);
  double aa = fact(j2+j4+j6+1)*fact(j4+j5+j3+1)/fact(j2+j4-j6)/fact(j6-j5+j1)/
              fact(j6+j5-j1)/fact(-j1+j2+j3)/fact(j1-j2+j3)/fact(j4+j5-j3); 
  double s0 = 0.0;

  for(uint z=0; z<=max_z; ++z) {
    s0 += pow(-1,z)*fact(2*j4-z)/fact(z)*fact(j4+j6-j1+j3-z)/
          fact(-j2+j4+j6-z)*fact(j4+j6+j1+j3+1-z)/
          fact(j4-j5+j3-z)/fact(j2+j4+j6+1-z)/fact(j4+j5+j3+1-z); 
  }
  return pow(-1,j4+j6+j1+j3)*triangle(j2,j4,j6)*triangle(j6,j5,j1)*
         triangle(j2,j1,j3)*triangle(j4,j5,j3)*aa*s0;
}

double wigner_6j_2e(uint L, uint la, uint lb, uint lc) {
  uint L2 = L+1;
  uint max_z = std::min(std::min(2*la,-L2+la+lc),std::min(la-lb+1,la+lc-L+1));

  double aa = fact(L2+la+lc+1)*fact(la+lb+2)/fact(L2+la-lc)/fact(lc-lb+L)/
              fact(lc+lb-L)/2/fact(la+lb-1);
  double s0 = 0.0;

  for (uint z=0; z<=max_z; ++z) {
    s0 += pow(-1,z)*fact(2*la-z)/fact(z)*fact(la+lc-L+1-z)/fact(-L2+la+lc-z)
          *fact(la+lc+L+2-z)/fact(la-lb+1-z)/fact(L2+la+lc+1-z)/fact(la+lb+2-z);
  }
  return pow(-1,la+lc+L2)*triangle(L2,la,lc)*triangle(L2,L,1)
        *triangle(lc,lb,L)*triangle(la,lb,1)*aa*s0;
}
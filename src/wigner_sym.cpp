#include "wigner_sym.h"

long int fact(int x) {
  constexpr long int tab[] = {1,1,2,6,24,120,720,5040,40320};

  if (x<=8)
    return tab[x];
  else {
    long int res = 40320;
    for(int i = 8; i <= x; ++i)
      res *= i;
    return res;
  }
}

double triangle(int a, int b, int c) {
  double res = sqrt(fact(a-b+c)*fact(a+b-c)*fact(-a+b+c)/(double)fact(a+c+b+1));
  return res;
}

double wigner_3j0(int j1, int j2, int j3) {
  int t_half = (j1+j2+j3)/2;

  return pow(-1,-t_half)*fact(t_half)*triangle(j1,j3,j2)/(double)
         fact(t_half-j1)/(double)fact(t_half-j2)/(double)fact(t_half-j3); 
}

double wigner_6j(int j1, int j2, int j3, 
                int j4, int j5, int j6) {
  int max_z = std::min(std::min(2*j4,-j2+j4+j6),j4-j5+j3);
  double aa = fact(j2+j4+j6+1)*fact(j4+j5+j3+1)/(double)fact(j2+j4-j6)/(double)fact(j6-j5+j1)/
              (double)fact(j6+j5-j1)/(double)fact(-j1+j2+j3)/(double)fact(j1-j2+j3)/(double)fact(j4+j5-j3); 
  double s0 = 0.0;

  for(int z=0; z<=max_z; ++z) {
    s0 += pow(-1,z)*fact(2*j4-z)/(double)fact(z)*fact(j4+j6-j1+j3-z)/
          (double)fact(-j2+j4+j6-z)*fact(j4+j6+j1+j3+1-z)/
          (double)fact(j4-j5+j3-z)/(double)fact(j2+j4+j6+1-z)/(double)fact(j4+j5+j3+1-z); 
  }
  return pow(-1,j4+j6+j1+j3)*triangle(j2,j4,j6)*triangle(j6,j5,j1)*
         triangle(j2,j1,j3)*triangle(j4,j5,j3)*aa*s0;
}

double wigner_6j_2e(int L, int la, int lb, int lc) {
  int L2 = L+1;
  int max_z = std::min(std::min(2*la,-L2+la+lc),std::min(la-lb+1,la+lc-L+1));

  double aa = fact(L2+la+lc+1)*fact(la+lb+2)/(double)fact(L2+la-lc)/(double)fact(lc-lb+L)/
              (double)fact(lc+lb-L)/2.0/(double)fact(la+lb-1);
  double s0 = 0.0;

  for (int z=0; z<=max_z; ++z) {
    s0 += pow(-1,z)*fact(2*la-z)/(double)fact(z)*fact(la+lc-L+1-z)/(double)fact(-L2+la+lc-z)
          *fact(la+lc+L+2-z)/(double)fact(la-lb+1-z)/(double)fact(L2+la+lc+1-z)/(double)fact(la+lb+2-z);
  }
  return pow(-1,la+lc+L2)*triangle(L2,la,lc)*triangle(L2,L,1)
        *triangle(lc,lb,L)*triangle(la,lb,1)*aa*s0;
}
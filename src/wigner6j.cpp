#include "wigner6j.h"

unsigned int fact(unsigned int x) {
  unsigned int tab[] = {1,1,2,6,24,120,720};

  if (x<=6)
    return tab[x];
  else {
    unsigned int res = 720;
    for(unsigned int i = 7; i <= x; ++i)
      res *= i;
    return res;
  }
}

double triangle(unsigned int a, unsigned int b, unsigned int c) {
  double res = sqrt(fact(a-b+c)*fact(a+b-c)*fact(-a+b+c)/(double)fact(a+c+b+1));
  return res;
}

double wigner_6j_2e(unsigned int L, unsigned int la, unsigned int lb, unsigned int lc) {
  unsigned int L2 = L+1;
  unsigned int max_z = std::min(std::min(2*la,-L2+la+lc),std::min(la-lb+1,la+lc-L+1));

  double aa = fact(L2+la+lc+1)*fact(la+lb+2)/fact(L2+la-lc)/fact(lc-lb+L)/fact(lc+lb-L)/2/fact(la+lb-1);
  double s0 = 0.0;

  for (unsigned int z=0; z<=max_z; ++z) {
    s0 += pow(-1,z)*fact(2*la-z)/fact(z)*fact(la+lc-L+1-z)/fact(-L2+la+lc-z)*fact(la+lc+L+2-z)/fact(la-lb+1-z)/fact(L2+la+lc+1-z)/fact(la+lb+2-z);
  }
  return pow(-1,la+lc+L2)*triangle(L2,la,lc)*triangle(L2,L,1)*triangle(lc,lb,L)*triangle(la,lb,1)*aa*s0;
}
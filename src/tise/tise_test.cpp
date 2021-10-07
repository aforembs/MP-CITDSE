#include <vector>
#include "fastgl.h"
#include "bsplines.h"
#include "tise.h"

int main() {
  int n=10000, k=9, r_max=5000;
  double z=2.0, mass=0.5;
  int l_max=5;

  std::vector<double> gl_x(k);
  std::vector<double> gl_w(k);
  fastgl::QuadPair gl_i;
  for(int i=1; i<=k; ++i) {
    gl_i = fastgl::GLPair(k, i); // generate GL nodes and weights over B-splines support
    gl_x[k-i] = gl_i.x(); 
    gl_w[k-i] = gl_i.weight;
  }

  double fk = 0.125;
  std::vector<double> kkn;
  std::vector<double> spl;
  std::vector<double> splp;
  bsp::GenKnots(n, k, r_max, fk, 'l', kkn);

  bsp::Splines(n, k, gl_x, kkn, spl);
  bsp::SplinesP(n, k, gl_x, kkn, splp);

  tise::GenCoeff(n, k, l_max, z, mass, "he", gl_w, gl_x, kkn, spl, splp, "dat/");
  return 0;
}
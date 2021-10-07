#include <unistd.h>
#include <cstdlib>
#include <filesystem>
#include "fastgl.h"
#include "bsplines.h"
#include "tise.h"

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
  std::string opt_file;
  auto n=400;
  auto k=9;
  auto r_max=200;
  auto l_max=4;
  auto z=1;
  std::string grid, k_file, pot;
  double mass=0.5, fkn=0.125;

  auto kkn = std::vector<double>();

  for(;;) {
    switch(getopt(argc, argv, "h")) {
      case 'h':
        return -1;
      case 'f':
        opt_file = optarg; 
        continue;
    }
    break;
  }

  tise::ReadConfig(opt_file, n, k, r_max, grid, k_file, pot, l_max, z, mass);

  if(grid.compare("linear")==0) {
    bsp::GenKnots(n, k, r_max, fkn, 'l', kkn);
  } else if(grid.compare("exponential")==0) {
    bsp::GenKnots(n, k, r_max, fkn, 'e', kkn);
  } else if(grid.compare("sine")==0) {
    bsp::GenKnots(n, k, r_max, fkn, 's', kkn);
  } else if(grid.compare("custom")==0) {
    if(std::to_string(fs::path(k_file).extension()).compare(".txt")==0) {
      bsp::GenKnots(n, k, r_max, k_file, 't', kkn);
    } else if (std::to_string(fs::path(k_file).extension()).compare(".bin")==0) {
      bsp::GenKnots(n, k, r_max, k_file, 'b', kkn);
    } else {
      std::cout << "Invalid knot file extension, use .txt or .bin\n"; 
      return -1;
    }
  } else {
    std::cout << "Invalid knot type, use 'linear', 'exponential', 'sine' or 'custom'\n"; 
    return -1;
  }

  // generate GL nodes and weights over B-splines support
  auto gl_x = std::vector<double>(k);
  auto gl_w = std::vector<double>(k);
  fastgl::QuadPair gl_i;
  for(int i=1; i<=k; ++i) {
    gl_i = fastgl::GLPair(k, i);
    gl_x[k-i] = gl_i.x(); 
    gl_w[k-i] = gl_i.weight;
  }

  // generate B_i(r) and B'_i(r)
  auto spl  = std::vector<double>();
  auto splp = std::vector<double>();
  bsp::Splines(n, k, gl_x, kkn, spl);
  bsp::SplinesP(n, k, gl_x, kkn, splp);

  tise::GenCoeff(n, k, l_max, z, mass, pot, gl_w, gl_x, kkn, spl, splp, "dat/");

  return 0;
}
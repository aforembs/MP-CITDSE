#include "fastgl.hpp"
#include "tise.hpp"
#include <cstdlib>
#include <filesystem>
#include <unistd.h>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
  std::string opt_file;
  auto n = 400;
  auto k = 9;
  auto glq_pt = 9;
  auto r_max = 200;
  auto l_max = 4;
  auto z = 1;
  std::string grid, k_file, pot;
  double mass = 0.5, fkn = 0.00012;

  auto kkn = std::vector<double>();

  for (;;) {
    switch (getopt(argc, argv, "hf:")) {
    case 'h':
      std::cout << "Program for calculating 1e eigenstates using B-splines\n"
                << "-f <path> yaml input file with the input settings\n";
      return -1;
    case 'f':
      opt_file = optarg;
      continue;
    }
    break;
  }

  tise::ReadConfig(opt_file, n, k, glq_pt, r_max, grid, k_file, pot, l_max, z,
                   mass);

  if (grid.compare("linear") == 0) {
    bsp::GenKnots(n, k, r_max, fkn, 'l', kkn);
  } else if (grid.compare("exponential") == 0) {
    bsp::GenKnots(n, k, r_max, fkn, 'e', kkn);
  } else if (grid.compare("sine") == 0) {
    bsp::GenKnots(n, k, r_max, fkn, 's', kkn);
  } else if (grid.compare("custom") == 0) {
    if (fs::path(k_file).extension().string().compare(".txt") == 0) {
      bsp::GenKnots(n, k, r_max, k_file, 't', kkn);
    } else if (fs::path(k_file).extension().string().compare(".bin") == 0) {
      bsp::GenKnots(n, k, r_max, k_file, 'b', kkn);
    } else {
      std::cout << "Invalid knot file extension, use .txt or .bin\n";
      return -1;
    }
  } else {
    std::cout << "Invalid knot type, use 'linear', 'exponential', 'sine' or "
                 "'custom'\n";
    return -1;
  }

  // generate GL nodes and weights over B-splines support
  auto gl_x = std::vector<double>(glq_pt);
  auto gl_w = std::vector<double>(glq_pt);
  fastgl::QuadPair gl_i;
  for (int i = 1; i <= glq_pt; ++i) {
    gl_i = fastgl::GLPair(glq_pt, i);
    gl_x[glq_pt - i] = gl_i.x();
    gl_w[glq_pt - i] = gl_i.weight;
  }

  // generate B_i(r) and B'_i(r)
  auto spl = std::vector<double>();
  auto splp = std::vector<double>();
  bsp::Splines(n, k, glq_pt, gl_x, kkn, spl);
  bsp::SplinesP(n, k, glq_pt, gl_x, kkn, splp);

  tise::GenCoeff(n, k, glq_pt, l_max, z, mass, pot, gl_w, gl_x, kkn, spl, splp,
                 "dat/");

  return 0;
}
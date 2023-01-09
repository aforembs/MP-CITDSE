#include "w1e.hpp"
#include <cstdlib>
#include <filesystem>
#include <unistd.h>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
  std::string opt_file;
  int qsz, R_max, l_max;
  std::string pot;
  std::string out_prefix;
  std::string quad_type;
  std::string quad_file;
  std::string in_quad_layout;
  std::string pt_file;

  for (;;) {
    switch (getopt(argc, argv, "hf:")) {
    case 'h':
      std::cout << "Program for calculating the inter-electronic interaction\n"
                << "coefficients <n1l1;n2l2|r_12|n'1l'1;n'2l'2>\n"
                << "-f <path> yaml input file with the input settings\n";
      return -1;
    case 'f':
      opt_file = optarg;
      continue;
    }
    break;
  }

  w1e::readConfig(opt_file, qsz, R_max, l_max, pot, quad_type, quad_file,
                  in_quad_layout, pt_file);

  out_prefix = "dat/" + pot;

  std::vector<double> q_x(qsz), q_w(qsz);

  if (quad_type.compare("gauss-legendre") == 0) {
    w1e::genGaussLegendre(qsz, R_max, q_x, q_w);
  } else if (quad_type.compare("user-defined") == 0) {
    if (fs::path(quad_file).extension().string().compare(".txt") == 0 ||
        fs::path(quad_file).extension().string().compare(".dat") == 0) {
      w1e::readQuad(qsz, R_max, quad_file, 't', q_x, q_w);
    } else if (fs::path(quad_file).extension().string().compare(".bin") == 0) {
      w1e::readQuad(qsz, R_max, quad_file, 'b', q_x, q_w);
    } else {
      std::cout << "Invalid quadrature file extension, use .txt/.dat or .bin\n";
      return -1;
    }
  } else {
    std::cout << "Invalid quadrature type!\n";
    return -1;
  }

  std::vector<uint8_t> pq_dx;
  int pti_sz = 0;
  if (in_quad_layout.compare("default") == 0) {
    w1e::defaultPointLayout(qsz, pq_dx, pti_sz);
  } else if (in_quad_layout.compare("user-defined") == 0) {
    if (fs::path(pt_file).extension().string().compare(".txt") == 0 ||
        fs::path(pt_file).extension().string().compare(".dat") == 0) {
      w1e::userPointLayout('t', pt_file, qsz, pq_dx, pti_sz);
    } else if (fs::path(quad_file).extension().string().compare(".bin") == 0) {
      w1e::userPointLayout('b', pt_file, qsz, pq_dx, pti_sz);
    } else {
      std::cout << "Invalid layout file extension, use .txt/.dat or .bin\n";
      return -1;
    }
  } else {
    std::cout << "Invalid inner quadrature layout!\n";
    return -1;
  }

  w1e::genWfn(out_prefix, qsz, pti_sz, l_max, q_x, q_w, pq_dx);

  return 0;
}
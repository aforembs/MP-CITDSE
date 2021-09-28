#include "input_r.h"

int inp::ReadConfig(std::string file,
                    int &n, int &k,
                    std::string &grid,
                    std::string &pot,
                    int &l_max) {

  YAML::Node settings = YAML::LoadFile(file);

  n = settings["Basis_Settings"]["state_no"].as<int>();
  std::cout << "Number of States:     "
            << n << std::endl;
  k = settings["Basis_Settings"]["max_spline_k"].as<int>();
  std::cout << "Maximum B-splines order:  "
            << k << std::endl;
  grid = settings["Basis_Settings"]["grid"].as<std::string>();
  std::cout << "Type of knot spacing:            "
            << grid << std::endl;
  pot = settings["Basis_Settings"]["potential"].as<std::string>();
  std::cout << "Core Potential:            "
            << grid << std::endl;
  l_max     = settings["Basis_Settings"]["l_max"].as<int>();
  std::cout << "Maximum l:                         "
            << l_max << std::endl;

  return 0;
}
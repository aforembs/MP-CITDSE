#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

static std::string dataDir() {
  return std::string(TEST_SOURCE_DIR) + "/dat";
}

#ifdef HAS_YAML_CPP
#include <yaml-cpp/yaml.h>

TEST_CASE("TDSE settings match input yaml", "[ref][tdse][yaml]") {
  auto d = dataDir();

  YAML::Node settings = YAML::LoadFile(d + "/H_tdat/settings.yaml");

  CHECK(settings["Global_Settings"]["potential"].as<std::string>() == "h");
  CHECK(settings["Global_Settings"]["gauge"].as<std::string>() == "v");
  CHECK(settings["Global_Settings"]["Outer_quadrature_size"].as<int>() ==
        600);

  CHECK(settings["Basis_Settings"]["state_no"].as<int>() == 322);
  CHECK(settings["Basis_Settings"]["atomic_no"].as<int>() == 1);
  CHECK(settings["Basis_Settings"]["R_max"].as<int>() == 240);

  CHECK(settings["Propagator_Settings"]["dt"].as<double>() == 0.1);

  CHECK(settings["Field_Parameters"]["shape"].as<std::string>() == "sine");
  CHECK(settings["Field_Parameters"]["w"].as<double>() == 8.5);
}

#else
TEST_CASE("TDSE reference tests require yaml-cpp", "[ref][tdse]") {
  SKIP("yaml-cpp not available — TDSE reference tests skipped");
}
#endif

#ifdef HAS_HDF5
#include <H5Cpp.h>
#endif

TEST_CASE("TDSE coefficient files", "[ref][tdse][text]") {
  auto tdir = dataDir() + "/H_tdat";

  SECTION("h_ct_0.000000.dat: fully populated initial state") {
    std::ifstream in(tdir + "/h_ct_0.000000.dat");
    REQUIRE(in.is_open());

    std::string line;
    std::getline(in, line); // header
    REQUIRE(line == "#Index, Re(c(t)), Im(c(t))");

    int idx;
    double re, im;
    in >> idx >> re >> im;
    CHECK(idx == 0);
    CHECK(re == 1.0);
    CHECK(im == 0.0);

    double pop = 0.0;
    while (in >> idx >> re >> im) {
      pop += re * re + im * im;
    }
    CHECK(pop == Catch::Approx(0.0));
  }

  SECTION("h_ct_0.100000.dat: populated after first step") {
    std::ifstream in(tdir + "/h_ct_0.100000.dat");
    REQUIRE(in.is_open());

    std::string line;
    std::getline(in, line);

    int idx0;
    double re0, im0;
    in >> idx0 >> re0 >> im0;
    CHECK(idx0 == 0);
    CHECK(re0 == Catch::Approx(0.99875).margin(1e-5));
    CHECK(im0 == Catch::Approx(0.0499792).margin(1e-5));
  }
}

TEST_CASE("TDSE population tracking", "[ref][tdse][text]") {
  auto tdir = dataDir() + "/H_tdat";
  std::ifstream in(tdir + "/h_pop10.dat");
  REQUIRE(in.is_open());

  std::string line;
  std::getline(in, line);
  REQUIRE(line == "#time (a.u.), population");

  double t, pop;
  in >> t >> pop;
  CHECK(t == 0.0);
  CHECK(pop == 1.0);
}

TEST_CASE("TDSE field envelope", "[ref][tdse][text]") {
  auto tdir = dataDir() + "/H_tdat";
  std::ifstream in(tdir + "/h_field.dat");
  REQUIRE(in.is_open());

  std::string line;
  std::getline(in, line);
  CHECK(line == "#time (a.u.), field (A(t) or E(t)) (a.u.)");

  double t, f;
  in >> t >> f;
  CHECK(t == 0.1);
  CHECK(f == Catch::Approx(2.85911e-09).margin(1e-14));

  int nrows = 1;
  while (in >> t >> f) {
    CHECK_FALSE(std::isnan(f));
    CHECK_FALSE(std::isinf(f));
    CHECK(std::abs(f) <= 1.0);
    ++nrows;
  }
  CHECK(nrows > 100);
}

TEST_CASE("TDSE PES data", "[ref][tdse][text]") {
  auto tdir = dataDir() + "/H_tdat";
  std::ifstream in(tdir + "/h_pes0.dat");
  REQUIRE(in.is_open());

  double e, prob;
  in >> e >> prob;
  CHECK(e > 0.0);
  CHECK(e == Catch::Approx(0.000490377).margin(1e-8));
  CHECK(prob > 0.0);
}

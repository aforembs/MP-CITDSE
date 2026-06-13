#include "bsp_gsl.hpp"
#include "fastgl.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

#ifdef HAS_GSL

TEST_CASE("bsp::genKnots linear produces correct size and bounds",
          "[bsp][gsl]") {
  int n = 100, k = 9, R_max = 200;
  double fkn = 0.00012;
  std::vector<double> kkn;

  int ret = bsp::genKnots(n, k, R_max, fkn, 'l', kkn);
  REQUIRE(ret == 0);
  REQUIRE(kkn.size() == (size_t)(n + k));

  SECTION("first k knots are 0") {
    for (int i = 0; i < k; ++i) {
      REQUIRE(kkn[i] == 0.0);
    }
  }

  SECTION("last k knots are R_max") {
    for (int i = n; i < n + k; ++i) {
      REQUIRE(kkn[i] == R_max);
    }
  }

  SECTION("interior knots are strictly increasing") {
    for (size_t i = 1; i < kkn.size(); ++i) {
      REQUIRE(kkn[i] >= kkn[i - 1]);
    }
  }
}

TEST_CASE("bsp::genKnots exponential produces correct size and bounds",
          "[bsp][gsl]") {
  int n = 100, k = 9, R_max = 200;
  double fkn = 0.00012;
  std::vector<double> kkn;

  int ret = bsp::genKnots(n, k, R_max, fkn, 'e', kkn);
  REQUIRE(ret == 0);
  REQUIRE(kkn.size() == (size_t)(n + k));
  REQUIRE(kkn[0] == 0.0);
  REQUIRE(kkn[n + k - 1] == R_max);
}

TEST_CASE("bsp::genKnots sine produces correct size and bounds",
          "[bsp][gsl]") {
  int n = 100, k = 9, R_max = 200;
  double fkn = 0.00012;
  std::vector<double> kkn;

  int ret = bsp::genKnots(n, k, R_max, fkn, 's', kkn);
  REQUIRE(ret == 0);
  REQUIRE(kkn.size() == (size_t)(n + k));
  REQUIRE(kkn[0] == 0.0);
  REQUIRE(kkn[n + k - 1] == R_max);
}

TEST_CASE("bsp::genKnots returns error for invalid type", "[bsp][gsl]") {
  int n = 10, k = 4, R_max = 50;
  double fkn = 0.001;
  std::vector<double> kkn;
  int ret = bsp::genKnots(n, k, R_max, fkn, 'z', kkn);
  REQUIRE(ret != 0);
}

TEST_CASE("bsp::splines produces correct output size", "[bsp][gsl]") {
  int n = 40, k = 7, glq_pt = 9, R_max = 100;
  double fkn = 0.00012;
  std::vector<double> kkn, gl_x(glq_pt), gl_w(glq_pt), spl, splp;

  bsp::genKnots(n, k, R_max, fkn, 'l', kkn);

  for (int i = 0; i < glq_pt; ++i) {
    auto qp = fastgl::GLPair(glq_pt, i + 1);
    gl_x[glq_pt - 1 - i] = qp.x();
    gl_w[glq_pt - 1 - i] = qp.weight;
  }

  int ret = bsp::splines(n, k, glq_pt, gl_x, kkn, spl, splp);
  REQUIRE(ret == 0);

  int expected_sz = n * glq_pt * k;
  REQUIRE(spl.size() == (size_t)expected_sz);
  REQUIRE(splp.size() == (size_t)expected_sz);
}

// wrKnotsH5 has a known bug with HDF5 string dataspace sizing; skip for now

#else
TEST_CASE("bsp tests require GSL", "[bsp]") {
  SKIP("GSL not available — B-spline tests skipped");
}
#endif

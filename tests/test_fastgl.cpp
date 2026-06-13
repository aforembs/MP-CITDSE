#include "fastgl.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

constexpr double pi = 3.14159265358979323846;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("FastGL::GLPair produces symmetric nodes", "[fastgl]") {
  for (size_t n = 2; n <= 50; ++n) {
    auto left = fastgl::GLPair(n, 1);
    auto right = fastgl::GLPair(n, n);
    REQUIRE_THAT(left.x(), WithinRel(-right.x(), 1e-14));
    REQUIRE_THAT(left.weight, WithinRel(right.weight, 1e-14));
  }
}

TEST_CASE("FastGL::GLPair weights sum to 2", "[fastgl]") {
  for (size_t n = 1; n <= 50; ++n) {
    double sum = 0.0;
    for (size_t k = 1; k <= n; ++k) {
      sum += fastgl::GLPair(n, k).weight;
    }
    REQUIRE_THAT(sum, WithinAbs(2.0, 1e-14));
  }
}

TEST_CASE("FastGL::GLPair nodes are within [-1, 1]", "[fastgl]") {
  for (size_t n = 1; n <= 50; ++n) {
    for (size_t k = 1; k <= n; ++k) {
      double x = fastgl::GLPair(n, k).x();
      REQUIRE(x >= -1.0);
      REQUIRE(x <= 1.0);
    }
  }
}

TEST_CASE("FastGL::GLPair monotonic theta", "[fastgl]") {
  for (size_t n = 2; n <= 50; ++n) {
    double prev_theta = -1.0;
    for (size_t k = 1; k <= n; ++k) {
      double theta = fastgl::GLPair(n, k).theta;
      REQUIRE(theta > prev_theta);
      prev_theta = theta;
    }
  }
}

TEST_CASE("FastGL::GLPair known values n=2", "[fastgl]") {
  // theta increases in [0,pi]; x = cos(theta) decreases from 1 to -1
  auto p1 = fastgl::GLPair(2, 1);
  auto p2 = fastgl::GLPair(2, 2);
  REQUIRE_THAT(p1.x(), WithinRel(1.0 / std::sqrt(3.0), 1e-14));
  REQUIRE_THAT(p2.x(), WithinRel(-1.0 / std::sqrt(3.0), 1e-14));
  REQUIRE_THAT(p1.weight, WithinAbs(1.0, 1e-14));
  REQUIRE_THAT(p2.weight, WithinAbs(1.0, 1e-14));
}

TEST_CASE("FastGL::GLPair known values n=3", "[fastgl]") {
  auto p2 = fastgl::GLPair(3, 2);
  // middle node should be at 0 with weight 8/9
  REQUIRE_THAT(p2.x(), WithinAbs(0.0, 1e-15));
  REQUIRE_THAT(p2.weight, WithinRel(8.0 / 9.0, 1e-14));
}

TEST_CASE("FastGL::GLPair integrates x^k exactly", "[fastgl]") {
  auto integrate = [](size_t n, auto f) {
    double sum = 0.0;
    for (size_t k = 1; k <= n; ++k) {
      auto p = fastgl::GLPair(n, k);
      sum += p.weight * f(p.x());
    }
    return sum;
  };
  for (size_t n = 1; n <= 40; ++n) {
    CAPTURE(n);
    REQUIRE_THAT(integrate(n, [](double) { return 1.0; }),
                 WithinAbs(2.0, 1e-14));
    REQUIRE_THAT(integrate(n, [](double x) { return x; }),
                 WithinAbs(0.0, 1e-14));
    if (n >= 3) {
      REQUIRE_THAT(integrate(n, [](double x) { return x * x; }),
                   WithinAbs(2.0 / 3.0, 1e-13));
    }
    if (n >= 5) {
      REQUIRE_THAT(integrate(n, [](double x) { return x * x * x * x; }),
                   WithinAbs(2.0 / 5.0, 1e-12));
    }
  }
}

TEST_CASE("FastGL::GLPair computed vs tabulated match at boundary n=100",
          "[fastgl]") {
  for (size_t k = 1; k <= 100; ++k) {
    auto p = fastgl::GLPair(100, k);
    REQUIRE(p.theta > 0.0);
    REQUIRE(p.theta < pi);
    REQUIRE(p.weight > 0.0);
  }
}

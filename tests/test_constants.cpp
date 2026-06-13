#include "au.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("au:: fundamental constants are correct", "[constants]") {
  SECTION("electron mass is 1 a.u.") {
    REQUIRE(au::m_e_ == 1.0);
  }
  SECTION("electron charge is 1 a.u.") {
    REQUIRE(au::e_ == 1.0);
  }
  SECTION("reduced Planck constant is 1 a.u.") {
    REQUIRE(au::h_ == 1.0);
  }
  SECTION("speed of light is ~137 a.u.") {
    REQUIRE_THAT(au::c_, WithinRel(137.035989561, 1e-12));
  }
  SECTION("fine structure constant is 1/c") {
    REQUIRE_THAT(au::a_, WithinRel(1.0 / au::c_, 1e-15));
  }
  SECTION("Rydberg energy is 0.5 a.u.") {
    REQUIRE(au::Ryd_ == 0.5);
  }
  SECTION("ground state hydrogen energy is -0.5 a.u.") {
    REQUIRE(au::GEn_H_ == -0.5);
  }
  SECTION("pi is accurate") {
    REQUIRE_THAT(au::m_pi_, WithinRel(3.14159265358979323846, 1e-15));
  }
}

TEST_CASE("conv:: conversion factors are internally consistent", "[constants]") {
  SECTION("Rydberg in a.u. is 0.5") {
    REQUIRE(conv::En_Ryd_au_ == 0.5);
  }
  SECTION("eV to a.u. conversion is self-consistent") {
    double ev_to_au = conv::En_ev_au_;
    double au_to_ev = conv::enau;
    REQUIRE_THAT(1.0 / ev_to_au, WithinRel(au_to_ev, 1e-10));
  }
  SECTION("time conversions are positive") {
    REQUIRE(conv::t_s_au_ > 0.0);
    REQUIRE(conv::t_au_fs_ > 0.0);
  }
}

TEST_CASE("au:: mass ratios are physically reasonable", "[constants]") {
  SECTION("proton/electron mass ratio ~1836") {
    REQUIRE_THAT(au::m_p_ / au::m_e_, WithinRel(1836.15270137, 1e-12));
  }
  SECTION("neutron/electron mass ratio ~1838") {
    REQUIRE_THAT(au::m_n_ / au::m_e_, WithinRel(1838.68366240, 1e-12));
  }
  SECTION("neutron is heavier than proton") {
    REQUIRE(au::m_n_ > au::m_p_);
  }
}

TEST_CASE("conv:: intensity conversion is finite and positive", "[constants]") {
  REQUIRE(conv::I_W_cm2_au_ > 0.0);
  REQUIRE(std::isfinite(conv::I_W_cm2_au_));
}

TEST_CASE("conv:: energy conversion factors are positive", "[constants]") {
  REQUIRE(conv::En_ev_au_ > 0.0);
  REQUIRE(conv::enau > 0.0);
  REQUIRE(conv::ev_J_ > 0.0);
}

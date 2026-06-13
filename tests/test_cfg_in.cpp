#include "cfg_in.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>

using Catch::Matchers::WithinRel;

static std::string inpDir() {
  std::string dir = std::filesystem::absolute("inp").string();
  if (std::filesystem::exists(dir + "/cfg-0.inp"))
    return dir;
  dir = std::filesystem::absolute("../inp").string();
  if (std::filesystem::exists(dir + "/cfg-0.inp"))
    return dir;
  return std::string(TEST_SOURCE_DIR) + "/inp";
}

TEST_CASE("cfg::readCfg parses L=0 config correctly", "[cfg_in]") {
  std::string dir = inpDir();
  REQUIRE(std::filesystem::exists(dir + "/cfg-0.inp"));

  int sym = 0, ncf = 0;
  std::vector<cfg::line> cfgs;

  int ret = cfg::readCfg(dir, 0, sym, ncf, cfgs);
  REQUIRE(ret == 0);
  REQUIRE(sym == 1);
  REQUIRE(ncf == 15);
  REQUIRE(cfgs.size() == 15);
}

TEST_CASE("cfg::readCfg first config line matches expected", "[cfg_in]") {
  std::string dir = inpDir();

  int sym, ncf;
  std::vector<cfg::line> cfgs;
  cfg::readCfg(dir, 0, sym, ncf, cfgs);

  auto const& first = cfgs[0];
  REQUIRE(first.n1 == 0);
  REQUIRE(first.l1 == 0);
  REQUIRE(first.l2 == 0);
  REQUIRE(first.n2min == 0);
  REQUIRE(first.n2max == 290);
}

TEST_CASE("cfg::readCfg parses L=1 config", "[cfg_in]") {
  std::string dir = inpDir();
  REQUIRE(std::filesystem::exists(dir + "/cfg-1.inp"));

  int sym, ncf;
  std::vector<cfg::line> cfgs;
  int ret = cfg::readCfg(dir, 1, sym, ncf, cfgs);
  REQUIRE(ret == 0);
  REQUIRE(ncf > 0);
}

TEST_CASE("cfg::readCfg returns -1 for missing file", "[cfg_in]") {
  std::string dir = inpDir();

  int sym, ncf;
  std::vector<cfg::line> cfgs;
  int ret = cfg::readCfg(dir, 99, sym, ncf, cfgs);
  REQUIRE(ret == -1);
}

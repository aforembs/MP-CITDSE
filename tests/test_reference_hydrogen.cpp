#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <vector>

#ifdef HAS_HDF5
#include <H5Cpp.h>

static std::string dataDir() {
  return std::string(TEST_SOURCE_DIR) + "/dat";
}

static std::string inpDir() {
  return std::string(TEST_SOURCE_DIR) + "/inp";
}

template <typename T>
static T readAttr(const H5::H5File& f, const std::string& name) {
  H5::Attribute a = f.openAttribute(name);
  T val;
  a.read(a.getDataType(), &val);
  return val;
}

template <typename T>
static std::vector<T> readDataset(const H5::H5File& f,
                                  const std::string& name) {
  H5::DataSet ds = f.openDataSet(name);
  H5::DataSpace sp = ds.getSpace();
  hsize_t dims[2];
  int nd = sp.getSimpleExtentDims(dims);
  hsize_t n = (nd == 1) ? dims[0] : dims[0] * dims[1];
  std::vector<T> buf(n);
  ds.read(buf.data(), ds.getDataType());
  return buf;
}

static std::vector<char> readFileBytes(const std::string& path) {
  std::ifstream in(path, std::ios::binary | std::ios::ate);
  auto sz = in.tellg();
  std::vector<char> buf(sz);
  in.seekg(0);
  in.read(buf.data(), sz);
  return buf;
}

TEST_CASE("H 1e dataset sizes and attributes", "[ref][h1e][hdf5]") {
  auto d = dataDir();

  for (int l = 0; l <= 5; ++l) {
    INFO("l = " << l);
    auto path = d + "/h" + std::to_string(l) + ".h5";
    REQUIRE(std::filesystem::exists(path));

    H5::H5File f(path, H5F_ACC_RDONLY);

    CHECK(readAttr<double>(f, "Z") == 1.0);
    CHECK(readAttr<double>(f, "M") == 0.5);
    CHECK(readAttr<int>(f, "N") == 322);
    CHECK(readAttr<int>(f, "K") == 9);
    CHECK(readAttr<double>(f, "R") == 240.0);
    CHECK(readAttr<int>(f, "l") == l);

    auto en = readDataset<double>(f, "En");
    REQUIRE(en.size() == 320);

    auto coeff = readDataset<double>(f, "Coeff");
    REQUIRE(coeff.size() == 320 * 322);

    auto knots = readDataset<double>(f, "Knots");
    REQUIRE(knots.size() == 331);
  }
}

TEST_CASE("H 1e energies are hydrogenic", "[ref][h1e][hdf5]") {
  auto d = dataDir();

  SECTION("h0: l=0 ground state -0.5") {
    H5::H5File f(d + "/h0.h5", H5F_ACC_RDONLY);
    auto en = readDataset<double>(f, "En");
    CHECK(en[0] == Catch::Approx(-0.5));
    CHECK(en[1] == Catch::Approx(-0.125));
    CHECK(en[2] == Catch::Approx(-0.5 / 9.0));
  }

  SECTION("h1: l=1 first p-state -0.125") {
    H5::H5File f(d + "/h1.h5", H5F_ACC_RDONLY);
    auto en = readDataset<double>(f, "En");
    CHECK(en[0] == Catch::Approx(-0.125));
  }

  SECTION("h2: l=2 first d-state -0.5/9") {
    H5::H5File f(d + "/h2.h5", H5F_ACC_RDONLY);
    auto en = readDataset<double>(f, "En");
    CHECK(en[0] == Catch::Approx(-0.5 / 9.0));
  }

  SECTION("h3: l=3 first f-state -0.5/16") {
    H5::H5File f(d + "/h3.h5", H5F_ACC_RDONLY);
    auto en = readDataset<double>(f, "En");
    CHECK(en[0] == Catch::Approx(-0.5 / 16.0));
  }
}

TEST_CASE("H dipole matrix dimensions", "[ref][dipole][hdf5]") {
  auto d = dataDir();

  for (int l = 0; l < 5; ++l) {
    INFO("l = " << l << " -> l+1");
    auto path =
        d + "/h" + std::to_string(l) + std::to_string(l + 1) + "v.h5";
    REQUIRE(std::filesystem::exists(path));

    H5::H5File f(path, H5F_ACC_RDONLY);
    auto dip = readDataset<double>(f, "d_if");
    CHECK(dip.size() == 320 * 320);
  }
}

TEST_CASE("H config files match inp originals", "[ref][cfg][hdf5]") {
  auto d = dataDir();
  auto inp = inpDir();

  for (int L = 0; L <= 3; ++L) {
    INFO("L = " << L);
    auto datPath =
        d + "/cfg-" + std::to_string(L) + ".inp";
    auto inpPath =
        inp + "/cfg-" + std::to_string(L) + ".inp";

    REQUIRE(std::filesystem::exists(datPath));
    REQUIRE(std::filesystem::exists(inpPath));

    auto datBytes = readFileBytes(datPath);
    auto inpBytes = readFileBytes(inpPath);

    CHECK(datBytes.size() == inpBytes.size());
    if (datBytes.size() == inpBytes.size()) {
      CHECK(std::equal(datBytes.begin(), datBytes.end(), inpBytes.begin()));
    }
  }
}

#else
TEST_CASE("H reference tests require HDF5", "[ref][h1e]") {
  SKIP("HDF5 not available — H reference tests skipped");
}
#endif

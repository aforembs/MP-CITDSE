#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>
#include <vector>

#ifdef HAS_HDF5
#include <H5Cpp.h>

static std::string dataDir() {
  return std::string(TEST_SOURCE_DIR) + "/dat";
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

TEST_CASE("He 1e dataset sizes and attributes", "[ref][he1e][hdf5]") {
  auto d = dataDir();

  for (int l = 0; l <= 7; ++l) {
    INFO("l = " << l);
    auto path = d + "/he" + std::to_string(l) + ".h5";
    REQUIRE(std::filesystem::exists(path));

    H5::H5File f(path, H5F_ACC_RDONLY);

    CHECK(readAttr<double>(f, "Z") == 2.0);
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

TEST_CASE("He 1e energies plausible", "[ref][he1e][hdf5]") {
  auto d = dataDir();

  auto scale = [](int l) {
    double Z = 2.0;
    int n = l + 1;
    return -0.5 * Z * Z / (n * n);
  };

  for (int l = 0; l <= 7; ++l) {
    INFO("l = " << l);
    H5::H5File f(d + "/he" + std::to_string(l) + ".h5", H5F_ACC_RDONLY);
    auto en = readDataset<double>(f, "En");

    CHECK(en[0] == Catch::Approx(scale(l)).margin(1e-8));

    for (auto e : en) {
      REQUIRE_FALSE(std::isnan(e));
      REQUIRE_FALSE(std::isinf(e));
    }
  }
}

TEST_CASE("He CI eigenstates", "[ref][heCI][hdf5]") {
  auto d = dataDir();

  // expected CI basis sizes for L=0..3
  std::vector<int> expected_dims = {3030, 2799, 2808, 2808};

  for (int L = 0; L <= 3; ++L) {
    INFO("L = " << L);
    auto path = d + "/heCI" + std::to_string(L) + ".h5";
    REQUIRE(std::filesystem::exists(path));

    H5::H5File f(path, H5F_ACC_RDONLY);

    auto en = readDataset<double>(f, "En_CI");
    auto vecs = readDataset<double>(f, "CI_vecs");

    CHECK(en.size() == (size_t)expected_dims[L]);
    CHECK(vecs.size() == (size_t)expected_dims[L] * expected_dims[L]);
  }
}

TEST_CASE("He CI ground state energy", "[ref][heCI][hdf5]") {
  auto d = dataDir();
  H5::H5File f(d + "/heCI0.h5", H5F_ACC_RDONLY);
  auto en = readDataset<double>(f, "En_CI");

  CHECK(en[0] == Catch::Approx(-2.88675).margin(5e-5));
}

TEST_CASE("He CI dipole dimensions", "[ref][heCI_dip][hdf5]") {
  auto d = dataDir();
  std::vector<int> dims = {3030, 2799, 2808, 2808};

  for (int L = 0; L <= 2; ++L) {
    INFO("L = " << L);
    auto path =
        d + "/heCI_" + std::to_string(L) + std::to_string(L + 1) + "v.h5";
    REQUIRE(std::filesystem::exists(path));

    H5::H5File f(path, H5F_ACC_RDONLY);
    auto dip = readDataset<double>(f, "CI_dip");
    CHECK(dip.size() == (size_t)dims[L] * dims[L + 1]);
  }
}

TEST_CASE("He 2e intermediate files", "[ref][he2e][hdf5]") {
  auto d = dataDir();

  SECTION("he2_0En.h5") {
    H5::H5File f(d + "/he2_0En.h5", H5F_ACC_RDONLY);
    auto e_2e = readDataset<double>(f, "e_2e");
    auto idx = readDataset<int>(f, "idx");
    CHECK(e_2e.size() == 3030);
    CHECK(idx.size() == 12120);
  }

  SECTION("he2_01v.h5") {
    H5::H5File f(d + "/he2_01v.h5", H5F_ACC_RDONLY);
    auto dip = readDataset<double>(f, "d_if");
    CHECK(dip.size() == 3030 * 2799);
  }

  SECTION("he2.h5 matches 1e structure") {
    H5::H5File f(d + "/he2.h5", H5F_ACC_RDONLY);
    CHECK(readAttr<int>(f, "l") == 2);
    auto en = readDataset<double>(f, "En");
    CHECK(en.size() == 320);
  }
}

TEST_CASE("He V12 matrices", "[ref][heV12][hdf5]") {
  auto d = dataDir();
  std::vector<hsize_t> expected_sizes = {4591965, 3918600, 3943836, 3943836};

  for (int L = 0; L <= 3; ++L) {
    INFO("L = " << L);
    auto path = d + "/heV12_" + std::to_string(L) + ".h5";
    REQUIRE(std::filesystem::exists(path));

    H5::H5File f(path, H5F_ACC_RDONLY);
    auto v12 = readDataset<double>(f, "V_12");
    CHECK(v12.size() == expected_sizes[L]);
  }
}

#else
TEST_CASE("He reference tests require HDF5", "[ref][he]") {
  SKIP("HDF5 not available — He reference tests skipped");
}
#endif

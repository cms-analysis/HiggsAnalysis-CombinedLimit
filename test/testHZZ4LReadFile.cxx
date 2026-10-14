#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <limits>

// These implementation helpers are not part of the public PDF header.
namespace RooFit {
  void readFile();
  double HiggsWidth(int id, double mass);
  double HiggsWidthSM4(int id, double mass);
  // Match the existing mutable global arrays; the preservation test modifies BR.
  // These declarations must retain the production arrays' types and mutability.
  // NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays,cppcoreguidelines-avoid-non-const-global-variables)
  extern double BR[26][217];
  extern double BR_SM4[18][102];
  // NOLINTEND(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays,cppcoreguidelines-avoid-non-const-global-variables)
}  // namespace RooFit

// FNV-1a over IEEE-754 values, least significant byte first on every platform.
// Reference checksums come from the original readFile assignments, including
// the implicit zeroes in SM channel 12, independently of the extracted files.
template <typename Table>
static std::uint64_t checksum(const Table& table) {
  static_assert(sizeof(double) == sizeof(std::uint64_t) && std::numeric_limits<double>::is_iec559);
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (const auto& channel : table) {
    for (double value : channel) {
      std::uint64_t bits = 0;
      std::memcpy(&bits, &value, sizeof(bits));
      for (unsigned int byte = 0; byte < 8; ++byte) {
        hash ^= (bits >> (8 * byte)) & 0xff;
        hash *= UINT64_C(1099511628211);
      }
    }
  }
  return hash;
}

namespace {

  TEST(HZZ4LReadFile, PreservesAllTableValues) {
    RooFit::readFile();
    EXPECT_EQ(checksum(RooFit::BR), UINT64_C(0x478ed166e4d5cf79));
    EXPECT_EQ(checksum(RooFit::BR_SM4), UINT64_C(0x500c4ee8b13e0f8a));
    // Calling the initializer again must leave both tables unchanged.
    RooFit::readFile();
    EXPECT_EQ(checksum(RooFit::BR), UINT64_C(0x478ed166e4d5cf79));
    EXPECT_EQ(checksum(RooFit::BR_SM4), UINT64_C(0x500c4ee8b13e0f8a));
  }

  TEST(HZZ4LReadFile, WidthsAndInterpolation) {
    RooFit::readFile();
    EXPECT_DOUBLE_EQ(RooFit::HiggsWidth(0, 90.), 0.0022);
    EXPECT_DOUBLE_EQ(RooFit::HiggsWidth(1, 90.), 0.0022 * 0.812);
    EXPECT_DOUBLE_EQ(RooFit::HiggsWidth(0, 92.5), (0.0022 + 0.00232) / 2.);
    EXPECT_DOUBLE_EQ(RooFit::HiggsWidth(1, 92.5), ((0.0022 * 0.812) + (0.00232 * 0.803)) / 2.);
    EXPECT_DOUBLE_EQ(RooFit::HiggsWidth(12, 125.), 0.);

    EXPECT_DOUBLE_EQ(RooFit::HiggsWidthSM4(0, 100.), 0.00552);
    EXPECT_DOUBLE_EQ(RooFit::HiggsWidthSM4(1, 100.), 0.00552 * 0.570);
    EXPECT_DOUBLE_EQ(RooFit::HiggsWidthSM4(0, 105.), (0.00552 + 0.00641) / 2.);
    EXPECT_DOUBLE_EQ(RooFit::HiggsWidthSM4(1, 105.), ((0.00552 * 0.570) + (0.00641 * 0.530)) / 2.);
  }

  TEST(HZZ4LReadFile, PreservesUnassignedSMChannel) {
    std::array<double, 217> saved{};
    unsigned int mass = 0;
    for (auto& value : RooFit::BR[12]) {
      saved.at(mass) = value;
      value = mass + 0.5;
      ++mass;
    }

    RooFit::readFile();
    mass = 0;
    for (auto& value : RooFit::BR[12]) {
      EXPECT_DOUBLE_EQ(value, mass + 0.5);
      value = saved.at(mass);
      ++mass;
    }
  }

  TEST(HZZ4LReadFile, RejectsInvalidRequests) {
    RooFit::readFile();
    EXPECT_EQ(RooFit::HiggsWidth(-1, 125.), -1.);
    EXPECT_EQ(RooFit::HiggsWidth(26, 125.), -1.);
    EXPECT_EQ(RooFit::HiggsWidth(0, 89.), -1.);
    EXPECT_EQ(RooFit::HiggsWidth(0, 1001.), -1.);
    EXPECT_EQ(RooFit::HiggsWidthSM4(-1, 125.), 0.);
    EXPECT_EQ(RooFit::HiggsWidthSM4(18, 125.), 0.);
    EXPECT_EQ(RooFit::HiggsWidthSM4(0, 99.), 0.);
    EXPECT_EQ(RooFit::HiggsWidthSM4(0, 1001.), 0.);
  }

}  // namespace

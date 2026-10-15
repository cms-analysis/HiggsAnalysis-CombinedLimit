#include "../interface/HZZ4LRooPdfs.h"
#include "RooRealVar.h"

#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>

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

namespace {

  template <typename Base, bool HasInterference>
  class HighMassProbe : public Base {
  public:
    using Base::Base;
    static constexpr std::size_t parameterCount = HasInterference ? 5 : 2;

    void initialize(bool eightTeV) {
      this->is8TeV = eightTeV;
      this->initMatrices();
    }

    auto matrixChecksums() const {
      if constexpr (HasInterference) {
        return std::array<std::uint64_t, 5>{{checksum(this->a_width),
                                             checksum(this->a_delta),
                                             checksum(this->a_alpha),
                                             checksum(this->a_r),
                                             checksum(this->a_beta)}};
      } else {
        return std::array<std::uint64_t, 2>{{checksum(this->a_width), checksum(this->a_delta)}};
      }
    }

    auto parameters(double mass, double coupling) const {
      if constexpr (HasInterference) {
        return std::array<double, 5>{{this->getWidth(mass, coupling),
                                      this->getDelta(mass, coupling),
                                      this->getAlpha(mass, coupling),
                                      this->getR(mass, coupling),
                                      this->getBeta(mass, coupling)}};
      } else {
        return std::array<double, 2>{{this->getWidth(mass, coupling), this->getDelta(mass, coupling)}};
      }
    }

    void overwriteMatrices() {
      const auto overwrite = [](auto& matrix) {
        for (auto& row : matrix) {
          for (auto& value : row)
            value = -999.;
        }
      };
      overwrite(this->a_width);
      overwrite(this->a_delta);
      if constexpr (HasInterference) {
        overwrite(this->a_alpha);
        overwrite(this->a_r);
        overwrite(this->a_beta);
      }
    }
  };

  class WidthProbe : public RooBWHighMassGGH {
  public:
    using RooBWHighMassGGH::H_width;
  };

}  // namespace

// Reference values below were captured from the pre-extraction library, not
// calculated from the new .inc files. Matrix order: width, delta, alpha, r, beta.
template <typename Probe>
static void checkHighMassMatrices(
    Probe& model,
    const std::array<std::array<std::uint64_t, Probe::parameterCount>, 2>& expectedHashes,
    const std::array<std::array<std::array<double, Probe::parameterCount>, 2>, 2>& expectedParameters) {
  constexpr std::array<std::array<double, 2>, 2> points = {{{{450., 0.3}}, {{775., 0.65}}}};
  EXPECT_EQ(model.matrixChecksums(), expectedHashes.at(0));  // Constructor initialization.
  for (std::size_t energy = 0; energy < 2; ++energy) {
    SCOPED_TRACE(energy);
    model.initialize(energy != 0);
    EXPECT_EQ(model.matrixChecksums(), expectedHashes.at(energy));
    const Probe copied(model);
    EXPECT_EQ(copied.matrixChecksums(), expectedHashes.at(energy));
    for (std::size_t point = 0; point < points.size(); ++point) {
      const auto actual = model.parameters(points.at(point).at(0), points.at(point).at(1));
      for (std::size_t parameter = 0; parameter < actual.size(); ++parameter) {
        SCOPED_TRACE(::testing::Message() << "point=" << point << ", parameter=" << parameter);
        EXPECT_DOUBLE_EQ(actual.at(parameter), expectedParameters.at(energy).at(point).at(parameter));
      }
    }
    model.overwriteMatrices();
    model.initialize(energy != 0);
    EXPECT_EQ(model.matrixChecksums(), expectedHashes.at(energy));
  }
  model.initialize(/*eightTeV=*/false);
  EXPECT_EQ(model.matrixChecksums(), expectedHashes.at(0));
}

template <typename Pdf, std::size_t Points>
static void checkSpline(const Pdf& model, const std::array<std::array<double, 2>, Points>& reference) {
  for (const auto& point : reference) {
    SCOPED_TRACE(point.at(0));
    EXPECT_DOUBLE_EQ(model.Spline(point.at(0)), point.at(1));
  }
}

namespace {
  class HZZ4LTables : public ::testing::Test {
  protected:
    RooRealVar m_x{/*name=*/"x", /*title=*/"x", /*value=*/600.};
    RooRealVar m_mass{/*name=*/"mass", /*title=*/"mass", /*value=*/600.};
    RooRealVar m_coupling{/*name=*/"coupling", /*title=*/"coupling", /*value=*/0.5};
    RooRealVar m_br{/*name=*/"br", /*title=*/"br", /*value=*/0.};
    RooRealVar m_strength{/*name=*/"strength", /*title=*/"strength", /*value=*/1.};
    RooRealVar m_scale{/*name=*/"scale", /*title=*/"scale", /*value=*/1.};
  };

  TEST_F(HZZ4LTables, GGHMatrices) {
    HighMassProbe<RooCPSHighMassGGH, true> model(
        "model", "model", m_x, m_mass, m_coupling, m_br, m_strength, m_scale, false);
    const std::array<std::array<std::uint64_t, 5>, 2> expectedHashes = {{
        {{UINT64_C(0xfd839cb8d527d7cc),
          UINT64_C(0x1acedc5bac6bbaf0),
          UINT64_C(0x8c8f4418f90a1098),
          UINT64_C(0x903f8f445a1f61ad),
          UINT64_C(0x2e06f587762a7c91)}},
        {{UINT64_C(0x33ebdc9f1e3a2946),
          UINT64_C(0x1acedc5bac6bbaf0),
          UINT64_C(0xd077419ba7a41c25),
          UINT64_C(0x67317da146161f02),
          UINT64_C(0x2e06f587762a7c91)}},
    }};
    const std::array<std::array<std::array<double, 5>, 2>, 2> expectedParameters = {{
        {{
            {{42.650330543500004, 0.15118292160332497, 1.0520236641170002, 1.546689145469e-05, 5.000000000000001}},
            {{217.17457580568748, 6.90548980236875, 2.042708046735625, 2.1872264426306253e-06, 5}},
        }},
        {{
            {{42.650330543500004, 0.15118292160332497, 0.96787841618, 1.5735916349530006e-05, 5.000000000000001}},
            {{217.16071701056248, 6.90548980236875, 2.283775582913125, 2.163183069115625e-06, 5}},
        }},
    }};
    checkHighMassMatrices(model, expectedHashes, expectedParameters);
  }

  TEST_F(HZZ4LTables, BWMatrices) {
    HighMassProbe<RooBWHighMassGGH, true> model(
        /*name=*/"model", /*title=*/"model", m_x, m_mass, m_coupling, m_br, m_strength, false);
    const std::array<std::array<std::uint64_t, 5>, 2> expectedHashes = {{
        {{UINT64_C(0xf23735af37408524),
          UINT64_C(0xa6430033e87a8f32),
          UINT64_C(0xd5668e4a6180f121),
          UINT64_C(0x3698ff7f928dcc7f),
          UINT64_C(0x2e06f587762a7c91)}},
        {{UINT64_C(0x136bba8d1b349962),
          UINT64_C(0x808d7241c900329f),
          UINT64_C(0x9cc6f379bedc0993),
          UINT64_C(0xb55b5e2f04a429e2),
          UINT64_C(0x2e06f587762a7c91)}},
    }};
    const std::array<std::array<std::array<double, 5>, 2>, 2> expectedParameters = {{
        {{
            {{48.742579936949994, 0.10911315772679997, 0.8817044496542499, 1.0507130582475002e-05, 5.000000000000001}},
            {{272.6959533690625, 4.554142951963125, 1.065753161906125, 1.6832410736375e-06, 5}},
        }},
        {{
            {{48.682844638825, 0.14059781748804995, 0.68090724945075, 1.0480038440607501e-05, 5.000000000000001}},
            {{254.47477340712499, 7.1606749892199995, 1.035383187236625, 1.62021111549375e-06, 5}},
        }},
    }};
    checkHighMassMatrices(model, expectedHashes, expectedParameters);
  }

  TEST_F(HZZ4LTables, GGHNoInterferenceMatrices) {
    HighMassProbe<RooCPSHighMassGGHNoInterf, false> model(
        /*name=*/"model", /*title=*/"model", m_x, m_mass, m_coupling, m_br, false);
    const std::array<std::array<std::uint64_t, 2>, 2> expectedHashes = {{
        {{UINT64_C(0xfd839cb8d527d7cc), UINT64_C(0x1acedc5bac6bbaf0)}},
        {{UINT64_C(0x33ebdc9f1e3a2946), UINT64_C(0x1acedc5bac6bbaf0)}},
    }};
    const std::array<std::array<std::array<double, 2>, 2>, 2> expectedParameters = {{
        {{
            {{42.650330543500004, 0.15118292160332497}},
            {{217.17457580568748, 6.90548980236875}},
        }},
        {{
            {{42.650330543500004, 0.15118292160332497}},
            {{217.16071701056248, 6.90548980236875}},
        }},
    }};
    checkHighMassMatrices(model, expectedHashes, expectedParameters);
  }

  TEST_F(HZZ4LTables, VBFMatrices) {
    HighMassProbe<RooCPSHighMassVBF, true> model(
        "model", "model", m_x, m_mass, m_coupling, m_br, m_strength, m_scale, false);
    const std::array<std::array<std::uint64_t, 5>, 2> expectedHashes = {{
        {{UINT64_C(0xf78621f296bfab54),
          UINT64_C(0x010155766ac08b8f),
          UINT64_C(0x366227c65dba739a),
          UINT64_C(0xc7cb5a3d5d9bd418),
          UINT64_C(0x2e06f587762a7c91)}},
        {{UINT64_C(0xc2b5209d8a7a8764),
          UINT64_C(0x0ba85d578bd7a97d),
          UINT64_C(0xc8037ac133718995),
          UINT64_C(0xa2e68eb90e89ea68),
          UINT64_C(0x2e06f587762a7c91)}},
    }};
    const std::array<std::array<std::array<double, 5>, 2>, 2> expectedParameters = {{
        {{
            {{42.786666870124996, 0.043665528064497486, 2.007501845247001, 9.559534987605e-06, 5.000000000000001}},
            {{207.50382900268747, 8.653739854693125, 1.0584758129004999, 1.904930812205e-06, 5}},
        }},
        {{
            {{41.847573757175, 0.0500111044384, 3.7646417878571516, 9.412172289557501e-06, 5.000000000000001}},
            {{215.2138080599375, 3.905521914365626, 0.679815314709875, 1.6529832507225e-06, 5}},
        }},
    }};
    checkHighMassMatrices(model, expectedHashes, expectedParameters);
  }

  TEST_F(HZZ4LTables, VBFNoInterferenceMatrices) {
    HighMassProbe<RooCPSHighMassVBFNoInterf, false> model(
        /*name=*/"model", /*title=*/"model", m_x, m_mass, m_coupling, m_br, false);
    const std::array<std::array<std::uint64_t, 2>, 2> expectedHashes = {{
        {{UINT64_C(0xf78621f296bfab54), UINT64_C(0x010155766ac08b8f)}},
        {{UINT64_C(0xc2b5209d8a7a8764), UINT64_C(0x0ba85d578bd7a97d)}},
    }};
    const std::array<std::array<std::array<double, 2>, 2>, 2> expectedParameters = {{
        {{
            {{42.786666870124996, 0.043665528064497486}},
            {{207.50382900268747, 8.653739854693125}},
        }},
        {{
            {{41.847573757175, 0.0500111044384}},
            {{215.2138080599375, 3.905521914365626}},
        }},
    }};
    checkHighMassMatrices(model, expectedHashes, expectedParameters);
  }

  TEST_F(HZZ4LTables, HDECAYWidthsAndIntegerTruncation) {
    const WidthProbe model;
    std::array<std::array<double, 901>, 1> widths{};
    for (std::size_t index = 0; index < widths.at(0).size(); ++index) {
      const double mass = 100. + static_cast<double>(index);
      SCOPED_TRACE(mass);
      widths.at(0).at(index) = model.H_width(mass);
      EXPECT_DOUBLE_EQ(model.H_width(mass + 0.75), widths.at(0).at(index));
    }
    EXPECT_EQ(checksum(widths), UINT64_C(0x03ae53ee31f02801));
    EXPECT_DOUBLE_EQ(model.H_width(1000.9), 639.7);
  }

  TEST_F(HZZ4LTables, HDECAYOutOfRangeDiagnostic) {
    const WidthProbe model;
    ::testing::internal::CaptureStdout();
    EXPECT_DOUBLE_EQ(model.H_width(99.9), 0.);
    EXPECT_DOUBLE_EQ(model.H_width(1001.), 0.);
    EXPECT_DOUBLE_EQ(model.H_width(0.), 0.);
    EXPECT_DOUBLE_EQ(model.H_width(-100.5), 0.);
    const auto output = ::testing::internal::GetCapturedStdout();
    const std::string message = "[RooBWHighMassGGH::H_width]: mass is out of range\n";
    EXPECT_EQ(output, message + message + message + message);
  }

  // All knots, all interval midpoints, and one point beyond each end.
  constexpr std::array<std::array<double, 2>, 41> splineGGH = {{
      {{199.0, 0.1397295381983673}},  {{200.0, 0.139988}},
      {{225.0, 0.1520421318359375}},  {{250.0, 0.17804079375}},
      {{275.0, 0.23456097949218752}}, {{300.0, 0.35250421875}},
      {{325.0, 0.5484333861328126}},  {{350.0, 0.80198245625}},
      {{375.0, 1.078102904296875}},   {{400.0, 1.351339125}},
      {{425.0, 1.5908943623046876}},  {{450.0, 1.75173965625}},
      {{475.0, 1.8197736357421876}},  {{500.0, 1.84162609375}},
      {{525.0, 1.8619363943359375}},  {{550.0, 1.86637056875}},
      {{575.0, 1.831480505078125}},   {{600.0, 1.7805203}},
      {{625.0, 1.7474705029296875}},  {{650.0, 1.72970253125}},
      {{700.0, 1.65961776875}},       {{750.0, 1.5303573}},
      {{800.0, 1.4093194281249999}},  {{850.0, 1.3124163}},
      {{900.0, 1.2163465625}},        {{950.0, 1.1175009500000002}},
      {{1000.0, 1.0265923753125}},    {{1050.0, 0.94667591}},
      {{1100.0, 0.875282065625}},     {{1150.0, 0.81252316}},
      {{1250.0, 0.715752397}},        {{1350.0, 0.640316104}},
      // Captured regression output, not a mathematical constant; retain its exact value.
      // NOLINTNEXTLINE(modernize-use-std-numbers)
      {{1450.0, 0.563636928}},        {{1550.0, 0.490842936}},
      {{1650.0, 0.436131748}},        {{1750.0, 0.396262416}},
      {{1850.0, 0.3583534992}},       {{1950.0, 0.3153438904}},
      {{2050.0, 0.26612762467}},      {{2150.0, 0.210836}},
      {{2151.0, 0.210252585121}},
  }};

  // All knots, all interval midpoints, and one point beyond each end.
  constexpr std::array<std::array<double, 2>, 35> splineVBF8TeV = {{
      {{224.0, 0.029034300568089998}}, {{225.0, 0.0306502}},
      {{325.0, 0.18176401}},           {{425.0, 0.32669308}},
      {{475.0, 0.40362806875}},        {{525.0, 0.48716722999999995}},
      {{575.0, 0.5746747125}},         {{625.0, 0.6421136000000001}},
      {{675.0, 0.672255225}},          {{725.0, 0.6965070999999999}},
      {{775.0, 0.75109185225}},        {{825.0, 0.842850518}},
      {{875.0, 0.9665444875}},         {{925.0, 1.0979954}},
      {{975.0, 1.21198784}},           {{1025.0, 1.29812952}},
      {{1075.0, 1.3555118875}},        {{1125.0, 1.4064209}},
      {{1175.0, 1.4645821250000002}},  {{1225.0, 1.4862188}},
      {{1275.0, 1.4340447875}},        {{1325.0, 1.3541557}},
      {{1375.0, 1.3029129775000001}},  {{1425.0, 1.29446032}},
      {{1475.0, 1.3196263875}},        {{1525.0, 1.3181398}},
      {{1575.0, 1.241918775}},         {{1625.0, 1.1427403}},
      {{1675.0, 1.076109125}},         {{1725.0, 1.012586}},
      {{1775.0, 0.9155380875000001}},  {{1825.0, 0.8044821000000001}},
      {{1875.0, 0.7129827625}},        {{1925.0, 0.674614}},
      {{1926.0, 147.08480165539999}},
  }};

  // All knots, all interval midpoints, and one point beyond each end.
  constexpr std::array<std::array<double, 2>, 29> splineVBF7TeV = {{
      {{249.0, 0.08441914315071}},    {{250.0, 0.0863593}},
      {{350.0, 0.28145649}},          {{450.0, 0.46558721999999997}},
      {{500.0, 0.54742756125}},       {{550.0, 0.61918889}},
      {{600.0, 0.6800992124999999}},  {{650.0, 0.7360762}},
      {{700.0, 0.7908971575}},        {{750.0, 0.8330946300000001}},
      {{800.0, 0.8592538750000001}},  {{850.0, 0.913433}},
      {{900.0, 1.03189925}},          {{950.0, 1.17228}},
      {{1000.0, 1.2805265875}},       {{1050.0, 1.3345225}},
      {{1100.0, 1.3296285374999999}}, {{1150.0, 1.2991986}},
      {{1200.0, 1.2764318499999998}}, {{1250.0, 1.2559464}},
      {{1350.0, 1.1843910599999998}}, {{1450.0, 1.07763328}},
      {{1550.0, 0.9499980000000001}}, {{1650.0, 0.81005156}},
      {{1750.0, 0.6687429300000001}}, {{1850.0, 0.55235056}},
      {{1900.0, 0.512910695}},        {{1950.0, 0.490991}},
      {{1951.0, 195.82395293200003}},
  }};

  TEST_F(HZZ4LTables, GGHSplineConsumers) {
    const RooSigPlusInt signalPlusInterference;
    checkSpline(signalPlusInterference, splineGGH);
    for (bool eightTeV : {false, true}) {
      SCOPED_TRACE(eightTeV);
      const RooCPSHighMassGGH ggh(
          /*name=*/"ggh", /*title=*/"ggh", m_x, m_mass, m_coupling, m_br, m_strength, m_scale, eightTeV);
      const RooCPSHighMassGGHNoInterf gghNo(
          /*name=*/"gghNo", /*title=*/"gghNo", m_x, m_mass, m_coupling, m_br, eightTeV);
      checkSpline(ggh, splineGGH);
      checkSpline(gghNo, splineGGH);
    }
  }

  TEST_F(HZZ4LTables, VBFSplineConsumers) {
    for (bool eightTeV : {false, true}) {
      SCOPED_TRACE(eightTeV);
      const RooCPSHighMassVBF vbf(
          /*name=*/"vbf", /*title=*/"vbf", m_x, m_mass, m_coupling, m_br, m_strength, m_scale, eightTeV);
      const RooCPSHighMassVBFNoInterf vbfNo(
          /*name=*/"vbfNo", /*title=*/"vbfNo", m_x, m_mass, m_coupling, m_br, eightTeV);
      if (eightTeV) {
        checkSpline(vbf, splineVBF8TeV);
        checkSpline(vbfNo, splineVBF8TeV);
      } else {
        checkSpline(vbf, splineVBF7TeV);
        checkSpline(vbfNo, splineVBF7TeV);
      }
    }
  }
}  // namespace

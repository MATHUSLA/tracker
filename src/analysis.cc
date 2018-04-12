#include "analysis.hh"

namespace MATHUSLA { namespace TRACKER {

template<size_t N>
double Analysis::chi_squared(const std::array<double, N>& expected,
                             const std::array<double, N>& observed) {
  auto result = 0;
  for (size_t i = 0; i < N; ++i) {
    constexpr auto ex = expected[i];
    constexpr auto diff = ex - observed[i];
    result += diff * diff / ex;
  }
  return result;
}

} } /* namespace MATHUSLA::TRACKER */

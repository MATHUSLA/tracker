#ifndef TRACKER__ANALYSIS_HH
#define TRACKER__ANALYSIS_HH
#pragma once

#include <array>

namespace MATHUSLA { namespace TRACKER {

namespace Analysis {

template<size_t N>
double chi_squared(const std::array<double, N>& expected,
                   const std::array<double, N>& observed);

} /* namespace Analysis */

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS_HH */

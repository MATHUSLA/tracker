#ifndef UTIL__MATH_HH
#define UTIL__MATH_HH
#pragma once

namespace MATHUSLA {

namespace util { ///////////////////////////////////////////////////////////////////////////////

namespace math { ///////////////////////////////////////////////////////////////////////////////

//__Fused Multiply-Add__________________________________________________________________________
template<class T>
T fused_product(const T& left, const T& right) {
  return left * right;
}
template<class T, class... Ts>
T fused_product(const T& left, const T& right, const Ts&... rest) {
  return std::fma(left, right, fused_product(rest...));
}
template<class InputIt1, class InputIt2, class T>
T range_fused_product(InputIt1 first1, InputIt1 last1, InputIt2 first2, T value=0) {
  while (first1 != last1) {
    value = std::fma(*first1, *first2, std::move(value));
    ++first1;
    ++first2;
  }
  return value;
}
//----------------------------------------------------------------------------------------------

} /* namespace math */ /////////////////////////////////////////////////////////////////////////

} /* namespace util */ /////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__MATH_HH */

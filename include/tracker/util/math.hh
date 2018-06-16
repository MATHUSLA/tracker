/*
 * include/tracker/util/math.hh
 *
 * Copyright 2018 Brandon Gomes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef UTIL__MATH_HH
#define UTIL__MATH_HH
#pragma once

#include <tracker/util/algorithm.hh>

#define FP_FAST_FMA
#define FP_FAST_FMAF
#define FP_FAST_FMAL

namespace MATHUSLA {

namespace util { namespace math { //////////////////////////////////////////////////////////////

//__Fused Multiply-Add__________________________________________________________________________
template<class T>
constexpr const T fused_product(const T& left,
                                const T& right) {
  return left * right;
}
template<class T, class... Ts>
constexpr const T fused_product(const T& left,
                                const T& right,
                                const Ts&... rest) {
  return std::fma(left, right, fused_product(rest...));
}
template<class T, class InputIt1, class InputIt2>
constexpr const T range_fused_product(InputIt1 first1,
                                      InputIt1 last1,
                                      InputIt2 first2,
                                      T value) {
  while (first1 != last1) {
    value = std::fma(*first1, *first2, std::move(value));
    ++first1;
    ++first2;
  }
  return value;
}
//----------------------------------------------------------------------------------------------

//__Sum of Squares using Fused Multiply-Add_____________________________________________________
template<class T, class... Ts>
constexpr const T square(const T& value) {
  return value * value;
}
template<class T, class... Ts>
constexpr const T square(const T& value,
                         const Ts&... rest) {
  return std::fma(value, value, square(rest...));
}
//----------------------------------------------------------------------------------------------

//__Hypotenuse using Fused Multiply-Add_________________________________________________________
template<class T, class... Ts>
constexpr const T hypot(const T& value,
                        const Ts&... rest) {
  return std::sqrt(square(value, rest...));
}
//----------------------------------------------------------------------------------------------

//__Constant Expression Absolute Value__________________________________________________________
template<class T>
constexpr const T abs(const T& value) {
  return std::max(value, -value);
}
//----------------------------------------------------------------------------------------------

//__Within Interval_____________________________________________________________________________
template<class T, class Compare>
constexpr bool within(const T& first,
                      const T& second,
                      const T& interval,
                      Compare comp) {
  return comp(util::math::abs(first - second), interval);
}
template<class T>
constexpr bool within(const T& first,
                      const T& second,
                      const T& interval) {
  return within(first, second, interval, std::less_equal<>());
}
template<class T, class Compare>
constexpr bool within(const T& first,
                      const T& second,
                      const T& low,
                      const T& high,
                      Compare comp) {
  return util::algorithm::between(util::math::abs(first - second), low, high, comp);
}
template<class T>
constexpr bool within(const T& first,
                      const T& second,
                      const T& low,
                      const T& high) {
  return within(first, second, low, high, std::less_equal<>());
}
//----------------------------------------------------------------------------------------------

} } /* namespace util::math */ /////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__MATH_HH */

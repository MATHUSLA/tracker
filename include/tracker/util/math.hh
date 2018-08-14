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
template<class T, class Input>
constexpr const T range_fused_product(Input begin,
                                      Input end) {
  T value{*begin * *(++begin)};
  while (++begin != end)
    value = std::fma(*begin, *++begin, std::move(value));
  return value;
}
template<class Range,
  typename Value = typename Range::value_type>
constexpr const Value range_fused_product(const Range& range) {
  return range_fused_product<Value>(std::cbegin(range), std::cend(range));
}
//----------------------------------------------------------------------------------------------

//__Sum of Squares using Fused Multiply-Add_____________________________________________________
template<class T, class... Ts>
constexpr const T sum_squares(const T& value) {
  return value * value;
}
template<class T, class... Ts>
constexpr const T sum_squares(const T& value,
                              const Ts&... rest) {
  return std::fma(value, value, sum_squares(rest...));
}
template<class T, class Input>
constexpr const T range_sum_squares(Input begin,
                                    Input end) {
  T value{*begin * *begin};
  while (++begin != end)
    value = std::fma(*begin, *begin, std::move(value));
  return value;
}
template<class Range,
  typename Value = typename Range::value_type>
constexpr const Value range_sum_squares(const Range& range) {
  return range_sum_squares<Value>(std::cbegin(range), std::cend(range));
}
//----------------------------------------------------------------------------------------------

//__Hypotenuse using Fused Multiply-Add_________________________________________________________
template<class T, class... Ts>
constexpr const T hypot(const T& value,
                        const Ts&... rest) {
  return std::sqrt(sum_squares(value, rest...));
}
template<class T, class Input>
constexpr const T range_hypot(Input begin,
                              Input end) {
  return std::sqrt(range_sum_squares<T>(begin, end));
}
template<class Range,
  typename Value = typename Range::value_type>
constexpr const Value range_hypot(const Range& range) {
  return range_hypot<Value>(std::cbegin(range), std::cend(range));
}
//----------------------------------------------------------------------------------------------

//__Constant Expression Absolute Value__________________________________________________________
template<class T>
constexpr const T abs(const T& value) {
  return value >= 0 ? value : -value;
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

//__Count Number of Digits in Number____________________________________________________________
template<class T>
std::size_t digits(T t) {
  if (t < 0) t = -t;
  if (t < 10) return 1UL;
  if (t < 100) return 2UL;
  if (t < 1000) return 3UL;
  if (t < 10000) return 4UL;
  if (t < 100000) return 5UL;
  if (t < 1000000) return 6UL;
  return 6UL + digits(t / 1000000);
}
//----------------------------------------------------------------------------------------------

} } /* namespace util::math */ /////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__MATH_HH */

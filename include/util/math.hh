/*
 * include/util/math.hh
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
T range_fused_product(InputIt1 first1, InputIt1 last1, InputIt2 first2, T value) {
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

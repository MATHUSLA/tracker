/*
 * include/util/algorithm.hh
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

#ifndef UTIL__ALGORITHM_HH
#define UTIL__ALGORITHM_HH
#pragma once

namespace MATHUSLA {

namespace util { ///////////////////////////////////////////////////////////////////////////////

namespace algorithm { //////////////////////////////////////////////////////////////////////////

//__Check If Value Within Range_________________________________________________________________
template<class T, class Compare>
constexpr bool between(const T& value,
                       const T& low,
                       const T& high,
                       Compare comp) {
  assert(!comp(high, low));
  return !(comp(value, low) || comp(high, value));
}
template<class T>
constexpr bool between(const T& value,
                       const T& low,
                       const T& high) {
  return between(value, low, high, std::less<>());
}
//----------------------------------------------------------------------------------------------

} /* namespace algorithm */ ////////////////////////////////////////////////////////////////////

} /* namespace util */ /////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__ALGORITHM_HH */

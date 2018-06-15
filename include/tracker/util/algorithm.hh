/*
 * include/tracker/util/algorithm.hh
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

#include <algorithm>
#include <cassert>

namespace MATHUSLA {

namespace util { namespace algorithm { /////////////////////////////////////////////////////////

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

//__Reverse Full Range__________________________________________________________________________
template<class Range>
constexpr void reverse(Range& range) {
  std::reverse(std::begin(range), std::end(range));
}
//----------------------------------------------------------------------------------------------

//__Dependent Copy Range into another Range Through Back Inserter_______________________________
template<class InputRange, class OutputRange, class UnaryFunction>
constexpr UnaryFunction back_insert_copy_if(const InputRange& in,
                                            OutputRange& out,
                                            UnaryFunction f) {
  std::copy_if(std::cbegin(in), std::cend(in), std::back_inserter(out), f);
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Transform Range into another Range Through Back Inserter____________________________________
template<class InputRange, class OutputRange, class UnaryFunction>
constexpr UnaryFunction back_insert_transform(const InputRange& in,
                                              OutputRange& out,
                                              UnaryFunction f) {
  std::transform(std::cbegin(in), std::cend(in), std::back_inserter(out), f);
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Range Subset Inclusion______________________________________________________________________
template<class Range1, class Range2>
constexpr bool includes(const Range1& range1,
                        const Range2& range2) {
  auto first1 = std::cbegin(range1);
  auto last1 = std::cend(range1);
  auto first2 = std::cbegin(range2);
  auto last2 = std::cend(range2);
  for (; first2 != last2; ++first1) {
    if (first1 == last1 || *first2 < *first1) return false;
    if (!(*first1 < *first2)) ++first2;
  }
  return true;
}
template<class Range1, class Range2, class Compare>
constexpr bool includes(const Range1& range1,
                        const Range2& range2,
                        Compare comp) {
  auto first1 = std::cbegin(range1);
  auto last1 = std::cend(range1);
  auto first2 = std::cbegin(range2);
  auto last2 = std::cend(range2);
  for (; first2 != last2; ++first1) {
    if (first1 == last1 || comp(*first2, *first1)) return false;
    if (!comp(*first1, *first2)) ++first2;
  }
  return true;
}
//----------------------------------------------------------------------------------------------

//__General Range Sorting Function______________________________________________________________
template<class Range, class Compare>
Range& sort_range(Range& range,
                  const Compare comp) {
  std::sort(range.begin(), range.end(), comp);
  return range;
}
//----------------------------------------------------------------------------------------------

//__General Range Stable Sorting Function_______________________________________________________
template<class Range, class Compare>
Range& stable_sort_range(Range& range,
                         const Compare comp) {
  std::stable_sort(range.begin(), range.end(), comp);
  return range;
}
//----------------------------------------------------------------------------------------------

//__General Range Copy Sorting Function_________________________________________________________
template<class Range, class Compare>
Range copy_sort_range(const Range& range,
                      const Compare comp) {
  auto copy = range;  //FIXME: unsure how to improve this (maybe: uninitialized_copy)
  std::partial_sort_copy(range.cbegin(), range.cend(), copy.begin(), copy.end(), comp);
  return copy;
}
//----------------------------------------------------------------------------------------------

//__General Range Stable Sorting Function_______________________________________________________
template<class Range, class Compare>
Range stable_copy_sort_range(const Range& range,
                             const Compare comp) {
  auto copy = range;  //FIXME: unsure how to improve this (maybe: uninitialized_copy)
  std::stable_sort(copy.begin(), copy.end(), comp);
  return copy;
}
//----------------------------------------------------------------------------------------------

} } /* namespace util::algorithm */ ////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__ALGORITHM_HH */

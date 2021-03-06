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
constexpr Range& reverse(Range& range) {
  std::reverse(std::begin(range), std::end(range));
  return range;
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

//__Dependent Copy Range into another Range Through Back Inserter_______________________________
template<class InputRange, class OutputRange, class UnaryFunction>
constexpr UnaryFunction back_insert_reverse_copy_if(const InputRange& in,
                                                    OutputRange& out,
                                                    UnaryFunction f) {
  std::copy_if(std::crbegin(in), std::crend(in), std::back_inserter(out), f);
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

//__Transform Range into another Range Through Back Inserter____________________________________
template<class InputRange, class OutputRange, class UnaryFunction>
constexpr UnaryFunction back_insert_reverse_transform(const InputRange& in,
                                                      OutputRange& out,
                                                      UnaryFunction f) {
  std::transform(std::crbegin(in), std::crend(in), std::back_inserter(out), f);
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Copy Elements to Destination Until Predicate is True________________________________________
template<class InputIt, class OutputIt, class UnaryFunction>
constexpr std::pair<InputIt, OutputIt> copy_until(InputIt first,
                                                  InputIt last,
                                                  OutputIt d_first,
                                                  UnaryFunction f) {
  while (first != last) {
    if (f(*first))
      break;
    *d_first++ = *first++;
  }
  return std::make_pair(first, d_first);
}
//----------------------------------------------------------------------------------------------

//__Copy Elements to Destination Until Predicate is True________________________________________
template<class InputRange, class OutputRange, class UnaryFunction>
constexpr UnaryFunction back_insert_copy_until(const InputRange& in,
                                               OutputRange& out,
                                               UnaryFunction f) {
  copy_until(std::cbegin(in), std::cend(out), std::back_inserter(out), f);
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Copy Elements to Destination Until Predicate is True________________________________________
template<class InputRange, class OutputRange, class UnaryFunction>
constexpr UnaryFunction back_insert_reverse_copy_until(const InputRange& in,
                                                       OutputRange& out,
                                                       UnaryFunction f) {
  copy_until(std::crbegin(in), std::crend(out), std::back_inserter(out), f);
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Range Subset Inclusion______________________________________________________________________
// [Implementation from cppreference]
template<class InputIt1, class InputIt2, class Compare=std::less<>>
constexpr bool includes(InputIt1 first1,
                        InputIt1 last1,
                        InputIt2 first2,
                        InputIt2 last2,
                        Compare comp={}) {
  for (; first2 != last2; ++first1) {
    if (first1 == last1 || comp(*first2, *first1))
      return false;
    if (!comp(*first1, *first2))
      ++first2;
  }
  return true;
}
//----------------------------------------------------------------------------------------------

//__Range Subset Inclusion______________________________________________________________________
template<class Range1, class Range2, class Compare=std::less<>>
constexpr bool range_includes(const Range1& range1,
                              const Range2& range2,
                              Compare comp={}) {
  return util::algorithm::includes(std::cbegin(range1), std::cend(range1),
                                   std::cbegin(range2), std::cend(range2), comp);
}
//----------------------------------------------------------------------------------------------

//__General Range Sorting Function______________________________________________________________
template<class Range, class Compare=std::less<>>
Range& sort_range(Range& range,
                  Compare comp={}) {
  std::sort(range.begin(), range.end(), comp);
  return range;
}
//----------------------------------------------------------------------------------------------

//__General Range Stable Sorting Function_______________________________________________________
template<class Range, class Compare=std::less<>>
Range& stable_sort_range(Range& range,
                         Compare comp={}) {
  std::stable_sort(range.begin(), range.end(), comp);
  return range;
}
//----------------------------------------------------------------------------------------------

//__General Range Copy Sorting Function_________________________________________________________
template<class Range, class Compare=std::less<>>
Range copy_sort_range(const Range& range,
                      Compare comp={}) {
  auto copy = range;
  std::partial_sort_copy(range.cbegin(), range.cend(), copy.begin(), copy.end(), comp);
  return copy;
}
//----------------------------------------------------------------------------------------------

//__General Range Stable Sorting Function_______________________________________________________
template<class Range, class Compare=std::less<>>
Range stable_copy_sort_range(const Range& range,
                             Compare comp={}) {
  auto copy = range;
  std::stable_sort(copy.begin(), copy.end(), comp);
  return copy;
}
//----------------------------------------------------------------------------------------------

//__Find Element from Binary Search_____________________________________________________________
// [Implementation from cppreference]
template<class ForwardIt, class T, class Compare=std::less<>>
constexpr ForwardIt binary_find_first(ForwardIt first,
                                      ForwardIt last,
                                      const T& value,
                                      Compare comp={}) {
  // TODO: should be upper bound or lower bound?
  first = std::lower_bound(first, last, value, comp);
  return first != last && !comp(value, *first) ? first : last;
}
//----------------------------------------------------------------------------------------------

//__Find Element from Binary Search in Range____________________________________________________
template<class Range, class T, class Compare=std::less<>>
constexpr typename Range::const_iterator range_binary_find_first(const Range& range,
                                                                 const T& value,
                                                                 Compare comp={}) {
  return binary_find_first(std::cbegin(range), std::cend(range), value, comp);
}
//----------------------------------------------------------------------------------------------

//__Find Element from Binary Search in Range____________________________________________________
template<class Range, class T, class Compare=std::less<>>
constexpr bool binary_search_range(const Range& range,
                                   const T& value,
                                   Compare comp={}) {
  return std::binary_search(std::cbegin(range), std::cend(range), value, comp);
}
//----------------------------------------------------------------------------------------------

} } /* namespace util::algorithm */ ////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__ALGORITHM_HH */

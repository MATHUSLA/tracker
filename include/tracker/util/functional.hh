/*
 * include/tracker/util/functional.hh
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

#ifndef UTIL__FUNCTIONAL_HH
#define UTIL__FUNCTIONAL_HH
#pragma once

#include <functional>

namespace MATHUSLA {

namespace util { namespace functional { ////////////////////////////////////////////////////////

//__Hashable Type Traits________________________________________________________________________
template<class C, typename = void>
struct has_hash_method : std::false_type {};
template<class C>
struct has_hash_method<C, decltype(&C::hash, void())> : std::true_type {};
template<class C>
constexpr bool has_hash_method_v = has_hash_method<C>::value;
//----------------------------------------------------------------------------------------------

//__Hash Function_______________________________________________________________________________
template<class T>
std::enable_if_t<has_hash_method_v<T>, std::size_t> hash(const T& t) {
  return t.hash();
}
template<class T>
std::enable_if_t<!has_hash_method_v<T>, std::size_t> hash(const T& t) {
  return std::hash<T>{}(t);
}
//----------------------------------------------------------------------------------------------

namespace detail { /////////////////////////////////////////////////////////////////////////////

//__General Hash Combiner_______________________________________________________________________
template<class T,
         uint_fast64_t Salt,
         uint_fast64_t LeftShift,
         uint_fast64_t RightShift>
constexpr std::size_t hash_combine(const std::size_t& seed,
                                   const T& t) {
  return seed ^ (hash(t) + Salt + (seed << LeftShift) + (seed >> RightShift));
}
template<class T,
         uint_fast64_t Salt,
         uint_fast64_t LeftShift,
         uint_fast64_t RightShift,
         class ...Args>
constexpr std::size_t hash_combine(const std::size_t& seed,
                                   const T& t,
                                   const Args& ...args) {
  return hash_combine<std::tuple_element_t<0UL, std::tuple<Args...>>, Salt, LeftShift, RightShift>(
           hash_combine<T, Salt, LeftShift, RightShift>(seed, t), args...);
}
//----------------------------------------------------------------------------------------------

//__BOOST Hash Combiner_________________________________________________________________________
template<class T, class ...Args>
constexpr std::size_t boost_hash_combine(const std::size_t& seed,
                                         const T& t,
                                         const Args& ...args) {
  return hash_combine<T, 0x9e3779b9ULL, 6ULL, 2ULL>(seed, t, args...);
}
//----------------------------------------------------------------------------------------------

} /* namespace detail */ ///////////////////////////////////////////////////////////////////////

//__Combine Hash Function_______________________________________________________________________
template<class T, class ...Args>
constexpr std::size_t hash_combine(const std::size_t& seed,
                                   const T& t,
                                   const Args& ...args) {
  return detail::boost_hash_combine<>(seed, t, args...);
}
template<class T, class ...Args>
constexpr std::size_t hash_combine(const T& t,
                                   const Args& ...args) {
  return hash_combine(0UL, t, args...);
}
//----------------------------------------------------------------------------------------------

//__Combine Hash Function and Change Seed_______________________________________________________
template<class T, class ...Args>
constexpr std::size_t hash_combine_modify(std::size_t& seed,
                                          const T& t,
                                          const Args& ...args) {
  return (seed = detail::boost_hash_combine<>(seed, t, args...));
}
//----------------------------------------------------------------------------------------------

//__Combine Hash Function over Range____________________________________________________________
template<class Iter>
constexpr std::size_t hash_combine_range(Iter begin,
                                         Iter end,
                                         const std::size_t& salt={}) {
  std::size_t out{salt};
  while (begin != end)
    hash_combine_modify(out, *begin++);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Combine Hash Function over Range____________________________________________________________
template<class Range>
constexpr std::size_t hash_combine_range(const Range& range,
                                         const std::size_t& salt={}) {
  return hash_combine_range(std::cbegin(range), std::cend(range), salt);
}
//----------------------------------------------------------------------------------------------

} } /* namespace util::functional */ ///////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__FUNCTIONAL_HH */

/*
 * include/tracker/util/type.hh
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

#ifndef UTIL__TYPE_HH
#define UTIL__TYPE_HH
#pragma once

#include <cctype>
#include <iterator>
#include <type_traits>

namespace MATHUSLA {

namespace util { namespace type { //////////////////////////////////////////////////////////////

//__Alternative to std::distance________________________________________________________________
template<class Iter>
constexpr std::size_t distance(const Iter& begin,
                               const Iter& end) {
  return static_cast<std::size_t>(std::distance(begin, end));
}
//----------------------------------------------------------------------------------------------

//__Size Method Typetrait_______________________________________________________________________
template<class C, typename = void>
struct has_size_method : std::false_type {};
template<class C>
struct has_size_method<C, decltype(&C::size, void())> : std::true_type {};
template<class C>
constexpr bool has_size_method_v = has_size_method<C>::value;
//----------------------------------------------------------------------------------------------

//__Get Constexpr Size__________________________________________________________________________
template<class T, std::size_t N>
constexpr std::size_t size(const T (&)[N]) noexcept { return N; }
template<class T>
constexpr std::enable_if_t<has_size_method_v<T>, std::size_t> size(const T& t) {
  return t.size();
}
template<class T>
constexpr std::enable_if_t<!has_size_method_v<T>, std::size_t> size(const T& t) noexcept {
  return distance(std::cbegin(t), std::cend(t));
}
template<class... Ts>
constexpr std::size_t count(const Ts& ...) noexcept { return sizeof...(Ts); }
//----------------------------------------------------------------------------------------------

//__Size Ordering for Sorter____________________________________________________________________
// TODO: use general size function
template<class C = void>
struct size_ordered {
  constexpr bool operator()(const C& a,
                            const C& b) const {
    return a.size() < b.size();
  }
};
template<class C = void>
struct size_less {
  constexpr bool operator()(const C& a,
                            const C& b) const {
    return a.size() < b.size();
  }
};
template<class C = void>
struct size_greater {
  constexpr bool operator()(const C& a,
                            const C& b) const {
    return a.size() > b.size();
  }
};
//----------------------------------------------------------------------------------------------

//__Check Type of Character_____________________________________________________________________
template<class CharT>
constexpr bool isalnum(CharT&& ch) {
  return std::isalnum(static_cast<unsigned char>(std::forward<CharT>(ch)));
}
template<class CharT>
constexpr bool isalpha(CharT&& ch) {
  return std::isalpha(static_cast<unsigned char>(std::forward<CharT>(ch)));
}
template<class CharT>
constexpr bool isblank(CharT&& ch) {
  return std::isblank(static_cast<unsigned char>(std::forward<CharT>(ch)));
}
template<class CharT>
constexpr bool iscntrl(CharT&& ch) {
  return std::iscntrl(static_cast<unsigned char>(std::forward<CharT>(ch)));
}
template<class CharT>
constexpr bool isdigit(CharT&& ch) {
  return std::isdigit(static_cast<unsigned char>(std::forward<CharT>(ch)));
}
template<class CharT>
constexpr bool isgraph(CharT&& ch) {
  return std::isgraph(static_cast<unsigned char>(std::forward<CharT>(ch)));
}
template<class CharT>
constexpr bool islower(CharT&& ch) {
  return std::islower(static_cast<unsigned char>(std::forward<CharT>(ch)));
}
template<class CharT>
constexpr bool isprint(CharT&& ch) {
  return std::isprint(static_cast<unsigned char>(std::forward<CharT>(ch)));
}
template<class CharT>
constexpr bool ispunct(CharT&& ch) {
  return std::ispunct(static_cast<unsigned char>(std::forward<CharT>(ch)));
}
template<class CharT>
constexpr bool isspace(CharT&& ch) {
  return std::isspace(static_cast<unsigned char>(std::forward<CharT>(ch)));
}
//----------------------------------------------------------------------------------------------

//__Primitive Range Class_______________________________________________________________________
template<class Begin, class End=Begin>
struct basic_range { Begin begin; End end; };
//----------------------------------------------------------------------------------------------

} } /* namespace util::type */ /////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__TYPE_HH */

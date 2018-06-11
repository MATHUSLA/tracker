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
#include <type_traits>

namespace MATHUSLA {

namespace util { namespace type { //////////////////////////////////////////////////////////////

//__Parameter Pack Count________________________________________________________________________
template<class... Ts>
struct count { static const std::size_t value = sizeof...(Ts); };
template<class ...Ts>
constexpr std::size_t count_v = count<Ts...>::value;
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

} } /* namespace util::type */ /////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__TYPE_HH */
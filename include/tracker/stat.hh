/*
 * include/tracker/stat.hh
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

#ifndef TRACKER__STAT_HH
#define TRACKER__STAT_HH
#pragma once

#include <type_traits>

#include <tracker/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace stat { ///////////////////////////////////////////////////////////////////////////////

using namespace type;

//__Chi^2 Type Reflection_______________________________________________________________________
template<class C, typename = void>
struct has_chi_squared_member : std::false_type {};
template<class C>
struct has_chi_squared_member<C, decltype(C::chi_squared, void())> : std::true_type {};
template<class C>
constexpr bool has_chi_squared_member_v = has_chi_squared_member<C>::value;

template<class C, typename = void>
struct has_chi_square_member : std::false_type {};
template<class C>
struct has_chi_square_member<C, decltype(C::chi_square, void())> : std::true_type {};
template<class C>
constexpr bool has_chi_square_member_v = has_chi_square_member<C>::value;

template<class C, typename = void>
struct has_chi2_member : std::false_type {};
template<class C>
struct has_chi2_member<C, decltype(C::chi2, void())> : std::true_type {};
template<class C>
constexpr bool has_chi2_member_v = has_chi2_member<C>::value;

template<class C,
  bool = has_chi_squared_member_v<C>
      || has_chi_square_member_v<C>
      || has_chi2_member_v<C>>
struct is_chi_squared_type : std::true_type {};
template<class C>
struct is_chi_squared_type<C, false> : std::false_type {};
template<class C>
constexpr bool is_chi_squared_type_v = is_chi_squared_type<C>::value;
//----------------------------------------------------------------------------------------------

//__Degrees of Freedom Type Reflection__________________________________________________________
template<class C, typename = void>
struct has_degree_of_freedom_member : std::false_type {};
template<class C>
struct has_degree_of_freedom_member<C, decltype(C::degree_of_freedom, void())> : std::true_type {};
template<class C>
constexpr bool has_degree_of_freedom_member_v = has_degree_of_freedom_member<C>::value;

template<class C, typename = void>
struct has_degrees_of_freedom_member : std::false_type {};
template<class C>
struct has_degrees_of_freedom_member<C, decltype(C::degrees_of_freedom, void())> : std::true_type {};
template<class C>
constexpr bool has_degrees_of_freedom_member_v = has_degrees_of_freedom_member<C>::value;

template<class C, typename = void>
struct has_dof_member : std::false_type {};
template<class C>
struct has_dof_member<C, decltype(C::dof, void())> : std::true_type {};
template<class C>
constexpr bool has_dof_member_v = has_dof_member<C>::value;

template<class C,
  bool = has_degree_of_freedom_member_v<C>
      || has_degrees_of_freedom_member_v<C>
      || has_dof_member_v<C>>
struct is_degree_of_freedom_type : std::true_type {};
template<class C>
struct is_degree_of_freedom_type<C, false> : std::false_type {};
template<class C>
constexpr bool is_degree_of_freedom_type_v = is_degree_of_freedom_type<C>::value;
//----------------------------------------------------------------------------------------------

//__Perform Chi^2 Cut on Data___________________________________________________________________
template<class Range,
  typename = std::enable_if_t<is_chi_squared_type_v<typename Range::value_type>>>
void chi_squared_cut(const Range& range,
                     const real min,
                     const real max,
                     Range& out) {
  // TODO: implement
}
//----------------------------------------------------------------------------------------------

} /* namespace stat */ /////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__STAT_HH */

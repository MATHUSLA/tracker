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

#include <tracker/type.hh>

#include <tracker/util/math.hh>
#include <tracker/util/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace stat { ///////////////////////////////////////////////////////////////////////////////

using namespace type;

//__Chi^2 Type Traits___________________________________________________________________________
template<class C, typename = void>
struct has_chi_squared_member : std::false_type {};
template<class C>
struct has_chi_squared_member<C, decltype(C::chi_squared, void())> : std::true_type {};
template<class C>
constexpr bool has_chi_squared_member_v = has_chi_squared_member<C>::value;

template<class C, typename = void>
struct has_chi_squared_method : std::false_type {};
template<class C>
struct has_chi_squared_method<C, decltype(&C::chi_squared, void())> : std::true_type {};
template<class C>
constexpr bool has_chi_squared_method_v = has_chi_squared_method<C>::value;
//----------------------------------------------------------------------------------------------

//__Degree of Freedom Type Traits_______________________________________________________________
template<class C, typename = void>
struct has_degrees_of_freedom_member : std::false_type {};
template<class C>
struct has_degrees_of_freedom_member<C, decltype(C::degrees_of_freedom, void())> : std::true_type {};
template<class C>
constexpr bool has_degrees_of_freedom_member_v = has_degrees_of_freedom_member<C>::value;

template<class C, typename = void>
struct has_degrees_of_freedom_method : std::false_type {};
template<class C>
struct has_degrees_of_freedom_method<C, decltype(&C::degrees_of_freedom, void())> : std::true_type {};
template<class C>
constexpr bool has_degrees_of_freedom_method_v = has_degrees_of_freedom_method<C>::value;
//----------------------------------------------------------------------------------------------

//__Chi^2 and DOF Type Reflection_______________________________________________________________
template<class C,
  bool = has_chi_squared_member_v<C>
      && has_degrees_of_freedom_member_v<C>>
struct has_chi2_and_dof_members : std::true_type {};
template<class C>
struct has_chi2_and_dof_members<C, false> : std::false_type {};
template<class C>
constexpr bool has_chi2_and_dof_members_v = has_chi2_and_dof_members<C>::value;

template<class C,
  bool = has_chi_squared_method_v<C>
      && has_degrees_of_freedom_method_v<C>>
struct has_chi2_and_dof_methods : std::true_type {};
template<class C>
struct has_chi2_and_dof_methods<C, false> : std::false_type {};
template<class C>
constexpr bool has_chi2_and_dof_methods_v = has_chi2_and_dof_methods<C>::value;

template<class C,
  bool =  has_chi2_and_dof_methods_v<C>
      ||  has_chi2_and_dof_members_v<C>
      || (has_chi_squared_member_v<C> && has_degrees_of_freedom_method_v<C>)
      || (has_chi_squared_method_v<C> && has_degrees_of_freedom_member_v<C>)>
struct is_chi2_dof_type : std::true_type {};
template<class C>
struct is_chi2_dof_type<C, false> : std::false_type {};
template<class C>
constexpr bool is_chi2_dof_type_v = is_chi2_dof_type<C>::value;
//----------------------------------------------------------------------------------------------

//__Free Function Alternatives to Memeber Functions_____________________________________________
template<class T>
real chi_squared(const T& t);
template<class T>
real degrees_of_freedom(const T& t);
//----------------------------------------------------------------------------------------------

//__Calculate P-Value from Chi^2________________________________________________________________
real chi_squared_p_value(const real chi2,
                         const size_t dof);
//----------------------------------------------------------------------------------------------

//__Calculate P-Value from Chi^2________________________________________________________________
template<class T>
std::enable_if_t<!is_chi2_dof_type_v<T>, real>
chi_squared_p_value(const T& t) {
  return chi_squared_p_value(chi_squared(t), degrees_of_freedom(t));
}
template<class T>
std::enable_if_t<has_chi2_and_dof_methods_v<T>, real>
chi_squared_p_value(const T& t) {
  return chi_squared_p_value(t.chi_squared(), t.degrees_of_freedom());
}
template<class T>
std::enable_if_t<has_chi2_and_dof_members_v<T>, real>
chi_squared_p_value(const T& t) {
  return chi_squared_p_value(t.chi_squared, t.degrees_of_freedom);
}
//----------------------------------------------------------------------------------------------

//__Perform Chi^2/DOF Cut on Range______________________________________________________________
template<class Range>
std::enable_if_t<!is_chi2_dof_type_v<typename Range::value_type>, Range>&
chi2_per_dof_cut(const Range& range,
                 const real min,
                 const real max,
                 Range& out) {
  util::algorithm::back_insert_copy_if(range, out, [&](const auto& value) {
    return util::algorithm::between(chi_squared(value) / degrees_of_freedom(value), min, max); });
  return out;
}
template<class Range>
std::enable_if_t<has_chi2_and_dof_methods_v<typename Range::value_type>, Range>&
chi2_per_dof_cut(const Range& range,
                 const real min,
                 const real max,
                 Range& out) {
  util::algorithm::back_insert_copy_if(range, out, [&](const auto& value) {
    return util::algorithm::between(value.chi_squared() / value.degrees_of_freedom(), min, max); });
  return out;
}
template<class Range>
std::enable_if_t<has_chi2_and_dof_members_v<typename Range::value_type>, Range>&
chi2_per_dof_cut(const Range& range,
                 const real min,
                 const real max,
                 Range& out) {
  util::algorithm::back_insert_copy_if(range, out, [&](const auto& value) {
    return util::algorithm::between(value.chi_squared / value.degrees_of_freedom, min, max); });
  return out;
}
//----------------------------------------------------------------------------------------------

namespace error { //////////////////////////////////////////////////////////////////////////////

//__Propagate Error_____________________________________________________________________________
template<std::size_t N>
real propagate(const real_array<N>& gradient,
               const real_array<N*N>& covariance) {
  return type::weighted_norm(gradient, covariance);
}
template<std::size_t N>
real propagate(const real_vector& gradient,
               const real_vector& covariance) {
  return propagate(to_array<N>(gradient), to_array<N>(covariance));
}
//----------------------------------------------------------------------------------------------

//__Propagate Error in Sum of Independent Errors________________________________________________
template<class ...Args>
constexpr real propagate_sum(const real error,
                             const Args... rest) {
  return util::math::hypot(error, rest...);
}
//----------------------------------------------------------------------------------------------

//__Propagate Error in Average of Independent Errors____________________________________________
template<class ...Args>
constexpr real propagate_average(const real error,
                                 const Args... rest) {
  return propagate_independent_sum(error, rest...) / std::sqrt(util::type::count_v<Args...>);
}
//----------------------------------------------------------------------------------------------

} /* namespace error */ ////////////////////////////////////////////////////////////////////////

} /* namespace stat */ /////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__STAT_HH */

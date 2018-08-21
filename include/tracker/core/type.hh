/*
 * include/tracker/core/type.hh
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

#ifndef TRACKER__CORE__TYPE_HH
#define TRACKER__CORE__TYPE_HH
#pragma once

#define FP_FAST_FMA
#define FP_FAST_FMAF
#define FP_FAST_FMAL

#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <ostream>
#include <type_traits>
#include <vector>

#include <tracker/util/algorithm.hh>
#include <tracker/util/functional.hh>
#include <tracker/util/math.hh>
#include <tracker/util/type.hh>

namespace MATHUSLA {

namespace type { ///////////////////////////////////////////////////////////////////////////////

//__Numerical Types_____________________________________________________________________________
using integer = long long;
using real = long double;
using integer_range = util::type::basic_range<integer>;
using real_range = util::type::basic_range<real>;
//----------------------------------------------------------------------------------------------

//__Custom Integer Containers___________________________________________________________________
template<std::size_t N>
using integer_array = std::array<real, N>;
using integer_vector = std::vector<integer>;
//----------------------------------------------------------------------------------------------

//__RN Point Types______________________________________________________________________________
struct r2_point { real     x, y;    };
struct r3_point { real     x, y, z; };
struct r4_point { real  t, x, y, z; };
enum class Coordinate { T, X, Y, Z  };
//----------------------------------------------------------------------------------------------

//__RN Point Type Reflection____________________________________________________________________
template<class C, typename = void>
struct has_t_member : std::false_type {};
template<class C>
struct has_t_member<C, decltype(C::t, void())> : std::true_type {};
template<class C>
constexpr bool has_t_member_v = has_t_member<C>::value;

template<class C, typename = void>
struct has_x_member : std::false_type {};
template<class C>
struct has_x_member<C, decltype(C::x, void())> : std::true_type {};
template<class C>
constexpr bool has_x_member_v = has_x_member<C>::value;

template<class C, typename = void>
struct has_y_member : std::false_type {};
template<class C>
struct has_y_member<C, decltype(C::y, void())> : std::true_type {};
template<class C>
constexpr bool has_y_member_v = has_y_member<C>::value;

template<class C, typename = void>
struct has_z_member : std::false_type {};
template<class C>
struct has_z_member<C, decltype(C::z, void())> : std::true_type {};
template<class C>
constexpr bool has_z_member_v = has_z_member<C>::value;

template<class C,
  bool = !has_t_member_v<C>
      &&  has_x_member_v<C>
      &&  has_y_member_v<C>
      && !has_z_member_v<C>>
struct is_r2_type : std::true_type {};
template<class C>
struct is_r2_type<C, false> : std::false_type {};
template<class C>
constexpr bool is_r2_type_v = is_r2_type<C>::value;

template<class C,
  bool = !has_t_member_v<C>
      &&  has_x_member_v<C>
      &&  has_y_member_v<C>
      &&  has_z_member_v<C>>
struct is_r3_type : std::true_type {};
template<class C>
struct is_r3_type<C, false> : std::false_type {};
template<class C>
constexpr bool is_r3_type_v = is_r3_type<C>::value;

template<class C,
  bool = has_t_member_v<C>
      && has_x_member_v<C>
      && has_y_member_v<C>
      && has_z_member_v<C>>
struct is_r4_type : std::true_type {};
template<class C>
struct is_r4_type<C, false> : std::false_type {};
template<class C>
constexpr bool is_r4_type_v = is_r4_type<C>::value;

template<class C,
  bool = is_r2_type_v<C>
      || is_r3_type_v<C>
      || is_r4_type_v<C>>
struct is_rN_type : std::true_type {};
template<class C>
struct is_rN_type<C, false> : std::false_type {};
template<class C>
constexpr bool is_rN_type_v = is_rN_type<C>::value;

template<class C1, class C2,
  bool = is_r2_type_v<C1>
      && is_r2_type_v<C2>>
struct are_both_r2_type : std::true_type {};
template<class C1, class C2>
struct are_both_r2_type<C1, C2, false> : std::false_type {};
template<class C1, class C2>
constexpr bool are_both_r2_type_v = are_both_r2_type<C1, C2>::value;

template<class C1, class C2,
  bool = is_r3_type_v<C1>
      && is_r3_type_v<C2>>
struct are_both_r3_type : std::true_type {};
template<class C1, class C2>
struct are_both_r3_type<C1, C2, false> : std::false_type {};
template<class C1, class C2>
constexpr bool are_both_r3_type_v = are_both_r3_type<C1, C2>::value;

template<class C1, class C2,
  bool = is_r4_type_v<C1>
      && is_r4_type_v<C2>>
struct are_both_r4_type : std::true_type {};
template<class C1, class C2>
struct are_both_r4_type<C1, C2, false> : std::false_type {};
template<class C1, class C2>
constexpr bool are_both_r4_type_v = are_both_r4_type<C1, C2>::value;

template<class C1, class C2,
  bool = are_both_r2_type_v<C1, C2>
      || are_both_r3_type_v<C1, C2>
      || are_both_r4_type_v<C1, C2>>
struct are_both_same_rN_type : std::true_type {};
template<class C1, class C2>
struct are_both_same_rN_type<C1, C2, false> : std::false_type {};
template<class C1, class C2>
constexpr bool are_both_same_rN_type_v = are_both_same_rN_type<C1, C2>::value;
//----------------------------------------------------------------------------------------------

//__Point-Wise Reduction of Dimension___________________________________________________________
template<class T, typename = std::enable_if_t<is_r3_type_v<T> || is_r4_type_v<T>>>
constexpr r2_point reduce_to_r2(const T& point) {
  return {point.x, point.y};
}
template<class T, typename = std::enable_if_t<is_r4_type_v<T>>>
constexpr r3_point reduce_to_r3(const T& point) {
  return {point.x, point.y, point.z};
}
template<class T, typename = std::enable_if_t<is_r4_type_v<T>>>
constexpr r4_point reduce_to_r4(const T& point) {
  return {point.t, point.x, point.y, point.z};
}
//----------------------------------------------------------------------------------------------

//__Select Coordinate Subset of Point___________________________________________________________
template<class T>
std::enable_if_t<is_r2_type_v<T>, real>&
select(T& point,
       const Coordinate c) {
  switch (c) {
    case Coordinate::T: return 0.0L;
    case Coordinate::X: return point.x;
    case Coordinate::Y: return point.y;
    case Coordinate::Z: return 0.0L;
  }
}
template<class T>
std::enable_if_t<is_r3_type_v<T>, real>&
select(T& point,
       const Coordinate c) {
  switch (c) {
    case Coordinate::X: return point.x;
    case Coordinate::Y: return point.y;
    case Coordinate::Z: return point.z;
  }
}
template<class T>
std::enable_if_t<is_r4_type_v<T>, real>&
select(T& point,
       const Coordinate c) {
  switch (c) {
    case Coordinate::T: return point.t;
    case Coordinate::X: return point.x;
    case Coordinate::Y: return point.y;
    case Coordinate::Z: return point.z;
  }
}
template<Coordinate C, class T>
std::enable_if_t<is_r2_type_v<T>, real>&
select(T& point) {
  switch (C) {
    case Coordinate::X: return point.x;
    case Coordinate::Y: return point.y;
  }
}
template<Coordinate C, class T>
std::enable_if_t<is_r3_type_v<T>, real>&
select(T& point) {
  switch (C) {
    case Coordinate::X: return point.x;
    case Coordinate::Y: return point.y;
    case Coordinate::Z: return point.z;
  }
}
template<Coordinate C, class T>
std::enable_if_t<is_r4_type_v<T>, real>&
select(T& point) {
  switch (C) {
    case Coordinate::T: return point.t;
    case Coordinate::X: return point.x;
    case Coordinate::Y: return point.y;
    case Coordinate::Z: return point.z;
  }
}
template<class T>
constexpr std::enable_if_t<is_r2_type_v<T>, real>
select_r1(const T& point,
          const Coordinate c) {
  switch (c) {
    case Coordinate::T: return 0.0L;
    case Coordinate::X: return point.x;
    case Coordinate::Y: return point.y;
    case Coordinate::Z: return 0.0L;
  }
}
template<class T>
constexpr std::enable_if_t<is_r3_type_v<T>, real>
select_r1(const T& point,
          const Coordinate c) {
  switch (c) {
    case Coordinate::T: return 0.0L;
    case Coordinate::X: return point.x;
    case Coordinate::Y: return point.y;
    case Coordinate::Z: return point.z;
  }
}
template<class T>
constexpr std::enable_if_t<is_r4_type_v<T>, real>
select_r1(const T& point,
          const Coordinate c) {
  switch (c) {
    case Coordinate::T: return point.t;
    case Coordinate::X: return point.x;
    case Coordinate::Y: return point.y;
    case Coordinate::Z: return point.z;
  }
}
template<class T>
constexpr std::enable_if_t<is_r3_type_v<T>, r2_point>
select_r2(const T& point,
          const Coordinate x1,
          const Coordinate x2) {
  r2_point out{};
  switch (x1) {
    case Coordinate::T: break;
    case Coordinate::X: out.x = point.x; break;
    case Coordinate::Y: out.x = point.y; break;
    case Coordinate::Z: out.x = point.z; break;
  }
  switch (x2) {
    case Coordinate::T: break;
    case Coordinate::X: out.y = point.x; break;
    case Coordinate::Y: out.y = point.y; break;
    case Coordinate::Z: out.y = point.z; break;
  }
  return out;
}
template<class T>
constexpr std::enable_if_t<is_r4_type_v<T>, r2_point>
select_r2(const T& point,
          const Coordinate x1,
          const Coordinate x2) {
  r2_point out{};
  switch (x1) {
    case Coordinate::T: out.x = point.t; break;
    case Coordinate::X: out.x = point.x; break;
    case Coordinate::Y: out.x = point.y; break;
    case Coordinate::Z: out.x = point.z; break;
  }
  switch (x2) {
    case Coordinate::T: out.y = point.t; break;
    case Coordinate::X: out.y = point.x; break;
    case Coordinate::Y: out.y = point.y; break;
    case Coordinate::Z: out.y = point.z; break;
  }
  return out;
}
template<class T>
constexpr std::enable_if_t<is_r3_type_v<T>, r3_point>
select_r3(const T& point,
          const Coordinate x1,
          const Coordinate x2,
          const Coordinate x3) {
  r3_point out{};
  switch (x1) {
    case Coordinate::T: break;
    case Coordinate::X: out.x = point.x; break;
    case Coordinate::Y: out.x = point.y; break;
    case Coordinate::Z: out.x = point.z; break;
  }
  switch (x2) {
    case Coordinate::T: break;
    case Coordinate::X: out.y = point.x; break;
    case Coordinate::Y: out.y = point.y; break;
    case Coordinate::Z: out.y = point.z; break;
  }
  switch (x3) {
    case Coordinate::T: break;
    case Coordinate::X: out.z = point.x; break;
    case Coordinate::Y: out.z = point.y; break;
    case Coordinate::Z: out.z = point.z; break;
  }
  return out;
}
template<class T>
constexpr std::enable_if_t<is_r4_type_v<T>, r3_point>
select_r3(const T& point,
          const Coordinate x1,
          const Coordinate x2,
          const Coordinate x3) {
  r3_point out{};
  switch (x1) {
    case Coordinate::T: out.x = point.t; break;
    case Coordinate::X: out.x = point.x; break;
    case Coordinate::Y: out.x = point.y; break;
    case Coordinate::Z: out.x = point.z; break;
  }
  switch (x2) {
    case Coordinate::T: out.y = point.t; break;
    case Coordinate::X: out.y = point.x; break;
    case Coordinate::Y: out.y = point.y; break;
    case Coordinate::Z: out.y = point.z; break;
  }
  switch (x3) {
    case Coordinate::T: out.z = point.t; break;
    case Coordinate::X: out.z = point.x; break;
    case Coordinate::Y: out.z = point.y; break;
    case Coordinate::Z: out.z = point.z; break;
  }
  return out;
}
template<class T>
constexpr std::enable_if_t<is_r4_type_v<T>, r4_point>
select_r4(const T& point,
          const Coordinate x1,
          const Coordinate x2,
          const Coordinate x3,
          const Coordinate x4) {
  r4_point out{};
  switch (x1) {
    case Coordinate::T: out.t = point.t; break;
    case Coordinate::X: out.t = point.x; break;
    case Coordinate::Y: out.t = point.y; break;
    case Coordinate::Z: out.t = point.z; break;
  }
  switch (x2) {
    case Coordinate::T: out.x = point.t; break;
    case Coordinate::X: out.x = point.x; break;
    case Coordinate::Y: out.x = point.y; break;
    case Coordinate::Z: out.x = point.z; break;
  }
  switch (x3) {
    case Coordinate::T: out.y = point.t; break;
    case Coordinate::X: out.y = point.x; break;
    case Coordinate::Y: out.y = point.y; break;
    case Coordinate::Z: out.y = point.z; break;
  }
  switch (x4) {
    case Coordinate::T: out.z = point.t; break;
    case Coordinate::X: out.z = point.x; break;
    case Coordinate::Y: out.z = point.y; break;
    case Coordinate::Z: out.z = point.z; break;
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Stream Convenience Printing_________________________________________________________________
inline std::ostream& operator<<(std::ostream& os,
                                const r2_point& point) {
  return os << '(' << point.x << ", " << point.y << ')';
}
inline std::ostream& operator<<(std::ostream& os,
                                const r3_point& point) {
  return os << '(' << point.x << ", " << point.y << ", " << point.z << ')';
}
inline std::ostream& operator<<(std::ostream& os,
                                const r4_point& point) {
  return os << '(' << point.t << ", " << point.x << ", " << point.y << ", " << point.z << ')';
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Negation_________________________________________________________________
template<class T>
std::enable_if_t<is_r2_type_v<T>, T>
operator-(T point) {
  point.x = -point.x;
  point.y = -point.y;
  return point;
}
template<class T>
std::enable_if_t<is_r3_type_v<T>, T>
operator-(T point) {
  point.x = -point.x;
  point.y = -point.y;
  point.z = -point.z;
  return point;
}
template<class T>
std::enable_if_t<is_r4_type_v<T>, T>
operator-(T point) {
  point.t = -point.t;
  point.x = -point.x;
  point.y = -point.y;
  point.z = -point.z;
  return point;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise In-Place Addition________________________________________________________
template<class T1, class T2>
std::enable_if_t<are_both_r2_type_v<T1, T2>, T1>&
operator+=(T1& left,
           const T2& right) {
  left.x += right.x;
  left.y += right.y;
  return left;
}
template<class T1, class T2>
std::enable_if_t<are_both_r3_type_v<T1, T2>, T1>&
operator+=(T1& left,
           const T2& right) {
  left.x += right.x;
  left.y += right.y;
  left.z += right.z;
  return left;
}
template<class T1, class T2>
std::enable_if_t<are_both_r4_type_v<T1, T2>, T1>&
operator+=(T1& left,
           const T2& right) {
  left.t += right.t;
  left.x += right.x;
  left.y += right.y;
  left.z += right.z;
  return left;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Addition_________________________________________________________________
template<class T1, class T2,
  typename = std::enable_if_t<are_both_same_rN_type_v<T1, T2>>>
T1 operator+(T1 left,
             const T2& right) {
  return left += right;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise In-Place Subtraction_____________________________________________________
template<class T1, class T2>
std::enable_if_t<are_both_r2_type_v<T1, T2>, T1>&
operator-=(T1& left,
           const T2& right) {
  left.x -= right.x;
  left.y -= right.y;
  return left;
}
template<class T1, class T2>
std::enable_if_t<are_both_r3_type_v<T1, T2>, T1>&
operator-=(T1& left,
           const T2& right) {
  left.x -= right.x;
  left.y -= right.y;
  left.z -= right.z;
  return left;
}
template<class T1, class T2>
std::enable_if_t<are_both_r4_type_v<T1, T2>, T1>&
operator-=(T1& left,
           const T2& right) {
  left.t -= right.t;
  left.x -= right.x;
  left.y -= right.y;
  left.z -= right.z;
  return left;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Subtraction______________________________________________________________
template<class T1, class T2,
  typename = std::enable_if_t<are_both_same_rN_type_v<T1, T2>>>
T1 operator-(T1 left,
             const T2& right) {
  return left -= right;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise In-Place Scalar Multiplication___________________________________________
template<class T, class A>
std::enable_if_t<is_r2_type_v<T> && std::is_arithmetic<A>::value, T>&
operator*=(T& left,
           const A right) {
  left.x *= right;
  left.y *= right;
  return left;
}
template<class T, class A>
std::enable_if_t<is_r3_type_v<T> && std::is_arithmetic<A>::value, T>&
operator*=(T& left,
           const A right) {
  left.x *= right;
  left.y *= right;
  left.z *= right;
  return left;
}
template<class T, class A>
std::enable_if_t<is_r4_type_v<T> && std::is_arithmetic<A>::value, T>&
operator*=(T& left,
           const A right) {
  left.t *= right;
  left.x *= right;
  left.y *= right;
  left.z *= right;
  return left;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Scalar Multiplication____________________________________________________
template<class T, class A,
  typename = std::enable_if_t<is_rN_type_v<T> && std::is_arithmetic<A>::value>>
T operator*(T left,
            const A right) {
  return left *= right;
}
template<class T, class A,
  typename = std::enable_if_t<is_rN_type_v<T> && std::is_arithmetic<A>::value>>
T operator*(const A left,
            T right) {
  return right *= left;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise In-Place Scalar Division_________________________________________________
template<class T, class A>
std::enable_if_t<is_r2_type_v<T> && std::is_arithmetic<A>::value, T>&
operator/=(T& left,
           const A right) {
  left.x /= right;
  left.y /= right;
  return left;
}
template<class T, class A>
std::enable_if_t<is_r3_type_v<T> && std::is_arithmetic<A>::value, T>&
operator/=(T& left,
           const A right) {
  left.x /= right;
  left.y /= right;
  left.z /= right;
  return left;
}
template<class T, class A>
std::enable_if_t<is_r4_type_v<T> && std::is_arithmetic<A>::value, T>&
operator/=(T& left,
           const A right) {
  left.t /= right;
  left.x /= right;
  left.y /= right;
  left.z /= right;
  return left;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Scalar Division__________________________________________________________
template<class T, class A,
  typename = std::enable_if_t<is_rN_type_v<T> && std::is_arithmetic<A>::value>>
T operator/(T left,
            const A right) {
  return left /= right;
}
template<class T, class A,
  typename = std::enable_if_t<is_rN_type_v<T> && std::is_arithmetic<A>::value>>
T operator/(const A left,
            T right) {
  return right /= left;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Equality_________________________________________________________________
inline constexpr bool operator==(const r2_point& left,
                                 const r2_point& right) {
  return left.x == right.x && left.y == right.y;
}
inline constexpr bool operator==(const r3_point& left,
                                 const r3_point& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}
inline constexpr bool operator==(const r4_point& left,
                                 const r4_point& right) {
  return left.t == right.t && left.x == right.x && left.y == right.y && left.z == right.z;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Inequality_________________________________________________________________
template<class T,
  typename = std::enable_if_t<is_rN_type_v<T>>>
constexpr bool operator!=(const T& left,
                          const T& right) {
  return !(left == right);
}
//----------------------------------------------------------------------------------------------

//__R2/R3/R4 Inner Product______________________________________________________________________
inline real operator*(const r2_point& left,
                      const r2_point& right) {
  return left.x * right.x + left.y * right.y;
}
inline real operator*(const r3_point& left,
                      const r3_point& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}
inline real operator*(const r4_point& left,
                      const r4_point& right) {
  return left.t * right.t + left.x * right.x + left.y * right.y + left.z * right.z;
}
//----------------------------------------------------------------------------------------------

//__R2/R3/R4 Length_____________________________________________________________________________
inline real norm(const r2_point& point) {
  return std::sqrt(point * point);
}
inline real norm(const r3_point& point) {
  return std::sqrt(point * point);
}
inline real norm(const r4_point& point) {
  return std::sqrt(point * point);
}
//----------------------------------------------------------------------------------------------

//__R2/R3 Cross Product_________________________________________________________________________
inline r3_point cross(const r2_point& left,
                      const r2_point& right) {
  return { 0, 0, left.x * right.y - left.y * right.x };
}
inline r3_point cross(const r3_point& left,
                      const r3_point& right) {
  return { left.y * right.z - left.z * right.y,
           left.z * right.x - left.x * right.z,
           left.x * right.y - left.y * right.x };
}
//----------------------------------------------------------------------------------------------

//__RN Point-Line Distance______________________________________________________________________
template<class T, typename = std::enable_if_t<is_rN_type_v<T>>>
real point_line_distance(const T& point,
                         const T& begin,
                         const T& end) {
  const auto delta = begin - point;
  const auto line = end - begin;
  const auto norm2_line = line * line;
  return !norm2_line ? std::numeric_limits<real>::max()
                     : norm(delta - (delta * line) * line / norm2_line);
}
template<class T, typename = std::enable_if_t<is_rN_type_v<T>>>
real point_line_distance(const T& point,
                         const T& begin,
                         const T& end,
                         const Coordinate x1,
                         const Coordinate x2) {
  return point_line_distance(select_r2(point, x1, x2),
                             select_r2(begin, x1, x2),
                             select_r2(end, x1, x2));
}
template<class T, typename = std::enable_if_t<is_rN_type_v<T>>>
real point_line_distance(const T& point,
                         const T& begin,
                         const T& end,
                         const Coordinate x1,
                         const Coordinate x2,
                         const Coordinate x3) {
  return point_line_distance(select_r3(point, x1, x2, x3),
                             select_r3(begin, x1, x2, x3),
                             select_r3(end, x1, x2, x3));
}
//----------------------------------------------------------------------------------------------

//__R3 Interval Check___________________________________________________________________________
template<class T, class I,
  typename = std::enable_if_t<(is_r3_type_v<T> || is_r4_type_v<T>)
                           && (is_r3_type_v<I> || is_r4_type_v<I>)>>
constexpr bool within_dr(const T& first,
                         const T& second,
                         const I& interval) {
  return util::math::within(first.x, second.x, interval.x)
      && util::math::within(first.y, second.y, interval.y)
      && util::math::within(first.z, second.z, interval.z);
}
//----------------------------------------------------------------------------------------------

//__R4 Interval Check___________________________________________________________________________
template<class T,
  typename = std::enable_if_t<is_r4_type_v<T>>>
constexpr bool within_ds(const T& first,
                         const T& second,
                         const T& interval) {
  return util::math::within(first.t, second.t, interval.t)
      && within_dr(first, second, interval);
}
//----------------------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////////////////

//__Array Point Constructs______________________________________________________________________
template<std::size_t N> using r2_point_array = typename std::array<r2_point, N>;
template<std::size_t N> using r3_point_array = typename std::array<r3_point, N>;
template<std::size_t N> using r4_point_array = typename std::array<r4_point, N>;
template<std::size_t N> using     real_array = typename std::array<real, N>;
template<std::size_t N> struct real_array2 { real_array<N>     xs, ys;     };
template<std::size_t N> struct real_array3 { real_array<N>     xs, ys, zs; };
template<std::size_t N> struct real_array4 { real_array<N> ts, xs, ys, zs; };
//----------------------------------------------------------------------------------------------

//__Turn RN Point to ArrayN_____________________________________________________________________
inline constexpr real_array<2> to_indicies(const r2_point& point) {
  return {point.x, point.y};
}
inline constexpr real_array<3> to_indicies(const r3_point& point) {
  return {point.x, point.y, point.z};
}
inline constexpr real_array<4> to_indicies(const r4_point& point) {
  return {point.t, point.x, point.y, point.z};
}
//----------------------------------------------------------------------------------------------

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N>
constexpr r2_point_array<N> reduce_to_r2(const r3_point_array<N>& arr) {
  r2_point_array<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    const auto& arr_i = arr[i];
    out[i] = {arr_i.x, arr_i.y};
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N>
constexpr r2_point_array<N> reduce_to_r2(const r4_point_array<N>& arr) {
  r2_point_array<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    const auto& arr_i = arr[i];
    out[i] = {arr_i.x, arr_i.y};
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N>
constexpr r3_point_array<N> reduce_to_r3(const r4_point_array<N>& arr) {
  r3_point_array<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    const auto& arr_i = arr[i];
    out[i] = {arr_i.x, arr_i.y, arr_i.z};
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N>
constexpr real_array2<N> reduce_to_r2(const real_array3<N>& arr) {
  return {arr.xs, arr.ys};
}
//----------------------------------------------------------------------------------------------

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N>
constexpr real_array2<N> reduce_to_r2(const real_array4<N>& arr) {
  return {arr.xs, arr.ys};
}
//----------------------------------------------------------------------------------------------

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N>
constexpr real_array3<N> reduce_to_r3(const real_array4<N>& arr) {
  return {arr.xs, arr.ys, arr.zs};
}
//----------------------------------------------------------------------------------------------

//__R2-Point Array Transposition________________________________________________________________
template<std::size_t N>
constexpr real_array2<N> transpose(const r2_point_array<N>& arr) {
  real_array2<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    const auto& arr_i = arr[i];
    out.xs[i] = arr_i.x;
    out.ys[i] = arr_i.y;
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Point Array Transposition________________________________________________________________
template<std::size_t N>
constexpr real_array3<N> transpose(const r3_point_array<N>& arr) {
  real_array3<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    const auto& arr_i = arr[i];
    out.xs[i] = arr_i.x;
    out.ys[i] = arr_i.y;
    out.zs[i] = arr_i.z;
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Point Array Transposition________________________________________________________________
template<std::size_t N>
constexpr real_array4<N> transpose(const r4_point_array<N>& arr) {
  real_array4<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    const auto& arr_i = arr[i];
    out.ts[i] = arr_i.t;
    out.xs[i] = arr_i.x;
    out.ys[i] = arr_i.y;
    out.zs[i] = arr_i.z;
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Array Point Transposition________________________________________________________________
template<std::size_t N>
constexpr r2_point_array<N> transpose(const real_array2<N>& arr) {
  r2_point_array<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    auto& out_i = out[i];
    out_i.x = arr.xs[i];
    out_i.y = arr.ys[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Array Point Transposition________________________________________________________________
template<std::size_t N>
constexpr r3_point_array<N> transpose(const real_array3<N>& arr) {
  r3_point_array<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    auto& out_i = out[i];
    out_i.x = arr.xs[i];
    out_i.y = arr.ys[i];
    out_i.z = arr.zs[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Array Point Transposition________________________________________________________________
template<std::size_t N>
constexpr r4_point_array<N> transpose(const real_array4<N>& arr) {
  r4_point_array<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    auto& out_i = out[i];
    out_i.t = arr.ts[i];
    out_i.x = arr.xs[i];
    out_i.y = arr.ys[i];
    out_i.z = arr.zs[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Real-Array Pointwise Sum____________________________________________________________________
template<std::size_t N>
constexpr real_array<N>& operator+=(real_array<N>& left,
                                    const real right) {
  for (auto& l : left)
    l += right;
  return left;
}
template<std::size_t N>
constexpr real_array<N> operator+(real_array<N> left,
                                   const real right) {
  return left += right;
}
template<std::size_t N>
constexpr real_array<N> operator+(const real left,
                                  real_array<N> right) {
  return right += left;
}
//----------------------------------------------------------------------------------------------

//__Real-Array Pointwise Difference_____________________________________________________________
template<std::size_t N>
constexpr real_array<N>& operator-=(real_array<N>& left,
                                    const real right) {
  for (auto& l : left)
    l -= right;
  return left;
}
template<std::size_t N>
constexpr real_array<N> operator-(real_array<N> left,
                                   const real right) {
  return left -= right;
}
template<std::size_t N>
constexpr real_array<N> operator-(const real left,
                                  real_array<N> right) {
  return right -= left;
}
//----------------------------------------------------------------------------------------------

//__Real-Array Pointwise Product________________________________________________________________
template<std::size_t N>
constexpr real_array<N>& operator*=(real_array<N>& left,
                                    const real right) {
  for (auto& l : left)
    l *= right;
  return left;
}
template<std::size_t N>
constexpr real_array<N> operator*(real_array<N> left,
                                   const real right) {
  return left *= right;
}
template<std::size_t N>
constexpr real_array<N> operator*(const real left,
                                  real_array<N> right) {
  return right *= left;
}
//----------------------------------------------------------------------------------------------

//__Real-Array Pointwise Quotient_______________________________________________________________
template<std::size_t N>
constexpr real_array<N>& operator/=(real_array<N>& left,
                                    const real right) {
  for (auto& l : left)
    l /= right;
  return left;
}
template<std::size_t N>
constexpr real_array<N> operator/(real_array<N> left,
                                   const real right) {
  return left /= right;
}
template<std::size_t N>
constexpr real_array<N> operator/(const real left,
                                  real_array<N> right) {
  return right /= left;
}
//----------------------------------------------------------------------------------------------

//__Real-Array Dot Product______________________________________________________________________
template<std::size_t N>
constexpr real operator*(const real_array<N>& left,
                         const real_array<N>& right) {
  return std::inner_product(left.cbegin(), left.cend(), right.cbegin(), 0.0L);
}
//----------------------------------------------------------------------------------------------

//__Real-Array Left Matrix Product______________________________________________________________
template<std::size_t N>
constexpr real_array<N> operator*(const real_array<N * N>& left,
                                  const real_array<N>& right) {
  real_array<N> out;
  for (std::size_t i = 0; i < N; ++i)
    for (std::size_t j = 0; j < N; ++j)
      out[i] = left[N*i+j] * right[j];
  return out;
}
//----------------------------------------------------------------------------------------------

//__Real-Array Right Matrix Product______________________________________________________________
template<std::size_t N>
constexpr real_array<N> operator*(const real_array<N>& left,
                                  const real_array<N * N>& right) {
  real_array<N> out;
  for (std::size_t i = 0; i < N; ++i)
    for (std::size_t j = 0; j < N; ++j)
      out[i] = right[N*j+i] * left[j];
  return out;
}
//----------------------------------------------------------------------------------------------

//__Real-Array Weighted Product_________________________________________________________________
template<std::size_t N>
constexpr real weighted_product(const real_array<N>& left,
                                const real_array<N * N>& weight,
                                const real_array<N>& right) {
  real out{};
  for (std::size_t i = 0; i < N; ++i)
    for (std::size_t j = 0; j < N; ++j)
      out = std::fma(left[i] * weight[N*i+j], right[j], out);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Real-Array Weighted Norm____________________________________________________________________
template<std::size_t N>
constexpr real weighted_norm(const real_array<N>& vector,
                             const real_array<N * N>& weight) {
  return weighted_product<>(vector, weight, vector);
}
//----------------------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////////////////

//__Vector Point Constructs_____________________________________________________________________
using r2_point_vector = typename std::vector<r2_point>;
using r3_point_vector = typename std::vector<r3_point>;
using r4_point_vector = typename std::vector<r4_point>;
using     real_vector = typename std::vector<real>;
struct real_vector2 { real_vector     xs, ys;     };
struct real_vector3 { real_vector     xs, ys, zs; };
struct real_vector4 { real_vector ts, xs, ys, zs; };
//----------------------------------------------------------------------------------------------

//__Vector-Wise Reduction of Dimension__________________________________________________________
inline r2_point_vector reduce_to_r2(const r3_point_vector& vec) {
  r2_point_vector out;
  out.reserve(vec.size());
  std::for_each(vec.cbegin(), vec.cend(),
    [&](const auto& point) { out.push_back(reduce_to_r2(point)); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Vector-Wise Reduction of Dimension__________________________________________________________
inline r2_point_vector reduce_to_r2(const r4_point_vector& vec) {
  r2_point_vector out;
  out.reserve(vec.size());
  std::for_each(vec.cbegin(), vec.cend(),
    [&](const auto& point) { out.push_back(reduce_to_r2(point)); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Vector-Wise Reduction of Dimension__________________________________________________________
inline r3_point_vector reduce_to_r3(const r4_point_vector& vec) {
  r3_point_vector out;
  out.reserve(vec.size());
  std::for_each(vec.cbegin(), vec.cend(),
    [&](const auto& point) { out.push_back(reduce_to_r3(point)); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Vector-Wise Reduction of Dimension__________________________________________________________
inline real_vector2 reduce_to_r2(const real_vector3& vec) {
  return {vec.xs, vec.ys};
}
//----------------------------------------------------------------------------------------------

//__Vector-Wise Reduction of Dimension__________________________________________________________
inline real_vector2 reduce_to_r2(const real_vector4& vec) {
  return {vec.xs, vec.ys};
}
//----------------------------------------------------------------------------------------------

//__Vector-Wise Reduction of Dimension__________________________________________________________
inline real_vector3 reduce_to_r3(const real_vector4& vec) {
  return {vec.xs, vec.ys, vec.zs};
}
//----------------------------------------------------------------------------------------------

//__R2-Point Vector Transposition_______________________________________________________________
inline real_vector2 transpose(const r2_point_vector& vec) {
  const auto size = vec.size();
  real_vector2 out;
  out.xs.reserve(size);
  out.ys.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    const auto& vec_i = vec[i];
    out.xs.push_back(vec_i.x);
    out.ys.push_back(vec_i.y);
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Point Vector Transposition_______________________________________________________________
inline real_vector3 transpose(const r3_point_vector& vec) {
  const auto size = vec.size();
  real_vector3 out;
  out.xs.reserve(size);
  out.ys.reserve(size);
  out.zs.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    const auto& vec_i = vec[i];
    out.xs.push_back(vec_i.x);
    out.ys.push_back(vec_i.y);
    out.zs.push_back(vec_i.z);
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Point Vector Transposition_______________________________________________________________
inline real_vector4 transpose(const r4_point_vector& vec) {
  const auto size = vec.size();
  real_vector4 out;
  out.ts.reserve(size);
  out.xs.reserve(size);
  out.ys.reserve(size);
  out.zs.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    const auto& vec_i = vec[i];
    out.ts.push_back(vec_i.t);
    out.xs.push_back(vec_i.x);
    out.ys.push_back(vec_i.y);
    out.zs.push_back(vec_i.z);
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Vector Point Transposition_______________________________________________________________
inline r2_point_vector transpose(const real_vector2& vec) {
  const auto size = vec.xs.size();
  r2_point_vector out;
  out.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    auto& out_i = out[i];
    out_i.x = vec.xs[i];
    out_i.y = vec.ys[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Vector Point Transposition_______________________________________________________________
inline r3_point_vector transpose(const real_vector3& vec) {
  const auto size = vec.xs.size();
  r3_point_vector out;
  out.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    auto& out_i = out[i];
    out_i.x = vec.xs[i];
    out_i.y = vec.ys[i];
    out_i.z = vec.zs[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Vector Point Transposition_______________________________________________________________
inline r4_point_vector transpose(const real_vector4& vec) {
  const auto size = vec.xs.size();
  r4_point_vector out;
  out.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    auto& out_i = out[i];
    out_i.t = vec.ts[i];
    out_i.x = vec.xs[i];
    out_i.y = vec.ys[i];
    out_i.z = vec.zs[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Real-Vector Pointwise Sum___________________________________________________________________
inline real_vector& operator+=(real_vector& left,
                               const real right) {
  for (auto& l : left)
    l += right;
  return left;
}
inline real_vector operator+(real_vector left,
                             const real right) {
  return left += right;
}
inline real_vector operator+(const real left,
                             real_vector right) {
  return right += left;
}
//----------------------------------------------------------------------------------------------

//__Real-Vector Pointwise Difference____________________________________________________________
inline real_vector& operator-=(real_vector& left,
                               const real right) {
  for (auto& l : left)
    l -= right;
  return left;
}
inline real_vector operator-(real_vector left,
                             const real right) {
  return left -= right;
}
inline real_vector operator-(const real left,
                             real_vector right) {
  return right -= left;
}
//----------------------------------------------------------------------------------------------

//__Real-Vector Pointwise Product_______________________________________________________________
inline real_vector& operator*=(real_vector& left,
                               const real right) {
  for (auto& l : left)
    l *= right;
  return left;
}
inline real_vector operator*(real_vector left,
                             const real right) {
  return left *= right;
}
inline real_vector operator*(const real left,
                             real_vector right) {
  return right *= left;
}
//----------------------------------------------------------------------------------------------

//__Real-Vector Pointwise Quotient______________________________________________________________
inline real_vector& operator/=(real_vector& left,
                               const real right) {
  for (auto& l : left)
    l /= right;
  return left;
}
inline real_vector operator/(real_vector left,
                             const real right) {
  return left /= right;
}
inline real_vector operator/(const real left,
                             real_vector right) {
  return right /= left;
}
//----------------------------------------------------------------------------------------------

//__Real-Vector Dot Product_____________________________________________________________________
inline real operator*(const real_vector& left,
                      const real_vector& right) {
  return std::inner_product(left.cbegin(), left.cend(), right.cbegin(), 0.0L);
}
//----------------------------------------------------------------------------------------------

//__Real-Vector Weighted Product________________________________________________________________
inline real weighted_product(const real_vector& left,
                             const real_vector& weight,
                             const real_vector& right) {
  const auto left_size = left.size();
  const auto right_size = right.size();
  const auto weight_size = weight.size();

  if (left_size != right_size
      || left_size * left_size != weight_size
      || right_size * right_size != weight_size)
    return 0;

  real out{};
  for (std::size_t i = 0; i < left_size; ++i)
    for (std::size_t j = 0; j < left_size; ++j)
      out += left[i] * weight[left_size*i+j] * right[j];
  return out;
}
//----------------------------------------------------------------------------------------------

//__Real-Array Weighted Norm____________________________________________________________________
inline real weighted_norm(const real_vector& vector,
                          const real_vector& weight) {
  return weighted_product(vector, weight, vector);
}
//----------------------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////////////////

//__Real Array to Vector Transformation_________________________________________________________
template<std::size_t N>
real_vector to_vector(const real_array<N>& arr) {
  real_vector out;
  out.reserve(N);
  std::copy(arr.cbegin(), arr.cend(), std::back_inserter(out));
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Point Array to Vector Transformation_____________________________________________________
template<std::size_t N>
r2_point_vector to_vector(const r2_point_array<N>& arr) {
  r2_point_vector out;
  out.reserve(N);
  std::copy(arr.cbegin(), arr.cend(), std::back_inserter(out));
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Point Array to Vector Transformation_____________________________________________________
template<std::size_t N>
r3_point_vector to_vector(const r3_point_array<N>& arr) {
  r3_point_vector out;
  out.reserve(N);
  std::copy(arr.cbegin(), arr.cend(), std::back_inserter(out));
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Point Array to Vector Transformation_____________________________________________________
template<std::size_t N>
r4_point_vector to_vector(const r4_point_array<N>& arr) {
  r4_point_vector out;
  out.reserve(N);
  std::copy(arr.cbegin(), arr.cend(), std::back_inserter(out));
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Array Point to Vector Point Transformation_______________________________________________
template<std::size_t N>
real_vector2 to_vector(const real_array2<N>& arr) {
  real_vector2 out;
  out.xs.reserve(N);
  out.ys.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    out.xs.push_back(arr.xs[i]);
    out.ys.push_back(arr.ys[i]);
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Array Point to Vector Point Transformation_______________________________________________
template<std::size_t N>
real_vector3 to_vector(const real_array3<N>& arr) {
  real_vector3 out;
  out.xs.reserve(N);
  out.ys.reserve(N);
  out.zs.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    out.xs.push_back(arr.xs[i]);
    out.ys.push_back(arr.ys[i]);
    out.zs.push_back(arr.zs[i]);
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Array Point to Vector Point Transformation_______________________________________________
template<std::size_t N>
real_vector4 to_vector(const real_array4<N>& arr) {
  real_vector4 out;
  out.ts.reserve(N);
  out.xs.reserve(N);
  out.ys.reserve(N);
  out.zs.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    out.ts.push_back(arr.ts[i]);
    out.xs.push_back(arr.xs[i]);
    out.ys.push_back(arr.ys[i]);
    out.zs.push_back(arr.zs[i]);
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Real Vector to Array Transformation_________________________________________________________
template<std::size_t N>
constexpr real_array<N> to_array(const real_vector& vec) {
  real_array<N> out;
  std::copy(vec.cbegin(), vec.cend(), out.begin());
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Point Vector to Array Transformation_____________________________________________________
template<std::size_t N>
constexpr r2_point_array<N> to_array(const r2_point_vector& vec) {
  r2_point_array<N> out;
  std::copy(vec.cbegin(), vec.cend(), out.begin());
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Point Vector to Array Transformation_____________________________________________________
template<std::size_t N>
constexpr r3_point_array<N> to_array(const r3_point_vector& vec) {
  r3_point_array<N> out;
  std::copy(vec.cbegin(), vec.cend(), out.begin());
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Point Vector to Array Transformation_____________________________________________________
template<std::size_t N>
constexpr r4_point_array<N> to_array(const r4_point_vector& vec) {
  r4_point_array<N> out;
  std::copy(vec.cbegin(), vec.cend(), out.begin());
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Vector Point to Array Point Transformation_______________________________________________
template<std::size_t N>
constexpr real_array2<N> to_array(const real_vector2& vec) {
  const auto size = std::min(vec.xs.size(), N);
  real_array2<N> out;
  for (std::size_t i = 0; i < size; ++i) {
    out.xs[i] = vec.xs[i];
    out.ys[i] = vec.ys[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Vector Point to Array Point Transformation_______________________________________________
template<std::size_t N>
constexpr real_array3<N> to_array(const real_vector3& vec) {
  const auto size = std::min(vec.xs.size(), N);
  real_array3<N> out;
  for (std::size_t i = 0; i < size; ++i) {
    out.xs[i] = vec.xs[i];
    out.ys[i] = vec.ys[i];
    out.zs[i] = vec.zs[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Vector Point to Array Point Transformation_______________________________________________
template<std::size_t N>
constexpr real_array4<N> to_array(const real_vector4& vec) {
  const auto size = std::min(vec.xs.size(), N);
  real_array4<N> out;
  for (std::size_t i = 0; i < size; ++i) {
    out.ts[i] = vec.ts[i];
    out.xs[i] = vec.xs[i];
    out.ys[i] = vec.ys[i];
    out.zs[i] = vec.zs[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////////////////

//__Coordinate-Wise Sorting Functors____________________________________________________________
template<class C = void>
struct t_ordered {
  constexpr bool operator()(const C& a, const C& b) const { return a.t < b.t; }
};
template<class C = void>
struct x_ordered {
  constexpr bool operator()(const C& a, const C& b) const { return a.x < b.x; }
};
template<class C = void>
struct y_ordered {
  constexpr bool operator()(const C& a, const C& b) const { return a.y < b.y; }
};
template<class C = void>
struct z_ordered {
  constexpr bool operator()(const C& a, const C& b) const { return a.z < b.z; }
};
template<Coordinate R, class C = void>
struct coordinate_ordered;
template<class C>
struct coordinate_ordered<Coordinate::T, C> {
  constexpr bool operator()(const C& a, const C& b) const { return a.t < b.t; }
};
template<class C>
struct coordinate_ordered<Coordinate::X, C> {
  constexpr bool operator()(const C& a, const C& b) const { return a.x < b.x; }
};
template<class C>
struct coordinate_ordered<Coordinate::Y, C> {
  constexpr bool operator()(const C& a, const C& b) const { return a.y < b.y; }
};
template<class C>
struct coordinate_ordered<Coordinate::Z, C> {
  constexpr bool operator()(const C& a, const C& b) const { return a.z < b.z; }
};
//----------------------------------------------------------------------------------------------

//__General Coordinate-Wise Sorting Functions___________________________________________________
template<class Range>
Range& t_sort(Range& range) {
  return util::algorithm::sort_range(range, t_ordered<typename Range::value_type>{});
}
template<class Range>
Range& x_sort(Range& range) {
  return util::algorithm::sort_range(range, x_ordered<typename Range::value_type>{});
}
template<class Range>
Range& y_sort(Range& range) {
  return util::algorithm::sort_range(range, y_ordered<typename Range::value_type>{});
}
template<class Range>
Range& z_sort(Range& range) {
  return util::algorithm::sort_range(range, z_ordered<typename Range::value_type>{});
}
template<Coordinate C, class Range>
Range& coordinate_sort(Range& range) {
  return util::algorithm::sort_range(range, coordinate_ordered<C, typename Range::value_type>{});
}
template<class Range>
Range& coordinate_sort(const Coordinate coordinate,
                       Range& range) {
  switch (coordinate) {
    case Coordinate::T: return t_sort(range);
    case Coordinate::X: return x_sort(range);
    case Coordinate::Y: return y_sort(range);
    case Coordinate::Z: return z_sort(range);
  }
}
//----------------------------------------------------------------------------------------------

//__General Coordinate-Wise Stable Sorting Functions____________________________________________
template<class Range>
Range& t_stable_sort(Range& range) {
  return util::algorithm::stable_sort_range(range, t_ordered<typename Range::value_type>{});
}
template<class Range>
Range& x_stable_sort(Range& range) {
  return util::algorithm::stable_sort_range(range, x_ordered<typename Range::value_type>{});
}
template<class Range>
Range& y_stable_sort(Range& range) {
  return util::algorithm::stable_sort_range(range, y_ordered<typename Range::value_type>{});
}
template<class Range>
Range& z_stable_sort(Range& range) {
  return util::algorithm::stable_sort_range(range, z_ordered<typename Range::value_type>{});
}
template<Coordinate C, class Range>
Range& coordinate_stable_sort(Range& range) {
  return util::algorithm::stable_sort_range(range, coordinate_ordered<C, typename Range::value_type>{});
}
template<class Range>
Range& coordinate_stable_sort(const Coordinate coordinate,
                              Range& range) {
  switch (coordinate) {
    case Coordinate::T: return t_stable_sort(range);
    case Coordinate::X: return x_stable_sort(range);
    case Coordinate::Y: return y_stable_sort(range);
    case Coordinate::Z: return z_stable_sort(range);
  }
}
//----------------------------------------------------------------------------------------------

//__General Coordinate-Wise Copy Sorting Functions______________________________________________
template<class Range>
Range t_copy_sort(const Range& range) {
  return util::algorithm::copy_sort_range(range, t_ordered<typename Range::value_type>{});
}
template<class Range>
Range x_copy_sort(const Range& range) {
  return util::algorithm::copy_sort_range(range, x_ordered<typename Range::value_type>{});
}
template<class Range>
Range y_copy_sort(const Range& range) {
  return util::algorithm::copy_sort_range(range, y_ordered<typename Range::value_type>{});
}
template<class Range>
Range z_copy_sort(const Range& range) {
  return util::algorithm::copy_sort_range(range, z_ordered<typename Range::value_type>{});
}
template<Coordinate C, class Range>
Range coordinate_copy_sort(const Range& range) {
  return util::algorithm::copy_sort_range(range, coordinate_ordered<C, typename Range::value_type>{});
}
template<class Range>
Range coordinate_copy_sort(const Coordinate coordinate,
                           const Range& range) {
  switch (coordinate) {
    case Coordinate::T: return t_copy_sort(range);
    case Coordinate::X: return x_copy_sort(range);
    case Coordinate::Y: return y_copy_sort(range);
    case Coordinate::Z: return z_copy_sort(range);
  }
}
//----------------------------------------------------------------------------------------------

//__General Coordinate-Wise Stable Sorting Functions____________________________________________
template<class Range>
Range t_stable_copy_sort(const Range& range) {
  return util::algorithm::stable_copy_sort_range(range, t_ordered<typename Range::value_type>{});
}
template<class Range>
Range x_stable_copy_sort(const Range& range) {
  return util::algorithm::stable_copy_sort_range(range, x_ordered<typename Range::value_type>{});
}
template<class Range>
Range y_stable_copy_sort(const Range& range) {
  return util::algorithm::stable_copy_sort_range(range, y_ordered<typename Range::value_type>{});
}
template<class Range>
Range z_stable_copy_sort(const Range& range) {
  return util::algorithm::stable_copy_sort_range(range, z_ordered<typename Range::value_type>{});
}
template<Coordinate C, class Range>
Range coordinate_stable_copy_sort(const Range& range) {
  return util::algorithm::stable_copy_sort_range(range, coordinate_ordered<C, typename Range::value_type>{});
}
template<class Range>
Range coordinate_stable_copy_sort(const Coordinate coordinate,
                                  const Range& range) {
  switch (coordinate) {
    case Coordinate::T: return t_stable_copy_sort(range);
    case Coordinate::X: return x_stable_copy_sort(range);
    case Coordinate::Y: return y_stable_copy_sort(range);
    case Coordinate::Z: return z_stable_copy_sort(range);
  }
}
//----------------------------------------------------------------------------------------------

} /* namespace type */ /////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

namespace std { ////////////////////////////////////////////////////////////////////////////////

// TODO: finish types

//__Hash Function for R2 Point__________________________________________________________________
template<>
struct hash<MATHUSLA::type::r2_point> {
  std::size_t operator()(const MATHUSLA::type::r2_point& in) const {
    return MATHUSLA::util::functional::hash_combine(in.x, in.y);
  }
};
//----------------------------------------------------------------------------------------------

//__Hash Function for R3 Point__________________________________________________________________
template<>
struct hash<MATHUSLA::type::r3_point> {
  std::size_t operator()(const MATHUSLA::type::r3_point& in) const {
    return MATHUSLA::util::functional::hash_combine(in.x, in.y, in.z);
  }
};
//----------------------------------------------------------------------------------------------

//__Hash Function for R4 Point__________________________________________________________________
template<>
struct hash<MATHUSLA::type::r4_point> {
  std::size_t operator()(const MATHUSLA::type::r4_point& in) const {
    return MATHUSLA::util::functional::hash_combine(in.t, in.x, in.y, in.z);
  }
};
//----------------------------------------------------------------------------------------------

//__Hash Function for Real Array________________________________________________________________
template<std::size_t N>
struct hash<MATHUSLA::type::real_array<N>> {
  std::size_t operator()(const MATHUSLA::type::real_array<N>& in) const {
    return MATHUSLA::util::functional::hash_combine_range(in);
  }
};
//----------------------------------------------------------------------------------------------

//__Hash Function for Real Vector_______________________________________________________________
template<>
struct hash<MATHUSLA::type::real_vector> {
  std::size_t operator()(const MATHUSLA::type::real_vector& in) const {
    return MATHUSLA::util::functional::hash_combine_range(in);
  }
};
//----------------------------------------------------------------------------------------------

} /* namespace std */ //////////////////////////////////////////////////////////////////////////

#endif /* TRACKER__CORE__TYPE_HH */

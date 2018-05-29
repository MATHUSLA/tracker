/*
 * include/point.hh
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

#ifndef TRACKER__POINT_HH
#define TRACKER__POINT_HH
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <ostream>
#include <numeric>
#include <vector>

namespace MATHUSLA {

namespace type { ///////////////////////////////////////////////////////////////////////////////

//__Numerical Types_____________________________________________________________________________
using integer = long long;
using real = long double;
//----------------------------------------------------------------------------------------------

//__RN Point Types______________________________________________________________________________
struct r2_point { real    x, y;    };
struct r3_point { real    x, y, z; };
struct r4_point { real t, x, y, z; };
enum class Coordinate { T, X, Y, Z };
//----------------------------------------------------------------------------------------------

//__Point-Wise Reduction of Dimension___________________________________________________________
inline r2_point reduce_to_r2(const r3_point& point) { return { point.x, point.y          }; }
inline r2_point reduce_to_r2(const r4_point& point) { return { point.x, point.y          }; }
inline r3_point reduce_to_r3(const r4_point& point) { return { point.x, point.y, point.z }; }
//----------------------------------------------------------------------------------------------

//__Stream Convenience Printing_________________________________________________________________
inline std::ostream& operator<<(std::ostream& os, const r2_point& point) {
  return os << '(' << point.x << ", " << point.y << ')';
}
inline std::ostream& operator<<(std::ostream& os, const r3_point& point) {
  return os << '(' << point.x << ", " << point.y << ", " << point.z << ')';
}
inline std::ostream& operator<<(std::ostream& os, const r4_point& point) {
  return os << '(' << point.t << ", " << point.x << ", " << point.y << ", " << point.z << ')';
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Negation_________________________________________________________________
inline r2_point operator-(const r2_point& point) { return {           -point.x, -point.y           }; }
inline r3_point operator-(const r3_point& point) { return {           -point.x, -point.y, -point.z }; }
inline r4_point operator-(const r4_point& point) { return { -point.t, -point.x, -point.y, -point.z }; }
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise In-Place Addition________________________________________________________
inline r2_point& operator+=(r2_point& left, const r2_point& right) {
  left.x += right.x;
  left.y += right.y;
  return left;
}
inline r3_point& operator+=(r3_point& left, const r3_point& right) {
  left.x += right.x;
  left.y += right.y;
  left.z += right.z;
  return left;
}
inline r4_point& operator+=(r4_point& left, const r4_point& right) {
  left.t += right.t;
  left.x += right.x;
  left.y += right.y;
  left.z += right.z;
  return left;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Addition_________________________________________________________________
inline r2_point operator+(r2_point left, const r2_point& right) { return left += right; }
inline r3_point operator+(r3_point left, const r3_point& right) { return left += right; }
inline r4_point operator+(r4_point left, const r4_point& right) { return left += right; }
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise In-Place Subtraction_____________________________________________________
inline r2_point& operator-=(r2_point& left, const r2_point& right) {
  left.x -= right.x;
  left.y -= right.y;
  return left;
}
inline r3_point& operator-=(r3_point& left, const r3_point& right) {
  left.x -= right.x;
  left.y -= right.y;
  left.z -= right.z;
  return left;
}
inline r4_point& operator-=(r4_point& left, const r4_point& right) {
  left.t -= right.t;
  left.x -= right.x;
  left.y -= right.y;
  left.z -= right.z;
  return left;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Subtraction______________________________________________________________
inline r2_point operator-(r2_point left, const r2_point& right) { return left -= right; }
inline r3_point operator-(r3_point left, const r3_point& right) { return left -= right; }
inline r4_point operator-(r4_point left, const r4_point& right) { return left -= right; }
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise In-Place Scalar Multiplication___________________________________________
inline r2_point& operator*=(r2_point& left, const real right) {
  left.x *= right;
  left.y *= right;
  return left;
}
inline r3_point& operator*=(r3_point& left, const real right) {
  left.x *= right;
  left.y *= right;
  left.z *= right;
  return left;
}
inline r4_point& operator*=(r4_point& left, const real right) {
  left.t *= right;
  left.x *= right;
  left.y *= right;
  left.z *= right;
  return left;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Scalar Multiplication____________________________________________________
inline r2_point operator*(r2_point left, const real right) { return left *= right; }
inline r2_point operator*(const real left, r2_point right) { return right *= left; }
inline r3_point operator*(r3_point left, const real right) { return left *= right; }
inline r3_point operator*(const real left, r3_point right) { return right *= left; }
inline r4_point operator*(r4_point left, const real right) { return left *= right; }
inline r4_point operator*(const real left, r4_point right) { return right *= left; }
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise In-Place Scalar Division_________________________________________________
inline r2_point& operator/=(r2_point& left, const real right) {
  left.x /= right;
  left.y /= right;
  return left;
}
inline r3_point& operator/=(r3_point& left, const real right) {
  left.x /= right;
  left.y /= right;
  left.z /= right;
  return left;
}
inline r4_point& operator/=(r4_point& left, const real right) {
  left.t /= right;
  left.x /= right;
  left.y /= right;
  left.z /= right;
  return left;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Scalar Division__________________________________________________________
inline r2_point operator/(r2_point left, const real right) { return left /= right; }
inline r2_point operator/(const real left, r2_point right) { return right /= left; }
inline r3_point operator/(r3_point left, const real right) { return left /= right; }
inline r3_point operator/(const real left, r3_point right) { return right /= left; }
inline r4_point operator/(r4_point left, const real right) { return left /= right; }
inline r4_point operator/(const real left, r4_point right) { return right /= left; }
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Equality_________________________________________________________________
inline bool operator==(const r2_point& left, const r2_point& right) {
  return left.x == right.x && left.y == right.y;
}
inline bool operator==(const r3_point& left, const r3_point& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}
inline bool operator==(const r4_point& left, const r4_point& right) {
  return left.t == right.t && left.x == right.x && left.y == right.y && left.z == right.z;
}
//----------------------------------------------------------------------------------------------

//__RN Coordinate-Wise Inequality_________________________________________________________________
inline bool operator!=(const r2_point& left, const r2_point& right) { return !(left == right); }
inline bool operator!=(const r3_point& left, const r3_point& right) { return !(left == right); }
inline bool operator!=(const r4_point& left, const r4_point& right) { return !(left == right); }
//----------------------------------------------------------------------------------------------

//__R2/R3 Inner Product_________________________________________________________________________
inline real operator*(const r2_point& left, const r2_point& right) {
  return left.x * right.x + left.y * right.y;
}
inline real operator*(const r3_point& left, const r3_point& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}
//----------------------------------------------------------------------------------------------

//__R2/R3 Length________________________________________________________________________________
inline real norm(const r2_point& point) {
  return std::sqrt(point * point);
}
inline real norm(const r3_point& point) {
  return std::sqrt(point * point);
}
//----------------------------------------------------------------------------------------------

//__R2/R3 Cross Product_________________________________________________________________________
inline r3_point cross(const r2_point& left, const r2_point& right) {
  return { 0, 0, left.x * right.y - left.y * right.x };
}
inline r3_point cross(const r3_point& left, const r3_point& right) {
  return { left.y * right.z - left.z * right.y,
           left.z * right.x - left.x * right.z,
           left.x * right.y - left.y * right.x };
}
//----------------------------------------------------------------------------------------------

//__R3 Point-Line Distance______________________________________________________________________
inline real point_line_distance(const r3_point& point, const r3_point& begin, const r3_point& end) {
  const auto denominator = norm(begin - end);
  return !denominator ? -1 : norm(cross(point - begin, point - end)) / denominator;
}
inline real point_line_distance(const r4_point& point, const r4_point& begin, const r4_point& end) {
  return point_line_distance(reduce_to_r3(point), reduce_to_r3(begin), reduce_to_r3(end));
}
//----------------------------------------------------------------------------------------------

//__R3 Interval Check___________________________________________________________________________
inline bool within_dr(const r3_point& a, const r3_point& b, const r3_point& dr) {
  return std::abs(a.x - b.x) <= dr.x && std::abs(a.y - b.y) <= dr.y && std::abs(a.z - b.z) <= dr.z;
}
inline bool within_dr(const r4_point& a, const r4_point& b, const r4_point& dr) {
  return std::abs(a.x - b.x) <= dr.x && std::abs(a.y - b.y) <= dr.y && std::abs(a.z - b.z) <= dr.z;
}
//----------------------------------------------------------------------------------------------

//__R4 Interval Check___________________________________________________________________________
inline bool within_ds(const r4_point& a, const r4_point& b, const r4_point& ds) {
  return std::abs(a.t - b.t) <= ds.t && within_dr(a, b, ds);
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

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N> inline r2_point_array<N> reduce_to_r2(const r3_point_array<N>& arr) {
  r2_point_array<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    const auto& arr_i = arr[i];
    out[i] = { arr_i.x, arr_i.y };
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N> inline r2_point_array<N> reduce_to_r2(const r4_point_array<N>& arr) {
  r2_point_array<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    const auto& arr_i = arr[i];
    out[i] = { arr_i.x, arr_i.y };
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N> inline r3_point_array<N> reduce_to_r3(const r4_point_array<N>& arr) {
  r3_point_array<N> out;
  for (std::size_t i = 0; i < N; ++i) {
    const auto& arr_i = arr[i];
    out[i] = { arr_i.x, arr_i.y, arr_i.z };
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N> inline real_array2<N> reduce_to_r2(const real_array3<N>& arr) {
  return { arr.xs, arr.ys };
}
//----------------------------------------------------------------------------------------------

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N> inline real_array2<N> reduce_to_r2(const real_array4<N>& arr) {
  return { arr.xs, arr.ys };
}
//----------------------------------------------------------------------------------------------

//__Array-Wise Reduction of Dimension___________________________________________________________
template<std::size_t N> inline real_array3<N> reduce_to_r3(const real_array4<N>& arr) {
  return { arr.xs, arr.ys, arr.zs };
}
//----------------------------------------------------------------------------------------------

//__R2-Point Array Transposition________________________________________________________________
template<std::size_t N> inline real_array2<N> transpose(const r2_point_array<N>& arr) {
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
template<std::size_t N> inline real_array3<N> transpose(const r3_point_array<N>& arr) {
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
template<std::size_t N> inline real_array4<N> transpose(const r4_point_array<N>& arr) {
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
template<std::size_t N> inline r2_point_array<N> transpose(const real_array2<N>& arr) {
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
template<std::size_t N> inline r3_point_array<N> transpose(const real_array3<N>& arr) {
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
template<std::size_t N> inline r4_point_array<N> transpose(const real_array4<N>& arr) {
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

//__Real-Array Dot Product______________________________________________________________________
template<std::size_t N> inline real operator*(const real_array<N>& left, const real_array<N>& right) {
  return std::inner_product(left.cbegin(), left.cend(), right.cbegin(), 0.0L);
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
    [&](const auto& point) { out.push_back({point.x, point.y}); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Vector-Wise Reduction of Dimension__________________________________________________________
inline r2_point_vector reduce_to_r2(const r4_point_vector& vec) {
  r2_point_vector out;
  out.reserve(vec.size());
  std::for_each(vec.cbegin(), vec.cend(),
    [&](const auto& point) { out.push_back({point.x, point.y}); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Vector-Wise Reduction of Dimension__________________________________________________________
inline r3_point_vector reduce_to_r3(const r4_point_vector& vec) {
  r3_point_vector out;
  out.reserve(vec.size());
  std::for_each(vec.cbegin(), vec.cend(),
    [&](const auto& point) { out.push_back({point.x, point.y, point.z}); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Vector-Wise Reduction of Dimension__________________________________________________________
inline real_vector2 reduce_to_r2(const real_vector3& vec) {
  return { vec.xs, vec.ys };
}
//----------------------------------------------------------------------------------------------

//__Vector-Wise Reduction of Dimension__________________________________________________________
inline real_vector2 reduce_to_r2(const real_vector4& vec) {
  return { vec.xs, vec.ys };
}
//----------------------------------------------------------------------------------------------

//__Vector-Wise Reduction of Dimension__________________________________________________________
inline real_vector3 reduce_to_r3(const real_vector4& vec) {
  return { vec.xs, vec.ys, vec.zs };
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

//__Real-Vector Dot Product_____________________________________________________________________
inline real operator*(const real_vector& left, const real_vector& right) {
  return std::inner_product(left.cbegin(), left.cend(), right.cbegin(), 0.0L);
}
//----------------------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////////////////

//__R2-Point Array to Vector Transformation_____________________________________________________
template<std::size_t N> inline r2_point_vector to_vector(const r2_point_array<N>& arr) {
  r2_point_vector out;
  out.reserve(N);
  std::copy(arr.cbegin(), arr.cend(), std::back_inserter(out));
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Point Array to Vector Transformation_____________________________________________________
template<std::size_t N> inline r3_point_vector to_vector(const r3_point_array<N>& arr) {
  r3_point_vector out;
  out.reserve(N);
  std::copy(arr.cbegin(), arr.cend(), std::back_inserter(out));
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Point Array to Vector Transformation_____________________________________________________
template<std::size_t N> inline r4_point_vector to_vector(const r4_point_array<N>& arr) {
  r4_point_vector out;
  out.reserve(N);
  std::copy(arr.cbegin(), arr.cend(), std::back_inserter(out));
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Array Point to Vector Point Transformation_______________________________________________
template<std::size_t N> inline real_vector2 to_vector(const real_array2<N>& arr) {
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
template<std::size_t N> inline real_vector3 to_vector(const real_array3<N>& arr) {
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
template<std::size_t N> inline real_vector4 to_vector(const real_array4<N>& arr) {
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

//__R2-Point Vector to Array Transformation_____________________________________________________
template<std::size_t N> inline r2_point_array<N> to_array(const r2_point_vector& vec) {
  r2_point_array<N> out;
  std::copy(vec.cbegin(), vec.cend(), out.begin());
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Point Vector to Array Transformation_____________________________________________________
template<std::size_t N> inline r3_point_array<N> to_array(const r3_point_vector& vec) {
  r3_point_array<N> out;
  std::copy(vec.cbegin(), vec.cend(), out.begin());
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Point Vector to Array Transformation_____________________________________________________
template<std::size_t N> inline r4_point_array<N> to_array(const r4_point_vector& vec) {
  r4_point_array<N> out;
  std::copy(vec.cbegin(), vec.cend(), out.begin());
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Vector Point to Array Point Transformation_______________________________________________
template<std::size_t N> inline real_array2<N> to_array(const real_vector2& vec) {
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
template<std::size_t N> inline real_array3<N> to_array(const real_vector3& vec) {
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
template<std::size_t N> inline real_array4<N> to_array(const real_vector4& vec) {
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

//__General Range Sorting Function______________________________________________________________
template<class Range, class Compare> inline Range& sort_range(Range& range, Compare comp) {
  std::sort(range.begin(), range.end(), comp);
  return range;
}
//----------------------------------------------------------------------------------------------

//__General Range Stable Sorting Function_______________________________________________________
template<class Range, class Compare> inline Range& stable_sort_range(Range& range, Compare comp) {
  std::stable_sort(range.begin(), range.end(), comp);
  return range;
}
//----------------------------------------------------------------------------------------------

//__General Range Copy Sorting Function_________________________________________________________
template<class Range, class Compare> inline Range copy_sort_range(const Range& range, Compare comp) {
  auto copy = range;  //FIXME: unsure how to improve this (maybe: uninitialized_copy)
  std::partial_sort_copy(range.cbegin(), range.cend(), copy.begin(), copy.end(), comp);
  return copy;
}
//----------------------------------------------------------------------------------------------

//__General Range Stable Sorting Function_______________________________________________________
template<class Range, class Compare> inline Range stable_copy_sort_range(const Range& range, Compare comp) {
  auto copy = range;  //FIXME: unsure how to improve this (maybe: uninitialized_copy)
  std::stable_sort(copy.begin(), copy.end(), comp);
  return copy;
}
//----------------------------------------------------------------------------------------------

//__Coordinate-Wise Sorting Functors____________________________________________________________
template<class T> struct t_ordered { constexpr bool operator()(const T& a, const T& b) const { return a.t < b.t; } };
template<class T> struct x_ordered { constexpr bool operator()(const T& a, const T& b) const { return a.x < b.x; } };
template<class T> struct y_ordered { constexpr bool operator()(const T& a, const T& b) const { return a.y < b.y; } };
template<class T> struct z_ordered { constexpr bool operator()(const T& a, const T& b) const { return a.z < b.z; } };
//----------------------------------------------------------------------------------------------

//__General Coordinate-Wise Sorting Functions___________________________________________________
template<class Range> inline Range& t_sort(Range& range) {
  return sort_range(range, t_ordered<typename Range::value_type>{});
}
template<class Range> inline Range& x_sort(Range& range) {
  return sort_range(range, x_ordered<typename Range::value_type>{});
}
template<class Range> inline Range& y_sort(Range& range) {
  return sort_range(range, y_ordered<typename Range::value_type>{});
}
template<class Range> inline Range& z_sort(Range& range) {
  return sort_range(range, z_ordered<typename Range::value_type>{});
}
template<class Range> inline Range& coordinate_sort(Range& range, const Coordinate& coordinate) {
  switch (coordinate) {
    case Coordinate::T: return t_sort(range);
    case Coordinate::X: return x_sort(range);
    case Coordinate::Y: return y_sort(range);
    case Coordinate::Z: return z_sort(range);
  }
}
//----------------------------------------------------------------------------------------------

//__General Coordinate-Wise Stable Sorting Functions____________________________________________
template<class Range> inline Range& t_stable_sort(Range& range) {
  return stable_sort_range(range, t_ordered<typename Range::value_type>{});
}
template<class Range> inline Range& x_stable_sort(Range& range) {
  return stable_sort_range(range, x_ordered<typename Range::value_type>{});
}
template<class Range> inline Range& y_stable_sort(Range& range) {
  return stable_sort_range(range, y_ordered<typename Range::value_type>{});
}
template<class Range> inline Range& z_stable_sort(Range& range) {
  return stable_sort_range(range, z_ordered<typename Range::value_type>{});
}
template<class Range> inline Range& coordinate_stable_sort(Range& range, const Coordinate& coordinate) {
  switch (coordinate) {
    case Coordinate::T: return t_stable_sort(range);
    case Coordinate::X: return x_stable_sort(range);
    case Coordinate::Y: return y_stable_sort(range);
    case Coordinate::Z: return z_stable_sort(range);
  }
}
//----------------------------------------------------------------------------------------------

//__General Coordinate-Wise Copy Sorting Functions______________________________________________
template<class Range> inline Range t_copy_sort(const Range& range) {
  return copy_sort_range(range, t_ordered<typename Range::value_type>{});
}
template<class Range> inline Range x_copy_sort(const Range& range) {
  return copy_sort_range(range, x_ordered<typename Range::value_type>{});
}
template<class Range> inline Range y_copy_sort(const Range& range) {
  return copy_sort_range(range, y_ordered<typename Range::value_type>{});
}
template<class Range> inline Range z_copy_sort(const Range& range) {
  return copy_sort_range(range, z_ordered<typename Range::value_type>{});
}
template<class Range> inline Range coordinate_copy_sort(const Range& range, const Coordinate& coordinate) {
  switch (coordinate) {
    case Coordinate::T: return t_copy_sort(range);
    case Coordinate::X: return x_copy_sort(range);
    case Coordinate::Y: return y_copy_sort(range);
    case Coordinate::Z: return z_copy_sort(range);
  }
}
//----------------------------------------------------------------------------------------------

//__General Coordinate-Wise Stable Sorting Functions____________________________________________
template<class Range> inline Range t_stable_copy_sort(const Range& range) {
  return stable_copy_sort_range(range, t_ordered<typename Range::value_type>{});
}
template<class Range> inline Range x_stable_copy_sort(const Range& range) {
  return stable_copy_sort_range(range, x_ordered<typename Range::value_type>{});
}
template<class Range> inline Range y_stable_copy_sort(const Range& range) {
  return stable_copy_sort_range(range, y_ordered<typename Range::value_type>{});
}
template<class Range> inline Range z_stable_copy_sort(const Range& range) {
  return stable_copy_sort_range(range, z_ordered<typename Range::value_type>{});
}
template<class Range> inline Range coordinate_stable_copy_sort(const Range& range, const Coordinate& coordinate) {
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

#endif /* TRACKER__POINT_HH */

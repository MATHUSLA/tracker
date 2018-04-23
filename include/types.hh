#ifndef TRACKER__TYPES_HH
#define TRACKER__TYPES_HH
#pragma once

#include <algorithm>
#include <array>
#include <string>
#include <vector>

namespace MATHUSLA {

//__Numerical Types_____________________________________________________________________________
using integer = long;
using real = double;
//----------------------------------------------------------------------------------------------

//__RN Point Types______________________________________________________________________________
struct r2_point { real    x, y;    };
struct r3_point { real    x, y, z; };
struct r4_point { real t, x, y, z; };
enum class Coordinate { T, X, Y, Z };
//----------------------------------------------------------------------------------------------

//__Array Point Constructs______________________________________________________________________
template<std::size_t N> using r2_point_array = typename std::array<r2_point, N>;
template<std::size_t N> using r3_point_array = typename std::array<r3_point, N>;
template<std::size_t N> using r4_point_array = typename std::array<r4_point, N>;
template<std::size_t N> using     real_array = typename std::array<real, N>;
template<std::size_t N> struct real_array2 { real_array<N>     xs, ys;     };
template<std::size_t N> struct real_array3 { real_array<N>     xs, ys, zs; };
template<std::size_t N> struct real_array4 { real_array<N> ts, xs, ys, zs; };
//----------------------------------------------------------------------------------------------

//__R2-Point Array Transposition________________________________________________________________
template<std::size_t N> inline real_array2<N> transpose(const r2_point_array<N>& arr) {
  auto out = real_array2<N>{};
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
  auto out = real_array3<N>{};
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
  auto out = real_array4<N>{};
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
  auto out = r2_point_array<N>{};
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
  auto out = r3_point_array<N>{};
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
  auto out = r4_point_array<N>{};
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

//__Vector Point Constructs_____________________________________________________________________
using r2_point_vector = typename std::vector<r2_point>;
using r3_point_vector = typename std::vector<r3_point>;
using r4_point_vector = typename std::vector<r4_point>;
using     real_vector = typename std::vector<real>;
struct real_vector2 { real_vector     xs, ys;     };
struct real_vector3 { real_vector     xs, ys, zs; };
struct real_vector4 { real_vector ts, xs, ys, zs; };
//----------------------------------------------------------------------------------------------

//__R2-Point Vector Transposition_______________________________________________________________
inline real_vector2 transpose(const r2_point_vector& vec) {
  const auto&& size = vec.size();
  auto out = real_vector2{};
  out.xs.reserve(size);
  out.ys.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    const auto& vec_i = vec[i];
    out.xs[i] = vec_i.x;
    out.ys[i] = vec_i.y;
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Point Vector Transposition_______________________________________________________________
inline real_vector3 transpose(const r3_point_vector& vec) {
  const auto&& size = vec.size();
  auto out = real_vector3{};
  out.xs.reserve(size);
  out.ys.reserve(size);
  out.zs.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    const auto& vec_i = vec[i];
    out.xs[i] = vec_i.x;
    out.ys[i] = vec_i.y;
    out.zs[i] = vec_i.z;
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Point Vector Transposition_______________________________________________________________
inline real_vector4 transpose(const r4_point_vector& vec) {
  const auto&& size = vec.size();
  auto out = real_vector4{};
  out.ts.reserve(size);
  out.xs.reserve(size);
  out.ys.reserve(size);
  out.zs.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    const auto& vec_i = vec[i];
    out.ts[i] = vec_i.t;
    out.xs[i] = vec_i.x;
    out.ys[i] = vec_i.y;
    out.zs[i] = vec_i.z;
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Vector Point Transposition_______________________________________________________________
inline r2_point_vector transpose(const real_vector2& vec) {
  const auto&& size = vec.xs.size();
  auto out = r2_point_vector(size);
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
  const auto&& size = vec.xs.size();
  auto out = r3_point_vector(size);
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
  const auto&& size = vec.xs.size();
  auto out = r4_point_vector(size);
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

//__R2-Point Array to Vector Transformation_____________________________________________________
template<std::size_t N> inline r2_point_vector to_vector(const r2_point_array<N>& arr) {
  auto out = r2_point_vector(N);
  for (std::size_t i = 0; i < N; ++i)
    out[i] = arr[i];
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Point Array to Vector Transformation_____________________________________________________
template<std::size_t N> inline r3_point_vector to_vector(const r3_point_array<N>& arr) {
  auto out = r3_point_vector(N);
  for (std::size_t i = 0; i < N; ++i)
    out[i] = arr[i];
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Point Array to Vector Transformation_____________________________________________________
template<std::size_t N> inline r4_point_vector to_vector(const r4_point_array<N>& arr) {
  auto out = r4_point_vector(N);
  for (std::size_t i = 0; i < N; ++i)
    out[i] = arr[i];
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Array Point to Vector Point Transformation_______________________________________________
template<std::size_t N> inline real_vector2 to_vector(const real_array2<N>& arr) {
  auto out = real_vector2{};
  out.xs.reserve(N);
  out.ys.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    out.xs[i] = arr.xs[i];
    out.ys[i] = arr.ys[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Array Point to Vector Point Transformation_______________________________________________
template<std::size_t N> inline real_vector3 to_vector(const real_array3<N>& arr) {
  auto out = real_vector3{};
  out.xs.reserve(N);
  out.ys.reserve(N);
  out.zs.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    out.xs[i] = arr.xs[i];
    out.ys[i] = arr.ys[i];
    out.zs[i] = arr.zs[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Array Point to Vector Point Transformation_______________________________________________
template<std::size_t N> inline real_vector4 to_vector(const real_array4<N>& arr) {
  auto out = real_vector4{};
  out.ts.reserve(N);
  out.xs.reserve(N);
  out.ys.reserve(N);
  out.zs.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    out.ts[i] = arr.ts[i];
    out.xs[i] = arr.xs[i];
    out.ys[i] = arr.ys[i];
    out.zs[i] = arr.zs[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Point Vector to Array Transformation_____________________________________________________
template<std::size_t N> inline r2_point_array<N> to_array(const r2_point_vector& vec) {
  const auto&& vec_size = vec.size();
  const auto&& size = vec_size < N ? vec_size : N;
  auto out = r2_point_array<N>{};
  for (std::size_t i = 0; i < size; ++i)
    out[i] = vec[i];
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Point Vector to Array Transformation_____________________________________________________
template<std::size_t N> inline r3_point_array<N> to_array(const r3_point_vector& vec) {
  const auto&& vec_size = vec.size();
  const auto&& size = vec_size < N ? vec_size : N;
  auto out = r3_point_array<N>{};
  for (std::size_t i = 0; i < size; ++i)
    out[i] = vec[i];
  return out;
}
//----------------------------------------------------------------------------------------------

//__R4-Point Vector to Array Transformation_____________________________________________________
template<std::size_t N> inline r4_point_array<N> to_array(const r4_point_vector& vec) {
  const auto&& vec_size = vec.size();
  const auto&& size = vec_size < N ? vec_size : N;
  auto out = r4_point_array<N>{};
  for (std::size_t i = 0; i < size; ++i)
    out[i] = vec[i];
  return out;
}
//----------------------------------------------------------------------------------------------

//__R2-Vector Point to Array Point Transformation_______________________________________________
template<std::size_t N> inline real_array2<N> to_array(const real_vector2& vec) {
  const auto&& vec_size = vec.xs.size();
  const auto&& size = vec_size < N ? vec_size : N;
  auto out = real_array2<N>{};
  for (std::size_t i = 0; i < size; ++i) {
    out.xs[i] = vec.xs[i];
    out.ys[i] = vec.ys[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__R3-Vector Point to Array Point Transformation_______________________________________________
template<std::size_t N> inline real_array3<N> to_array(const real_vector3& vec) {
  const auto&& vec_size = vec.xs.size();
  const auto&& size = vec_size < N ? vec_size : N;
  auto out = real_array3<N>{};
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
  const auto&& vec_size = vec.xs.size();
  const auto&& size = vec_size < N ? vec_size : N;
  auto out = real_array4<N>{};
  for (std::size_t i = 0; i < size; ++i) {
    out.ts[i] = vec.ts[i];
    out.xs[i] = vec.xs[i];
    out.ys[i] = vec.ys[i];
    out.zs[i] = vec.zs[i];
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__General Range Sorting Function______________________________________________________________
template<class Range, class Compare> inline void sort_range(Range& range, Compare comp) {
  std::sort(range.begin(), range.end(), comp);
}
//----------------------------------------------------------------------------------------------

//__General Range Stable Sorting Function_______________________________________________________
template<class Range, class Compare> inline void stable_sort_range(Range& range, Compare comp) {
  std::stable_sort(range.begin(), range.end(), comp);
}
//----------------------------------------------------------------------------------------------

//__General Range Copy Sorting Function_________________________________________________________
template<class Range, class Compare> inline Range copy_sort_range(const Range& range, Compare comp) {
  auto copy = range;  //FIXME: unsure how to improve this
  std::partial_sort_copy(range.cbegin(), range.cend(), copy.begin(), copy.end(), comp);
  return copy;
}
//----------------------------------------------------------------------------------------------

//__Coordinate-Wise Sorting Functors____________________________________________________________
template<class T> struct t_sorter { constexpr bool operator()(const T& a, const T& b) const { return a.t < b.t; } };
template<class T> struct x_sorter { constexpr bool operator()(const T& a, const T& b) const { return a.x < b.x; } };
template<class T> struct y_sorter { constexpr bool operator()(const T& a, const T& b) const { return a.y < b.y; } };
template<class T> struct z_sorter { constexpr bool operator()(const T& a, const T& b) const { return a.z < b.z; } };
//----------------------------------------------------------------------------------------------

//__General Coordinate-Wise Sorting Functions___________________________________________________
template<class Range> inline void t_sort(Range& range) {
  sort_range(range, t_sorter<typename Range::value_type>{});
}
template<class Range> inline void x_sort(Range& range) {
  sort_range(range, x_sorter<typename Range::value_type>{});
}
template<class Range> inline void y_sort(Range& range) {
  sort_range(range, y_sorter<typename Range::value_type>{});
}
template<class Range> inline void z_sort(Range& range) {
  sort_range(range, z_sorter<typename Range::value_type>{});
}
template<class Range> inline void coordinate_sort(Range& range, Coordinate coordinate) {
  switch (coordinate) {
    case Coordinate::T: t_sort(range); return;
    case Coordinate::X: x_sort(range); return;
    case Coordinate::Y: y_sort(range); return;
    case Coordinate::Z: z_sort(range); return;
  }
}
//----------------------------------------------------------------------------------------------

//__General Coordinate-Wise Stable Sorting Functions____________________________________________
template<class Range> inline void t_stable_sort(Range& range) {
  stable_sort_range(range, t_sorter<typename Range::value_type>{});
}
template<class Range> inline void x_stable_sort(Range& range) {
  stable_sort_range(range, x_sorter<typename Range::value_type>{});
}
template<class Range> inline void y_stable_sort(Range& range) {
  stable_sort_range(range, y_sorter<typename Range::value_type>{});
}
template<class Range> inline void z_stable_sort(Range& range) {
  stable_sort_range(range, z_sorter<typename Range::value_type>{});
}
template<class Range> inline void coordinate_stable_sort(Range& range, Coordinate coordinate) {
  switch (coordinate) {
    case Coordinate::T: t_stable_sort(range); return;
    case Coordinate::X: x_stable_sort(range); return;
    case Coordinate::Y: y_stable_sort(range); return;
    case Coordinate::Z: z_stable_sort(range); return;
  }
}
//----------------------------------------------------------------------------------------------

//__General Coordinate-Wise Copy Sorting Functions______________________________________________
template<class Range> inline Range t_copy_sort(const Range& range) {
  return copy_sort_range(range, t_sorter<typename Range::value_type>{});
}
template<class Range> inline Range x_copy_sort(const Range& range) {
  return copy_sort_range(range, x_sorter<typename Range::value_type>{});
}
template<class Range> inline Range y_copy_sort(const Range& range) {
  return copy_sort_range(range, y_sorter<typename Range::value_type>{});
}
template<class Range> inline Range z_copy_sort(const Range& range) {
  return copy_sort_range(range, z_sorter<typename Range::value_type>{});
}
template<class Range> inline Range coordinate_copy_sort(const Range& range, Coordinate coordinate) {
  switch (coordinate) {
    case Coordinate::T: return t_copy_sort(range);
    case Coordinate::X: return x_copy_sort(range);
    case Coordinate::Y: return y_copy_sort(range);
    case Coordinate::Z: return z_copy_sort(range);
  }
}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

#endif /* TRACKER__TYPES_HH */

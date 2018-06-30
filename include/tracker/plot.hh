/*
 * include/tracker/plot.hh
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

#ifndef TRACKER__PLOT_HH
#define TRACKER__PLOT_HH
#pragma once

#include <memory>
#include <tuple>
#include <unordered_map>

#include <tracker/core/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace plot { ///////////////////////////////////////////////////////////////////////////////

using namespace type;

//__Plot System_________________________________________________________________________________
void init(bool on=true);
void end();
bool is_on();
//----------------------------------------------------------------------------------------------

//__Simple Color Type___________________________________________________________________________
struct color {
  uint_fast8_t r, g, b;
  static const color WHITE;
  static const color BLACK;
  static const color RED;
  static const color GREEN;
  static const color BLUE;
  static const color CYAN;
  static const color MAGENTA;
  static const color YELLOW;
};
//----------------------------------------------------------------------------------------------

//__Color Equality______________________________________________________________________________
inline bool operator==(const color& left,
                       const color& right) {
  return left.r == right.r && left.g == right.g && left.b == right.b;
}
//----------------------------------------------------------------------------------------------

//__Color Inequality____________________________________________________________________________
inline bool operator!=(const color& left,
                       const color& right) {
  return !(left == right);
}
//----------------------------------------------------------------------------------------------

//__Plotting Histogram Type_____________________________________________________________________
class histogram {
public:
  histogram();
  histogram(const std::string& name,
            const size_t bins,
            const real min,
            const real max);
  histogram(const std::string& name,
            const std::string& title,
            const size_t bins,
            const real min,
            const real max);
  histogram(const std::string& name,
            const std::string& title,
            const std::string& x_title,
            const std::string& y_title,
            const size_t bins,
            const real min,
            const real max);
  histogram(const histogram& other);
  histogram(histogram&& other) noexcept;
  histogram& operator=(const histogram& other);
  histogram& operator=(histogram&& other) noexcept;

  ~histogram();

  const std::string name() const;
  const std::string title() const;
  const std::string x_title() const;
  const std::string y_title() const;

  void name(const std::string& new_name);
  void title(const std::string& new_title);
  void x_title(const std::string& new_x_title);
  void y_title(const std::string& new_y_title);

  bool empty() const;
  size_t size() const;
  real min_x() const;
  real max_x() const;
  real mean() const;

  real bin_value(const size_t index) const;
  real value(const real point) const;

  size_t insert(const real point);
  void increment(const size_t index,
                 const real weight=1);

  void scale(const real weight);

  void draw();
  void clear();

  bool save(const std::string& path) const;

  static histogram load(const std::string& path,
                        const std::string& name="histogram");

  template<class ...Histograms>
  static void draw_all(histogram& h,
                       Histograms& ...hs) {
    draw_all(h);
    draw_all(hs...);
  }
  static void draw_all(histogram& h) { h.draw(); }

  template<class ...Histograms>
  static void clear_all(histogram& h,
                        Histograms& ...hs) {
    clear_all(h);
    clear_all(hs...);
  }
  static void clear_all(histogram& h) { h.clear(); }

  template<class ...Histograms>
  static bool save_all(const std::string& path,
                       const histogram& h,
                       const Histograms& ...hs) {
    return save_all(path, h) && save_all(path, hs...);
  }
  static bool save_all(const std::string& path,
                       const histogram& h) {
    return h.save(path);
  }

private:
  struct impl;
  std::unique_ptr<impl> _impl;
};
//----------------------------------------------------------------------------------------------

//__Histogram Collection Type___________________________________________________________________
class histogram_collection {
public:
  histogram_collection() = default;
  histogram_collection(const std::string& prefix);
  histogram_collection(std::initializer_list<histogram>&& histograms);
  histogram_collection(const std::string& prefix,
                       std::initializer_list<histogram>&& histograms);
  histogram_collection(const histogram_collection& other) = default;
  histogram_collection(histogram_collection&& other) = default;

  histogram_collection& operator=(const histogram_collection& other) = default;
  histogram_collection& operator=(histogram_collection&& other) = default;

  ~histogram_collection() = default;

  histogram& emplace(const histogram& hist);
  histogram& emplace(histogram&& hist);
  histogram& emplace(const std::string& name,
                     const size_t bins,
                     const real min,
                     const real max);
  histogram& emplace(const std::string& name,
                     const std::string& title,
                     const size_t bins,
                     const real min,
                     const real max);
  histogram& emplace(const std::string& name,
                     const std::string& title,
                     const std::string& x_title,
                     const std::string& y_title,
                     const size_t bins,
                     const real min,
                     const real max);

  histogram& load(const std::string& path,
                  const std::string& name="histogram");
  histogram& load_with_prefix(const std::string& path,
                              const std::string& name="histogram");

  histogram& operator[](const std::string& name);

  void draw_all();
  void clear_all();
  bool save_all(const std::string& path) const;

private:
  std::string _prefix;
  std::unordered_map<std::string, histogram> _histograms;
};
//----------------------------------------------------------------------------------------------

//__Plotting Canvas Type________________________________________________________________________
class canvas {
public:
  canvas(const std::string& name="canvas",
         const size_t width=900,
         const size_t height=600);
  canvas(canvas&& other) = default;
  ~canvas();

  canvas& operator=(canvas&& other) = default;

  const std::string name() const;
  size_t width() const;
  size_t height() const;
  bool empty() const;

  void add_point(const real x,
                 const real y,
                 const real z,
                 const real size=1,
                 const color& color=color::BLACK);

  void add_point(const r3_point& point,
                 const real size=1,
                 const color& color=color::BLACK);

  void add_point(const r4_point& point,
                 const real size=1,
                 const color& color=color::BLACK);

  void add_points(const r4_point_vector& points,
                  const real width=1,
                  const color& color=color::BLACK);

  void add_line(const real x1,
                const real y1,
                const real z1,
                const real x2,
                const real y2,
                const real z2,
                const real width=1,
                const color& color=color::BLACK);

  void add_line(const r3_point& first,
                const r3_point& second,
                const real width=1,
                const color& color=color::BLACK);

  void add_line(const r4_point& first,
                const r4_point& second,
                const real width=1,
                const color& color=color::BLACK);

  void add_polyline(const r4_point_vector& points,
                    const real width=1,
                    const color& color=color::BLACK);

  void add_box(const real min_x,
               const real min_y,
               const real min_z,
               const real max_x,
               const real max_y,
               const real max_z,
               const real width=1,
               const color& color=color::BLACK);

  void add_box(const r3_point& min,
               const r3_point& max,
               const real width=1,
               const color& color=color::BLACK);

  void add_box(const r4_point& min,
               const r4_point& max,
               const real width=1,
               const color& color=color::BLACK);

  void add_box(const r3_point& center,
               const real width_x,
               const real width_y,
               const real width_z,
               const real width=1,
               const color& color=color::BLACK);

  void add_box(const r4_point& center,
               const real width_x,
               const real width_y,
               const real width_z,
               const real width=1,
               const color& color=color::BLACK);

  void draw();
  void clear();

  bool save(const std::string& path) const;

  static canvas load(const std::string& path,
                     const std::string& name="canvas");

  template<class ...Canvases>
  static void draw_all(canvas& c,
                       Canvases& ...cs) {
    draw_all(c);
    draw_all(cs...);
  }
  static void draw_all(canvas& c) { c.draw(); }

  template<class ...Canvases>
  static void clear_all(canvas& c,
                        Canvases& ...cs) {
    clear_all(c);
    clear_all(cs...);
  }
  static void clear_all(canvas& c) { c.clear(); }

  template<class ...Canvases>
  static bool save_all(const std::string& path,
                       const canvas& c,
                       const Canvases& ...cs) {
    return save_all(path, c) && save_all(path, cs...);
  }
  static bool save_all(const std::string& path,
                       const canvas& c) {
    return c.save(path);
  }

private:
  struct impl;
  std::unique_ptr<impl> _impl;
};
//----------------------------------------------------------------------------------------------

//__Draw All Drawable Objects___________________________________________________________________
template<class T>
auto draw_all(T& t) -> decltype(t.draw()) {
  t.draw();
}
template<class T>
auto draw_all(T& t) -> decltype (t.draw_all()){
  t.draw_all();
}
template<class T, class ...Ts>
void draw_all(T& t,
              Ts& ...ts) {
  draw_all(t);
  draw_all(ts...);
}
//----------------------------------------------------------------------------------------------

//__Clear All Drawable Objects__________________________________________________________________
template<class T>
auto clear_all(T& t) -> decltype(t.clear()) {
  t.clear();
}
template<class T>
auto clear_all(T& t) -> decltype(t.clear_all()) {
  t.clear_all();
}
template<class T, class ...Ts>
void clear_all(T& t,
               Ts& ...ts) {
  clear_all(t);
  clear_all(ts...);
}
//----------------------------------------------------------------------------------------------

//__Save All Saveable Objects___________________________________________________________________
template<class T>
auto save_all(const std::string& path,
              const T& t) -> decltype(t.save(path)) {
  return t.save(path);
}
template<class T>
auto save_all(const std::string& path,
              const T& t) -> decltype(t.save_all(path)) {
  return t.save_all(path);
}
template<class T, class ...Ts>
bool save_all(const std::string& path,
              const T& t,
              const Ts& ...ts) {
  return save_all(path, t) && save_all(path, ts...);
}
//----------------------------------------------------------------------------------------------

} /* namespace plot */ /////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__PLOT_HH */

#ifndef TRACKER__PLOT_HH
#define TRACKER__PLOT_HH
#pragma once

#include <memory>

#include "point.hh"

namespace MATHUSLA { namespace TRACKER {

namespace plot { ///////////////////////////////////////////////////////////////////////////////

using namespace type;

//__Plot System_________________________________________________________________________________
void init();
void end();
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
inline bool operator==(const color& left, const color& right) {
  return left.r == right.r && left.g == right.g && left.b == right.b;
}
//----------------------------------------------------------------------------------------------

//__Color Inequality____________________________________________________________________________
inline bool operator!=(const color& left, const color& right) {
  return !(left == right);
}
//----------------------------------------------------------------------------------------------

//__Plotting Canvas Type________________________________________________________________________
class canvas {
public:
  canvas(const std::string& name="canvas",
         const integer width=800,
         const integer height=500);
  canvas(canvas&& other) = default;
  ~canvas();

  canvas& operator=(canvas&& other) = default;

  const std::string name() const;
  integer width() const;
  integer height() const;

  void draw();
  void clear();

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

private:
  struct canvas_impl;
  std::unique_ptr<canvas_impl> _impl;
};
//----------------------------------------------------------------------------------------------

namespace root { ///////////////////////////////////////////////////////////////////////////////

//__Export Plot Canvas to ROOT File_____________________________________________________________
void to_file(const canvas& canvas,
             const std::string& path);
//----------------------------------------------------------------------------------------------

//__Import Plot Canvas from ROOT File___________________________________________________________
canvas from_file(const std::string& path);
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

} /* namespace plot */ /////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__PLOT_HH */

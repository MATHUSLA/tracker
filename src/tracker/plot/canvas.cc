/*
 * src/tracker/plot/canvas.cc
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

#include <tracker/plot.hh>

#include <iomanip>
#include <sstream>
#include <unordered_map>

#include <ROOT/TCanvas.h>
#include <ROOT/TPolyLine3D.h>
#include <ROOT/TPolyMarker3D.h>
#include <ROOT/TView3D.h>
#include <ROOT/TAxis3D.h>
#include <ROOT/TColor.h>
#include <ROOT/TFile.h>

#include <tracker/core/units.hh>

namespace MATHUSLA { namespace TRACKER {

namespace plot { ///////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Convert RGB Color to TColor_________________________________________________________________
Int_t _to_TColor_id(const color& color) {
  return TColor::GetColor(color.r, color.g, color.b);
}
//----------------------------------------------------------------------------------------------

//__Style Type__________________________________________________________________________________
struct _style {
  color rgb;
  real size;
};
//----------------------------------------------------------------------------------------------

//__Style Hash Functor__________________________________________________________________________
struct _style_hash {
  std::size_t operator()(const _style& style) const {
    static std::ostringstream s;
    s.str("");
    s.clear();
    s << std::hex
      << std::setfill('0')
      << std::uppercase
      << std::setw(6)
      << ((style.rgb.r << 16) | (style.rgb.g << 8) | style.rgb.b)
      << '_'
      << std::dec
      << std::setw(17)
      << std::setprecision(17)
      << style.size;
    return std::hash<std::string>{}(s.str());
  }
};
//----------------------------------------------------------------------------------------------

//__Style Equality Functor______________________________________________________________________
struct _style_equal {
  bool operator()(const _style& left, const _style& right) const {
    return left.rgb == right.rgb && left.size == right.size;
  }
};
//----------------------------------------------------------------------------------------------

//__Style to Points Map_________________________________________________________________________
using _style_point_map = std::unordered_multimap<_style, r3_point, _style_hash, _style_equal>;
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Canvas Implementation Definition____________________________________________________________
struct canvas::impl {
  TCanvas* _canvas;
  TView3D* _view;
  std::vector<TPolyLine3D*> _poly_lines;
  _style_point_map _polymarker_map;
  bool _has_updated = false;

  void reset_view() {
    _view = dynamic_cast<TView3D*>(TView::CreateView());
    _view->SetAutoRange(true);
  }

  impl(const std::string& name,
       const std::string& title,
       const size_t width,
       const size_t height)
      : _canvas(new TCanvas(name.c_str(), title.c_str(), width, height)) {
    reset_view();
  }

  impl(const impl& other) = default;
  impl(impl&& other) = default;
  impl& operator=(const impl& other) = default;
  impl& operator=(impl&& other) = default;
  ~impl() = default;
};
//----------------------------------------------------------------------------------------------

//__Canvas Constructor__________________________________________________________________________
canvas::canvas(const std::string& name,
               const size_t width,
               const size_t height)
    : canvas(name, name, width, height) {}
//----------------------------------------------------------------------------------------------

//__Canvas Constructor__________________________________________________________________________
canvas::canvas(const std::string& name,
               const std::string& title,
               const size_t width,
               const size_t height)
    : _impl(std::make_unique<impl>(name, title, width, height)) {}
//----------------------------------------------------------------------------------------------

//__Canvas Destructor___________________________________________________________________________
canvas::~canvas() = default;
//----------------------------------------------------------------------------------------------

//__Get Canvas Name_____________________________________________________________________________
const std::string canvas::name() const {
  return _impl->_canvas->GetName();
}
//----------------------------------------------------------------------------------------------

//__Get Canvas Title____________________________________________________________________________
const std::string canvas::title() const {
  return _impl->_canvas->GetTitle();
}
//----------------------------------------------------------------------------------------------

//__Get Canvas Width____________________________________________________________________________
std::size_t canvas::width() const {
  return _impl->_canvas->GetWindowWidth();
}
//----------------------------------------------------------------------------------------------

//__Get Canvas Height___________________________________________________________________________
std::size_t canvas::height() const {
  return _impl->_canvas->GetWindowHeight();
}
//----------------------------------------------------------------------------------------------

//__Set Canvas Name_____________________________________________________________________________
void canvas::name(const std::string& name) {
  _impl->_canvas->SetName(name.c_str());
}
//----------------------------------------------------------------------------------------------

//__Set Canvas Title____________________________________________________________________________
void canvas::title(const std::string& title) {
  _impl->_canvas->SetTitle(title.c_str());
}
//----------------------------------------------------------------------------------------------

//__Set Canvas Width____________________________________________________________________________
void canvas::width(const std::size_t width) {
  // TODO: check sizing
  _impl->_canvas->SetCanvasSize(width, height());
}
//----------------------------------------------------------------------------------------------

//__Set Canvas Height___________________________________________________________________________
void canvas::height(const std::size_t height) {
  // TODO: check sizing
  _impl->_canvas->SetCanvasSize(width(), height);
}
//----------------------------------------------------------------------------------------------

//__Set Canvas Shape____________________________________________________________________________
void canvas::set_shape(const std::size_t width,
                       const std::size_t height) {
  _impl->_canvas->SetCanvasSize(width, height);
}
//----------------------------------------------------------------------------------------------

//__Check if Canvas Has Objects_________________________________________________________________
bool canvas::empty() const {
  return _impl->_poly_lines.empty() && _impl->_polymarker_map.empty();
}
//----------------------------------------------------------------------------------------------

//__Add Point to Canvas_________________________________________________________________________
void canvas::add_point(const real x,
                       const real y,
                       const real z,
                       const real size,
                       const color& color) {
  _impl->_polymarker_map.insert({{color, size}, r3_point{x, y, z} / units::length});
}
//----------------------------------------------------------------------------------------------

//__Add Point to Canvas_________________________________________________________________________
void canvas::add_point(const r3_point& point,
                       const real size,
                       const color& color) {
  add_point(point.x, point.y, point.z, size, color);
}
//----------------------------------------------------------------------------------------------

//__Add Point to Canvas_________________________________________________________________________
void canvas::add_point(const r4_point& point,
                       const real size,
                       const color& color) {
  add_point(reduce_to_r3(point), size, color);
}
//----------------------------------------------------------------------------------------------

//__Add Points to Canvas________________________________________________________________________
void canvas::add_points(const r4_point_vector& points,
                        const real width,
                        const color& color) {
  for (const auto& point : points)
    add_point(point, width, color);
}
//----------------------------------------------------------------------------------------------

//__Add Line to Canvas__________________________________________________________________________
void canvas::add_line(const real x1,
                      const real y1,
                      const real z1,
                      const real x2,
                      const real y2,
                      const real z2,
                      const real width,
                      const color& color) {
  _impl->_canvas->cd();
  auto& lines = _impl->_poly_lines;
  lines.push_back(new TPolyLine3D);
  auto& line = lines.back();
  line->SetNextPoint(x1 / units::length, y1 / units::length, z1 / units::length);
  line->SetNextPoint(x2 / units::length, y2 / units::length, z2 / units::length);
  line->SetLineWidth(width);
  line->SetLineColor(_to_TColor_id(color));
}
//----------------------------------------------------------------------------------------------

//__Add Line to Canvas__________________________________________________________________________
void canvas::add_line(const r3_point& first,
                      const r3_point& second,
                      const real width,
                      const color& color) {
  add_line(first.x, first.y, first.z, second.x, second.y, second.z, width, color);
}
//----------------------------------------------------------------------------------------------

//__Add Line to Canvas__________________________________________________________________________
void canvas::add_line(const r4_point& first,
                      const r4_point& second,
                      const real width,
                      const color& color) {
  add_line(reduce_to_r3(first), reduce_to_r3(second), width, color);
}
//----------------------------------------------------------------------------------------------

//__Add Polyline to Canvas______________________________________________________________________
void canvas::add_polyline(const r4_point_vector& points,
                          const real width,
                          const color& color) {
  for (const auto& point : points)
    add_point(point, width, color);
}
//----------------------------------------------------------------------------------------------

//__Add Box to Canvas___________________________________________________________________________
void canvas::add_box(const real min_x,
                     const real min_y,
                     const real min_z,
                     const real max_x,
                     const real max_y,
                     const real max_z,
                     const real width,
                     const color& color) {
  add_line(min_x, min_y, min_z, max_x, min_y, min_z, width, color);
  add_line(min_x, max_y, min_z, max_x, max_y, min_z, width, color);
  add_line(min_x, min_y, max_z, max_x, min_y, max_z, width, color);
  add_line(min_x, max_y, max_z, max_x, max_y, max_z, width, color);
  add_line(min_x, min_y, min_z, min_x, max_y, min_z, width, color);
  add_line(max_x, min_y, min_z, max_x, max_y, min_z, width, color);
  add_line(min_x, min_y, max_z, min_x, max_y, max_z, width, color);
  add_line(max_x, min_y, max_z, max_x, max_y, max_z, width, color);
  add_line(min_x, min_y, min_z, min_x, min_y, max_z, width, color);
  add_line(max_x, min_y, min_z, max_x, min_y, max_z, width, color);
  add_line(min_x, max_y, min_z, min_x, max_y, max_z, width, color);
  add_line(max_x, max_y, min_z, max_x, max_y, max_z, width, color);
}
//----------------------------------------------------------------------------------------------

//__Add Box to Canvas___________________________________________________________________________
void canvas::add_box(const r3_point& min,
                     const r3_point& max,
                     const real width,
                     const color& color) {
  add_box(min.x, min.y, min.z, max.x, max.y, max.z, width, color);
}
//----------------------------------------------------------------------------------------------

//__Add Box to Canvas___________________________________________________________________________
void canvas::add_box(const r4_point& min,
                     const r4_point& max,
                     const real width,
                     const color& color) {
  add_box(reduce_to_r3(min), reduce_to_r3(max), width, color);
}
//----------------------------------------------------------------------------------------------

//__Add Box to Canvas___________________________________________________________________________
void canvas::add_box(const r3_point& center,
                     const real width_x,
                     const real width_y,
                     const real width_z,
                     const real width,
                     const color& color) {
  const auto half_widths = r3_point{width_x, width_y, width_z} / 2.0L;
  add_box(center - half_widths, center + half_widths, width, color);
}
//----------------------------------------------------------------------------------------------

//__Add Box to Canvas___________________________________________________________________________
void canvas::add_box(const r4_point& center,
                     const real width_x,
                     const real width_y,
                     const real width_z,
                     const real width,
                     const color& color) {
  add_box(reduce_to_r3(center), width_x, width_y, width_z, width, color);
}
//----------------------------------------------------------------------------------------------

//__Draw Canvas_________________________________________________________________________________
void canvas::draw() {
  _impl->_canvas->cd();

  const auto& marker_map = _impl->_polymarker_map;
  const auto marker_map_size = marker_map.bucket_count();
  for (size_t i = 0; i < marker_map_size; ++i) {
    auto polymarker = new TPolyMarker3D(marker_map.bucket_size(i), 20);
    const auto& begin = marker_map.cbegin(i);
    const auto& end = marker_map.cend(i);
    std::for_each(begin, end, [&](const auto& entry) {
      const auto& point = entry.second;
      polymarker->SetNextPoint(point.x, point.y, point.z);
    });
    if (begin != end) {
      const auto& style = (*begin).first;
      polymarker->SetMarkerSize(style.size);
      polymarker->SetMarkerColor(_to_TColor_id(style.rgb));
      polymarker->Draw();
    }
  }

  for (const auto& poly_line : _impl->_poly_lines)
    poly_line->Draw();

  if (!_impl->_has_updated) {
    _impl->_view->ShowAxis();
    _impl->_canvas->cd();
    auto axis = TAxis3D::GetPadAxis();
    if (axis) {
      axis->SetLabelColor(kBlack);
      axis->SetAxisColor(kBlack);
      axis->SetTitleOffset(2);
      axis->SetXTitle(("X (" + units::length_string + ")").c_str());
      axis->SetYTitle(("Y (" + units::length_string + ")").c_str());
      axis->SetZTitle(("Z (" + units::length_string + ")").c_str());
    }
  }

  _impl->_canvas->Modified();
  _impl->_canvas->Update();
  _impl->_has_updated = true;
}
//----------------------------------------------------------------------------------------------

//__Clear A Canvas______________________________________________________________________________
void canvas::clear() {
  if (_impl->_has_updated) {
    _impl->_canvas->cd();
    _impl->_canvas->Clear();
    _impl->_canvas->Modified();
    _impl->_canvas->Update();
    _impl->reset_view();
    _impl->_polymarker_map.clear();
    _impl->_poly_lines.clear();
    _impl->_has_updated = false;
  }
}
//----------------------------------------------------------------------------------------------

//__Save Canvas to ROOT File____________________________________________________________________
bool canvas::save(const std::string& path) const {
  TFile file(path.c_str(), "UPDATE");
  if (!file.IsZombie()) {
    file.cd();
    file.WriteTObject(_impl->_canvas->Clone());
    file.Close();
    return true;
  }
  return false;
}
//----------------------------------------------------------------------------------------------

//__Construct Canvas from File__________________________________________________________________
canvas canvas::load(const std::string& path,
                    const std::string& name) {
  canvas out;
  TFile file(path.c_str(), "READ");
  if (!file.IsZombie()) {
    TCanvas* test = nullptr;
    file.GetObject(name.c_str(), test);
    if (test) {
      out._impl->_canvas = test;
      out._impl->reset_view();
    }
    file.Close();
  }
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace plot */  ////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */
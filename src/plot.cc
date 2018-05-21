#include "plot.hh"

#include <iomanip>
#include <sstream>
#include <unordered_map>

#include "ROOT/TApplication.h"
#include "ROOT/TCanvas.h"
#include "ROOT/TPolyLine3D.h"
#include "ROOT/TPolyMarker3D.h"
#include "ROOT/TView3D.h"
#include "ROOT/TAxis3D.h"
#include "ROOT/TColor.h"

#include <iostream>

namespace MATHUSLA { namespace TRACKER {

namespace plot { ///////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Plotting Application________________________________________________________________________
TApplication* _app = nullptr;
//----------------------------------------------------------------------------------------------

//__Convert RGB Color to TColor_________________________________________________________________
Int_t _color_to_TColor_free_index(const color& color) {
  const auto&& free_index = TColor::GetFreeColorIndex();
  new TColor(free_index, color.r/255.0L, color.g/255.0L, color.b/255.0L);
  return free_index;
}
//----------------------------------------------------------------------------------------------

//__Style Type__________________________________________________________________________________
struct _style {
  color color;
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
      << ((style.color.r << 16) | (style.color.g << 8) | style.color.b)
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
    return left.color == right.color && left.size == right.size;
  }
};
//----------------------------------------------------------------------------------------------

//__Style to Points Map_________________________________________________________________________
using _style_point_map = std::unordered_multimap<_style, r3_point, _style_hash, _style_equal>;
//----------------------------------------------------------------------------------------------

//__Canvas Names________________________________________________________________________________
std::unordered_map<std::string, size_t> _canvas_names;
//----------------------------------------------------------------------------------------------

//__Create a Unique Name for Canvas_____________________________________________________________
const std::string _make_unique_name(const std::string& name) {
  const auto& search = _canvas_names.find(name);
  if (search != _canvas_names.end()) {
    return name + std::to_string(search->second++);
  } else {
    _canvas_names.insert({name, 1});
    return name;
  }
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Start Plotting Environment__________________________________________________________________
void init() {
  int argc = 1;
  std::array<char*, 1> argv{strdup("app")};
  if (_app == nullptr)
    _app = new TApplication("app", &argc, argv.data());
}
//----------------------------------------------------------------------------------------------

//__End Plotting Environment____________________________________________________________________
void end() {
  _app->Run(true);
}
//----------------------------------------------------------------------------------------------

//__Default Colors______________________________________________________________________________
const color color::WHITE   = {255, 255, 255};
const color color::BLACK   = {  0,   0,   0};
const color color::RED     = {255,   0,   0};
const color color::GREEN   = {  0, 255,   0};
const color color::BLUE    = {  0,   0, 255};
const color color::CYAN    = {  0, 255, 255};
const color color::MAGENTA = {255,   0, 255};
const color color::YELLOW  = {255, 255,   0};
//----------------------------------------------------------------------------------------------

//__Canvas Impl Constructor_____________________________________________________________________
struct canvas::canvas_impl {
  TCanvas* _canvas;
  TView3D* _view;
  std::vector<TPolyLine3D*> _poly_lines;
  _style_point_map _polymarker_map, _polyline_map;
  bool _has_updated = false;

  canvas_impl(const std::string& name, const integer width, const integer height)
      : _canvas(new TCanvas(name.c_str(), name.c_str(), width, height)),
        _view(static_cast<TView3D*>(TView::CreateView())), _poly_lines({}) {
    _view->SetAutoRange(true);
  }

  explicit canvas_impl(const canvas_impl& other) = default;
};
//----------------------------------------------------------------------------------------------

//__Canvas Constructor__________________________________________________________________________
canvas::canvas(const std::string& name, const integer width, const integer height)
    : _impl(std::make_unique<canvas_impl>(_make_unique_name(name), width, height)) {}
//----------------------------------------------------------------------------------------------

//__Canvas Destructor___________________________________________________________________________
canvas::~canvas() = default;
//----------------------------------------------------------------------------------------------

//__Canvas Name_________________________________________________________________________________
const std::string canvas::name() const {
  return _impl->_canvas->GetName();
}
//----------------------------------------------------------------------------------------------

//__Canvas Width________________________________________________________________________________
integer canvas::width() const {
  return _impl->_canvas->GetWindowWidth();
}
//----------------------------------------------------------------------------------------------

//__Canvas Height_______________________________________________________________________________
integer canvas::height() const {
  return _impl->_canvas->GetWindowHeight();
}
//----------------------------------------------------------------------------------------------

//__Draw Canvas_________________________________________________________________________________
void canvas::draw() {
  _impl->_canvas->cd();

  const auto& marker_map = _impl->_polymarker_map;
  const auto&& marker_map_size = marker_map.bucket_count();
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
      polymarker->SetMarkerColor(_color_to_TColor_free_index(style.color));
      polymarker->Draw();
    }
  }

  for (auto poly_line : _impl->_poly_lines)
    poly_line->Draw();


  if (!_impl->_has_updated) {
    _impl->_view->ShowAxis();
    auto axis = TAxis3D::GetPadAxis();
    if (axis) {
      axis->SetLabelColor(kBlack);
      axis->SetAxisColor(kBlack);
    }
  }

  _impl->_canvas->Modified();
  _impl->_canvas->Update();
  _impl->_has_updated = true;
}
//----------------------------------------------------------------------------------------------

//__Clear A Canvas______________________________________________________________________________
void canvas::clear() {
  // TODO: fix
  _impl->_polymarker_map.clear();
  _impl->_polyline_map.clear();
  _impl->_poly_lines.clear();
  _impl->_canvas->GetPad(0)->Clear();
  //_impl->_has_updated = false;
}
//----------------------------------------------------------------------------------------------

//__Add Point to Canvas_________________________________________________________________________
void canvas::add_point(const real x,
                       const real y,
                       const real z,
                       const real size,
                       const color& color) {
  _impl->_polymarker_map.insert({{color, size}, r3_point{x, y, z}});
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
  lines.push_back(new TPolyLine3D());
  auto& line = lines.back();
  line->SetNextPoint(x1, y1, z1);
  line->SetNextPoint(x2, y2, z2);
  line->SetLineWidth(width);
  line->SetLineColor(_color_to_TColor_free_index(color));
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

namespace root { ///////////////////////////////////////////////////////////////////////////////

//__Export Plot Canvas to ROOT File_____________________________________________________________
void to_file(const canvas& canvas,
             const std::string& path) {
  // TODO: implement
}
//----------------------------------------------------------------------------------------------

//__Import Plot Canvas from ROOT File___________________________________________________________
canvas from_file(const std::string& path) {
  // TODO: implement
  return std::move(canvas());
}
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

} /* namespace plot */  ////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

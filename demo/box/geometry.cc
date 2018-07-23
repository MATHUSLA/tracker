/*
 * demo/box/geometry.cc
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

#include "geometry.hh"

namespace MATHUSLA {

namespace box_geometry { ///////////////////////////////////////////////////////////////////////

//__Total Geometry of the Box Detector__________________________________________________________
const geometry::structure_vector& full() {
  static geometry::structure_vector out;
  if (out.empty()) {
    out.reserve(constants::x_total_count * constants::y_total_count * constants::layer_count);
    for (std::size_t z{}; z < constants::layer_count; ++z) {
      const auto z_fullname = std::to_string(1 + z);
      for (std::size_t y{}; y < constants::y_total_count; ++y) {
        const auto y_name = std::to_string(y);
        const auto y_fullname = y < 10UL ? "00" + y_name : (y < 100UL ? "0" + y_name : y_name);
        for (std::size_t x{}; x < constants::x_total_count; ++x) {
          const auto x_name = std::to_string(x);
          out.push_back(z_fullname
            + (x < 10UL ? "00" + x_name : (x < 100UL ? "0" + x_name : x_name))
            + y_fullname);
        }
      }
    }
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Hits per Total Geometry_____________________________________________________________________
type::real event_density(const analysis::full_event& event) {
  return event.size() / static_cast<type::real>(full().size());
}
//----------------------------------------------------------------------------------------------

//__Add Widths to Point_________________________________________________________________________
const analysis::full_hit add_widths(const analysis::hit& point) {
  const auto limits = limits_of_volume(point);
  return analysis::full_hit{
    point.t, limits.center.x, limits.center.y, limits.center.z,
    {time_resolution_of_volume(point),
     limits.max.x - limits.min.x,
     limits.max.y - limits.min.y,
     limits.max.z - limits.min.z}};
}
const analysis::full_event add_widths(const analysis::event& points) {
  analysis::full_event out;
  for (const auto& point : points)
    out.push_back(add_widths(point));
  return out;
}
//----------------------------------------------------------------------------------------------

//__Limits of Volume____________________________________________________________________________
const geometry::box_volume limits_of(const geometry::structure_value& name) {
  namespace cn = constants;
  static const auto half_y_edge = 0.5L * cn::y_edge_length;
  const auto z = std::stoul(name.substr(0, 1));
  const auto x = std::stoul(name.substr(1, 3));
  const auto y = std::stoul(name.substr(4, 3));

  return limits_of_volume(type::r3_point{
    cn::x_displacement + (x+0.5L) * cn::scintillator_x_width,
    cn::y_displacement - half_y_edge + (y+0.5L) * cn::scintillator_y_width,
    -z * cn::layer_spacing - 0.5L * cn::scintillator_height
    });
}
//----------------------------------------------------------------------------------------------

//__Limits of Point_____________________________________________________________________________
const geometry::box_volume limits_of_volume(const type::r3_point point) {
  namespace cn = constants;
  static const auto half_y_edge = 0.5L * cn::y_edge_length;

  const auto local_position = point - type::r3_point{cn::x_displacement, -half_y_edge, cn::steel_height};
  const auto layer_z = local_position.z + cn::steel_height + cn::layer_count * cn::scintillator_height - cn::layer_spacing;
  const auto x_index = static_cast<std::size_t>(std::ceil(local_position.x / cn::scintillator_x_width));
  const auto y_index = static_cast<std::size_t>(std::ceil(local_position.y / cn::scintillator_y_width));
  const auto z_index = static_cast<std::size_t>(std::abs(std::ceil(layer_z / cn::layer_spacing)));

  geometry::box_volume out;
  out.min.x = cn::x_displacement + cn::scintillator_x_width * x_index;
  out.max.x = out.min.x + cn::scintillator_x_width;
  out.min.y = cn::y_displacement - half_y_edge + cn::scintillator_y_width * y_index;
  out.max.y = out.min.y + cn::scintillator_y_width;
  out.max.z = -(cn::scintillator_height + cn::layer_spacing) * z_index;
  out.min.z = out.max.z - cn::scintillator_height;
  out.center = 0.5L * (out.min + out.max);
  return out;
}
const geometry::box_volume limits_of_volume(const type::r4_point point) {
  return limits_of_volume(type::reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

//__Time Resolution of Point____________________________________________________________________
type::real time_resolution_of_volume(const type::r3_point point) {
  return 1.5L;
}
type::real time_resolution_of_volume(const type::r4_point point) {
  return time_resolution_of_volume(type::reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

} /* namespace box_geometry */ /////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

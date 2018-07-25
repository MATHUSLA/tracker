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

#include <iostream>

namespace MATHUSLA {

namespace box { ////////////////////////////////////////////////////////////////////////////////

//__Index Triple Constructor____________________________________________________________________
geometry::index_triple::index_triple(const type::r3_point point) {
  namespace cn = constants;
  const auto local_position = point - type::r3_point{cn::x_displacement, cn::y_displacement, 0.0L};
  x = static_cast<std::size_t>(std::floor(+local_position.x / cn::scintillator_x_width));
  y = static_cast<std::size_t>(std::floor(+local_position.y / cn::scintillator_y_width));
  z = 1UL + static_cast<std::size_t>(std::floor(-local_position.z / (cn::layer_spacing + cn::scintillator_height)));
}
//----------------------------------------------------------------------------------------------

//__Limits of Index Triple Volume_______________________________________________________________
const tracker_geometry::box_volume geometry::index_triple::limits() const {
  namespace cn = constants;
  tracker_geometry::box_volume out;
  out.min.x = cn::x_displacement + cn::scintillator_x_width * x;
  out.max.x = out.min.x + cn::scintillator_x_width;
  out.min.y = cn::y_displacement + cn::scintillator_y_width * y;
  out.max.y = out.min.y + cn::scintillator_y_width;
  out.max.z = -(cn::scintillator_height + cn::layer_spacing) * (z - 1UL);
  out.min.z = out.max.z - cn::scintillator_height;
  out.center = 0.5L * (out.min + out.max);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Name of Index Triple Volume_________________________________________________________________
const tracker_geometry::structure_value geometry::index_triple::name() const {
  const auto x_name = std::to_string(x);
  const auto y_name = std::to_string(y);
  return std::to_string(z)
       + (x < 10UL ? "00" + x_name : (x < 100UL ? "0" + x_name : x_name))
       + (y < 10UL ? "00" + y_name : (y < 100UL ? "0" + y_name : y_name));
}
//----------------------------------------------------------------------------------------------

//__Total Geometry of the Box Detector__________________________________________________________
const tracker_geometry::structure_vector& geometry::full() {
  static tracker_geometry::structure_vector out;
  if (out.empty()) {
    out.reserve(constants::x_total_count * constants::y_total_count * constants::layer_count);
    for (std::size_t z{}; z < constants::layer_count; ++z) {
      const auto z_fullname = std::to_string(1UL + z);
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
type::real geometry::event_density(const analysis::full_event& event) {
  return event.size() / static_cast<type::real>(full().size());
}
//----------------------------------------------------------------------------------------------

//__Get Volume from Point_______________________________________________________________________
const tracker_geometry::structure_value geometry::volume(const type::r3_point point) {
  return index_triple{point}.name();
}
const tracker_geometry::structure_value geometry::volume(const type::r4_point point) {
  return volume(type::reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

//__Limits of Volume____________________________________________________________________________
const tracker_geometry::box_volume geometry::limits_of(const tracker_geometry::structure_value& name) {
  return index_triple{
    std::stoul(name.substr(1, 3)),
    std::stoul(name.substr(4, 3)),
    std::stoul(name.substr(0, 1))}.limits();
}
//----------------------------------------------------------------------------------------------

//__Limits of Point_____________________________________________________________________________
const tracker_geometry::box_volume geometry::limits_of_volume(const type::r3_point point) {
  return index_triple{point}.limits();
}
const tracker_geometry::box_volume geometry::limits_of_volume(const type::r4_point point) {
  return limits_of_volume(type::reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

//__Default Time Resolution_____________________________________________________________________
type::real geometry::default_time_resolution() {
  return constants::scintillator_time_resolution;
}
//----------------------------------------------------------------------------------------------

//__Time Resolution of Volume___________________________________________________________________
type::real geometry::time_resolution_of(const tracker_geometry::structure_value& name) {
  return default_time_resolution();
}
//----------------------------------------------------------------------------------------------

//__Time Resolution of Point____________________________________________________________________
type::real geometry::time_resolution_of_volume(const type::r3_point point) {
  return default_time_resolution();
}
type::real geometry::time_resolution_of_volume(const type::r4_point point) {
  return time_resolution_of_volume(type::reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

} /* namespace box */ //////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

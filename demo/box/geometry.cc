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

#include <tracker/util/string.hh>

namespace MATHUSLA {

namespace box { ////////////////////////////////////////////////////////////////////////////////

//__Default Geometry Components_________________________________________________________________
std::size_t geometry::layer_count = constants::layer_count;
type::real geometry::scintillator_x_width = constants::scintillator_x_width;
type::real geometry::scintillator_y_width = constants::scintillator_y_width;
type::real geometry::scintillator_height = constants::scintillator_height;
type::real geometry::layer_spacing = constants::layer_spacing;
type::real geometry::x_displacement = constants::x_displacement;
type::real geometry::y_displacement = constants::y_displacement;
type::real geometry::x_edge_length = constants::x_edge_length;
type::real geometry::y_edge_length = constants::y_edge_length;
//----------------------------------------------------------------------------------------------

//__Total Scintillator Count in X Direction_____________________________________________________
type::real geometry::x_total_count() {
  return static_cast<std::size_t>(std::ceil(x_edge_length / scintillator_x_width));
}
//----------------------------------------------------------------------------------------------

//__Total Scintillator Count in Y Direction_____________________________________________________
type::real geometry::y_total_count() {
  return static_cast<std::size_t>(std::ceil(y_edge_length / scintillator_y_width));
}
//----------------------------------------------------------------------------------------------

//__Total Scintillator Count____________________________________________________________________
type::real geometry::total_count() {
  return x_total_count() * y_total_count() * layer_count;
}
//----------------------------------------------------------------------------------------------

//__Index Triple Constructor____________________________________________________________________
geometry::index_triple::index_triple(const type::r3_point point) {
  const auto local_position = point - type::r3_point{x_displacement, y_displacement, 0.0L};
  x = static_cast<std::size_t>(std::floor(+local_position.x / scintillator_x_width));
  y = static_cast<std::size_t>(std::floor(+local_position.y / scintillator_y_width));
  z = 1UL + static_cast<std::size_t>(std::floor(-local_position.z / (layer_spacing + scintillator_height)));
}
//----------------------------------------------------------------------------------------------

//__Index Triple Constructor____________________________________________________________________
geometry::index_triple::index_triple(const std::string& name,
                                     const std::string& delimeter) {
  util::string_vector tokens;
  util::string::split(name, tokens, delimeter);
  x = std::stoul(tokens[1]);
  y = std::stoul(tokens[2]);
  z = std::stoul(tokens[0]);
}
//----------------------------------------------------------------------------------------------

//__Limits of Index Triple Volume_______________________________________________________________
const tracker_geometry::box_volume geometry::index_triple::limits() const {
  tracker_geometry::box_volume out;
  out.min.x = x_displacement + scintillator_x_width * x;
  out.max.x = out.min.x + scintillator_x_width;
  out.min.y = y_displacement + scintillator_y_width * y;
  out.max.y = out.min.y + scintillator_y_width;
  out.max.z = -(scintillator_height + layer_spacing) * (z - 1UL);
  out.min.z = out.max.z - scintillator_height;
  out.center = 0.5L * (out.min + out.max);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Name of Index Triple Volume_________________________________________________________________
const tracker_geometry::structure_value geometry::index_triple::name() const {
  return std::to_string(z) + "_" + std::to_string(x) + "_" + std::to_string(y);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Insert Volumes into Layers__________________________________________________________________
void _insert_volumes(const tracker_geometry::structure_value& prefix,
                     const tracker_geometry::structure_value& suffix,
                     const std::size_t count,
                     tracker_geometry::structure_vector& out) {
  for (std::size_t x{}; x < count; ++x)
    out.push_back(prefix + std::to_string(x) + suffix);
}
//----------------------------------------------------------------------------------------------

//__Create Scintillator Layers__________________________________________________________________
const tracker_geometry::structure_vector _make_layer(const std::size_t layer_index,
                                                     const std::size_t x_count,
                                                     const std::size_t y_count) {
  tracker_geometry::structure_vector out;
  out.reserve(x_count * y_count);
  const auto layer_string = std::to_string(1UL + layer_index);
  for (std::size_t y{}; y < y_count; ++y)
    _insert_volumes(layer_string + "_", "_" + std::to_string(y), x_count, out);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Add Layers for Full Geometry________________________________________________________________
void _add_layers(const std::size_t begin,
                 const std::size_t end,
                 const std::size_t x_count,
                 const std::size_t y_count,
                 std::vector<tracker_geometry::structure_vector>& layers) {
  for (std::size_t z{}; z < end; ++z)
    layers.push_back(_make_layer(z + begin, x_count, y_count));
}
//----------------------------------------------------------------------------------------------

//__Load Layers into Full Structure_____________________________________________________________
void _load_layers(const std::size_t begin,
                  const std::size_t end,
                  std::vector<tracker_geometry::structure_vector>& layers,
                  tracker_geometry::structure_vector& out) {
  for (std::size_t i{begin}; i < end; ++i)
    for (const auto& volume : layers[i])
      out.push_back(volume);
}
//----------------------------------------------------------------------------------------------

//__Unload Layers from Full Structure___________________________________________________________
void _unload_layers(const std::size_t count,
                    const std::size_t x_count,
                    const std::size_t y_count,
                    tracker_geometry::structure_vector& out) {
  out.erase(out.cend() - count * x_count * y_count, out.cend());
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Total Geometry of the Box Detector__________________________________________________________
const tracker_geometry::structure_vector& geometry::full(const std::size_t count) {
  static tracker_geometry::structure_vector out;
  static std::vector<tracker_geometry::structure_vector> layers;
  static std::size_t current_layer_count{}, current_x_count{}, current_y_count{};

  if (current_x_count != x_total_count() || current_y_count != y_total_count()) {
    current_x_count = x_total_count();
    current_y_count = y_total_count();
    layers.clear();
    out.clear();
  }

  if (current_layer_count < count) {
    const auto layers_size = layers.size();
    if (layers_size < count)
      _add_layers(layers_size, count - layers_size, current_x_count, current_y_count, layers);
    _load_layers(current_layer_count, count - current_layer_count, layers, out);
  } else if (current_layer_count > count) {
    _unload_layers(current_layer_count - count, current_x_count, current_y_count, out);
  }

  current_layer_count = count;
  return out;
}
//----------------------------------------------------------------------------------------------

//__Hits per Total Geometry_____________________________________________________________________
type::real geometry::event_density(const analysis::full_event& event) {
  return event.size() / static_cast<type::real>(total_count());
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
  return index_triple{name}.limits();
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

//__Restrict Number of Layers in Geometry_______________________________________________________
const analysis::full_event geometry::restrict_layer_count(const analysis::full_event& event,
                                                          const std::size_t layers) {
  analysis::full_event out;
  util::algorithm::back_insert_copy_if(event, out, [&](const auto& hit) {
    const auto name = volume(type::reduce_to_r3(hit));
    return std::stoul(name.substr(0, name.find_first_of("_"))) <= layers;
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Vector of Value Tags of Geometrical Quantities______________________________________________
const plot::value_tag_vector geometry::value_tags() {
  return plot::value_tag_vector{
    {"LAYER_COUNT",          std::to_string(layer_count)},
    {"SCINTILLATOR_X_WIDTH", std::to_string(scintillator_x_width / units::length) + " " + units::length_string},
    {"SCINTILLATOR_Y_WIDTH", std::to_string(scintillator_y_width / units::length) + " " + units::length_string},
    {"SCINTILLATOR_HEIGHT",  std::to_string(scintillator_height  / units::length) + " " + units::length_string},
    {"LAYER_SPACING",        std::to_string(layer_spacing        / units::length) + " " + units::length_string},
    {"X_DISPLACEMENT",       std::to_string(x_displacement       / units::length) + " " + units::length_string},
    {"Y_DISPLACEMENT",       std::to_string(y_displacement       / units::length) + " " + units::length_string},
    {"X_EDGE_LENGTH",        std::to_string(x_edge_length        / units::length) + " " + units::length_string},
    {"Y_EDGE_LENGTH",        std::to_string(y_edge_length        / units::length) + " " + units::length_string}
  };
}
//----------------------------------------------------------------------------------------------

} /* namespace box */ //////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

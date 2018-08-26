/*
 * demo/box/geometry.hh
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

#ifndef TRACKER__BOX__GEOMETRY_HH
#define TRACKER__BOX__GEOMETRY_HH
#pragma once

#include <tracker/analysis/type.hh>
#include <tracker/geometry.hh>
#include <tracker/plot.hh>

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace tracker_geometry = MATHUSLA::TRACKER::geometry;
namespace plot = MATHUSLA::TRACKER::plot;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

namespace box { ////////////////////////////////////////////////////////////////////////////////

namespace constants { //////////////////////////////////////////////////////////////////////////

static const auto scintillator_time_resolution  = 1.5L*units::ns;

static const auto x_edge_length                 = 100.00L*units::m;
static const auto y_edge_length                 = 100.00L*units::m;
static const auto x_displacement                = 100.00L*units::m;
static const auto y_displacement                = -50.00L*units::m;
static const auto steel_height                  =   3.00L*units::cm;
static const auto air_gap                       =  20.00L*units::m;
static const auto scintillator_x_width          =   0.25L*units::m;
static const auto scintillator_y_width          =   0.25L*units::m;
static const auto scintillator_height           =   1.00L*units::cm;
static const auto scintillator_casing_thickness =   0.10L*units::cm;

static const auto layer_spacing                 = 1.50L*units::m;
static const auto layer_count                   = 5UL;
static const auto x_total_count                 = static_cast<std::size_t>(std::ceil(x_edge_length / scintillator_x_width));
static const auto y_total_count                 = static_cast<std::size_t>(std::ceil(y_edge_length / scintillator_y_width));
static const auto total_count                   = x_total_count * y_total_count * layer_count;

static const auto half_x_edge_length            = 0.5L * x_edge_length;
static const auto half_y_edge_length            = 0.5L * y_edge_length;
static const auto full_detector_height          = steel_height + layer_count * (layer_spacing + scintillator_height) - layer_spacing;
static const auto half_detector_height          = 0.5L * full_detector_height;

} /* namespace constants */ ////////////////////////////////////////////////////////////////////

//__Geometry Structure__________________________________________________________________________
struct geometry {
  static std::size_t layer_count;
  static type::real scintillator_x_width;
  static type::real scintillator_y_width;
  static type::real scintillator_height;
  static type::real layer_spacing;
  static type::real x_displacement;
  static type::real y_displacement;
  static type::real x_edge_length;
  static type::real y_edge_length;

  static type::real x_total_count();
  static type::real y_total_count();
  static type::real total_count();

  template<class Geometry>
  static void import(const Geometry& g) {
    layer_count = g.layer_count;
    scintillator_x_width = g.scintillator_x_width;
    scintillator_y_width = g.scintillator_y_width;
    scintillator_height = g.scintillator_height;
    layer_spacing = g.layer_spacing;
    x_displacement = g.x_displacement;
    y_displacement = g.y_displacement;
    x_edge_length = g.x_edge_length;
    y_edge_length = g.y_edge_length;
  }

  struct index_triple {
    std::size_t x, y, z;
    index_triple() = default;
    index_triple(std::size_t x_index,
                 std::size_t y_index,
                 std::size_t z_index) : x(x_index), y(y_index), z(z_index) {}
    index_triple(const type::r3_point point);
    index_triple(const std::string& name,
                 const std::string& delimeter="_");
    const tracker_geometry::box_volume limits() const;
    const tracker_geometry::structure_value name() const;
  };

  static const tracker_geometry::structure_vector& full(const std::size_t count=box::constants::layer_count);
  static type::real event_density(const analysis::full_event& event);
  static const tracker_geometry::structure_value volume(const type::r3_point point);
  static const tracker_geometry::structure_value volume(const type::r4_point point);
  static const tracker_geometry::box_volume limits_of(const tracker_geometry::structure_value& name);
  static const tracker_geometry::box_volume limits_of_volume(const type::r3_point point);
  static const tracker_geometry::box_volume limits_of_volume(const type::r4_point point);
  static type::real default_time_resolution();
  static type::real time_resolution_of(const tracker_geometry::structure_value& name);
  static type::real time_resolution_of_volume(const type::r3_point point);
  static type::real time_resolution_of_volume(const type::r4_point point);

  static const analysis::full_event restrict_layer_count(const analysis::full_event& event,
                                                         const std::size_t layers);
  static const plot::value_tag_vector value_tags();
};
//----------------------------------------------------------------------------------------------

} /* namespace box */ //////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* TRACKER__BOX__GEOMETRY_HH */

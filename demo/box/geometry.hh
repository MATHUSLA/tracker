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

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace geometry = MATHUSLA::TRACKER::geometry;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

namespace box_geometry { ///////////////////////////////////////////////////////////////////////

namespace constants { //////////////////////////////////////////////////////////////////////////

static const auto x_edge_length                 = 100*units::m;
static const auto y_edge_length                 = 100*units::m;
static const auto x_displacement                = 100*units::m;
static const auto y_displacement                = 0;
static const auto steel_height                  = 3.0L*units::cm;
static const auto air_gap                       = 20*units::m;
static const auto scintillator_x_width          = 0.25L*units::m;
static const auto scintillator_y_width          = 0.25L*units::m;
static const auto scintillator_height           = 1*units::cm;
static const auto scintillator_casing_thickness = 0.1*units::cm;
static const auto layer_spacing                 = 1.5L*units::m;
static const auto layer_count                   = 5UL;
static const auto x_total_count                 = static_cast<std::size_t>(std::ceil(x_edge_length / scintillator_x_width));
static const auto y_total_count                 = static_cast<std::size_t>(std::ceil(y_edge_length / scintillator_y_width));
static const auto full_detector_height          = steel_height + layer_count * (layer_spacing + scintillator_height);
static const auto half_detector_height          = 0.5L * full_detector_height;

} /* namespace constants */ ////////////////////////////////////////////////////////////////////

//__Total Geometry of the Box Detector__________________________________________________________
const geometry::structure_vector& full();
//----------------------------------------------------------------------------------------------

//__Hits per Total Geometry_____________________________________________________________________
type::real event_density(const analysis::full_event& event);
//----------------------------------------------------------------------------------------------

//__Add Widths to Point_________________________________________________________________________
const analysis::full_hit add_widths(const analysis::hit& point);
const analysis::full_event add_widths(const analysis::event& points);
//----------------------------------------------------------------------------------------------

//__Limits of Volume____________________________________________________________________________
const geometry::box_volume limits_of(const geometry::structure_value& name);
//----------------------------------------------------------------------------------------------

//__Limits of Point_____________________________________________________________________________
const geometry::box_volume limits_of_volume(const type::r3_point point);
const geometry::box_volume limits_of_volume(const type::r4_point point);
//----------------------------------------------------------------------------------------------

//__Time Resolution of Point____________________________________________________________________
type::real time_resolution_of_volume(const type::r3_point point);
type::real time_resolution_of_volume(const type::r4_point point);
//----------------------------------------------------------------------------------------------

} /* namespace box_geometry */ /////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* TRACKER__BOX__GEOMETRY_HH */

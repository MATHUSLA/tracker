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

//__Total Geometry of the Box Detector__________________________________________________________
const geometry::structure_vector box_geometry();
//----------------------------------------------------------------------------------------------

//__Hits per Total Geometry_____________________________________________________________________
type::real modified_geometry_event_density(const analysis::event& event);
//----------------------------------------------------------------------------------------------

//__Limits of Point_____________________________________________________________________________
const geometry::box_volume limits_of(const type::r3_point point);
const geometry::box_volume limits_of(const type::r4_point point);
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

#endif /* TRACKER__BOX__GEOMETRY_HH */

/*
 * include/util/geometry.hh
 *
 * Copyright 2018 Brandon Gomes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef TRACKER__GEOMETRY_HH
#define TRACKER__GEOMETRY_HH
#pragma once

#include "point.hh"

namespace MATHUSLA { namespace TRACKER {

namespace geometry { ///////////////////////////////////////////////////////////////////////////

using namespace type;

//__Geometry Navigation System__________________________________________________________________
void open(const std::string& path);
void close();
//----------------------------------------------------------------------------------------------

//__Volume Containment Check____________________________________________________________________
bool is_inside_volume(const r3_point& point, const std::string& name);
bool is_inside_volume(const r4_point& point, const std::string& name);
//----------------------------------------------------------------------------------------------

//__Volume Hierarchy Search_____________________________________________________________________
const std::vector<std::string> volume_hierarchy(const r3_point& point);
const std::vector<std::string> volume_hierarchy(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Volume Search_______________________________________________________________________________
const std::string volume(const r3_point& point);
const std::string volume(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Box Volume__________________________________________________________________________________
struct box_volume { r3_point center, min, max; };
//----------------------------------------------------------------------------------------------

//__Limit Box of a Volume_______________________________________________________________________
const box_volume limits_of(const std::string& name);
//----------------------------------------------------------------------------------------------

//__Limit Box of the Volume with Point__________________________________________________________
const box_volume limits_of_volume(const r3_point& point);
const box_volume limits_of_volume(const r4_point& point);
//----------------------------------------------------------------------------------------------

} /* namespace geometry */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__GEOMETRY_HH */

/*
 * include/tracker/geometry.hh
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

#ifndef TRACKER__GEOMETRY_HH
#define TRACKER__GEOMETRY_HH
#pragma once

#include <unordered_map>

#include <tracker/core/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace geometry { ///////////////////////////////////////////////////////////////////////////

using namespace type;

//__Geometry Structure Types____________________________________________________________________
using structure_value = std::string;
using structure_vector = std::vector<structure_value>;
//----------------------------------------------------------------------------------------------

//__Geometry Detector Map_______________________________________________________________________
using detector_map = std::unordered_map<integer, structure_value>;
//----------------------------------------------------------------------------------------------

//__Geometry Time Resolution Map________________________________________________________________
using time_resolution_map = std::unordered_map<structure_value, real>;
//----------------------------------------------------------------------------------------------

//__Geometry Navigation System__________________________________________________________________
void open(const std::string& path,
          const real default_time_error,
          const time_resolution_map& map={});
void close();
//----------------------------------------------------------------------------------------------

//__Current Geometry File_______________________________________________________________________
const std::string& current_geometry_path();
//----------------------------------------------------------------------------------------------

//__List of all Geometry Volumes________________________________________________________________
const structure_vector full_structure();
//----------------------------------------------------------------------------------------------

//__List of all Geometry Volumes Except those in the List_______________________________________
const structure_vector full_structure_except(const structure_vector& names);
//----------------------------------------------------------------------------------------------

//__Get Default Time Resolution for Detector Volumes____________________________________________
real default_time_resolution();
//----------------------------------------------------------------------------------------------

//__Volume Hierarchy of a Volume________________________________________________________________
// TODO: finish
// const std::vector<std::string> volume_hierarchy(const std::string& name);
//----------------------------------------------------------------------------------------------

//__Volume Containment Check____________________________________________________________________
bool is_inside_volume(const r3_point& point,
                      const structure_value& name);
bool is_inside_volume(const r4_point& point,
                      const structure_value& name);
//----------------------------------------------------------------------------------------------

//__Volume Hierarchy Search_____________________________________________________________________
const structure_vector volume_hierarchy(const r3_point& point);
const structure_vector volume_hierarchy(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Volume Search_______________________________________________________________________________
const structure_value volume(const r3_point& point);
const structure_value volume(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Box Volume__________________________________________________________________________________
struct box_volume { r3_point center, min, max; };
using box_volume_vector = std::vector<box_volume>;
//----------------------------------------------------------------------------------------------

//__Box Volume Stream Overload__________________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const box_volume& box);
//----------------------------------------------------------------------------------------------

//__Limit Box of a Volume_______________________________________________________________________
const box_volume limits_of(const structure_value& name);
//----------------------------------------------------------------------------------------------

//__Time Error on a Detector Component__________________________________________________________
real time_resolution_of(const structure_value& name);
//----------------------------------------------------------------------------------------------

//__Compute Intersection of Box Volumes_________________________________________________________
const box_volume coordinatewise_intersection(const box_volume& first,
                                             const box_volume& second);
//----------------------------------------------------------------------------------------------

//__Compute Union of Box Volumes________________________________________________________________
const box_volume coordinatewise_union(const box_volume& first,
                                      const box_volume& second);
//----------------------------------------------------------------------------------------------

//__Limit Box of the Volume with Point__________________________________________________________
const box_volume limits_of_volume(const r3_point& point);
const box_volume limits_of_volume(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Time Error on a Detector Component from a Point_____________________________________________
real time_resolution_of_volume(const r3_point& point);
real time_resolution_of_volume(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Box Volume Containment Check________________________________________________________________
bool is_inside_volume(const r3_point& point,
                      const box_volume& box);
bool is_inside_volume(const r4_point& point,
                      const box_volume& box);
//----------------------------------------------------------------------------------------------

//__Find Center of Geometry_____________________________________________________________________
const r3_point find_center(const structure_value& name);
const r3_point find_center(const r3_point& point);
const r4_point find_center(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Get Volumes from Set of Points______________________________________________________________
const structure_vector volumes(const r3_point_vector& points);
const structure_vector volumes(const r4_point_vector& points);
//----------------------------------------------------------------------------------------------

//__Get Unique Volumes from Set of Points_______________________________________________________
const structure_vector unique_volumes(const r3_point_vector& points);
const structure_vector unique_volumes(const r4_point_vector& points);
//----------------------------------------------------------------------------------------------

//__Add Volume to Geometry Using Local Coordinates______________________________________________
std::size_t add_to_volume_local(const structure_value& parent,
                                const structure_value& name,
                                const box_volume& box);
//----------------------------------------------------------------------------------------------

//__Add Volumes to Geometry Using Local Coordinates_____________________________________________
std::size_t add_to_volume_local(const structure_value& parent,
                                const structure_vector& names,
                                const box_volume_vector& boxes);
//----------------------------------------------------------------------------------------------

//__Add Volume to Geometry Using Global Coordinates_____________________________________________
std::size_t add_to_volume_global(const structure_value& parent,
                                 const structure_value& name,
                                 const box_volume& box);
//----------------------------------------------------------------------------------------------

//__Add Volumes to Geometry Using Global Coordinates____________________________________________
std::size_t add_to_volume_global(const structure_value& parent,
                                 const structure_vector& names,
                                 const box_volume_vector& boxes);
//----------------------------------------------------------------------------------------------

namespace custom { /////////////////////////////////////////////////////////////////////////////

//__Current Geometry File_______________________________________________________________________
template<class Geometry=void>
const std::string& current_geometry_path() {
  return Geometry::current_geometry_path();
}
template<>
inline const std::string& current_geometry_path<void>() {
  return geometry::current_geometry_path();
}
//----------------------------------------------------------------------------------------------

//__List of all Geometry Volumes________________________________________________________________
template<class Geometry=void>
const structure_vector full_structure() {
  return Geometry::full_structure();
}
template<>
inline const structure_vector full_structure<void>() {
  return geometry::full_structure();
}
//----------------------------------------------------------------------------------------------

//__List of all Geometry Volumes Except those in the List_______________________________________
template<class Geometry=void>
const structure_vector full_structure_except(const structure_vector& names) {
  return Geometry::full_structure_except(names);
}
template<>
inline const structure_vector full_structure_except<void>(const structure_vector& names) {
  return geometry::full_structure_except(names);
}
//----------------------------------------------------------------------------------------------

//__Get Default Time Resolution for Detector Volumes____________________________________________
template<class Geometry=void>
real default_time_resolution() {
  return Geometry::default_time_resolution();
}
template<>
inline real default_time_resolution<void>() {
  return geometry::default_time_resolution();
}
//----------------------------------------------------------------------------------------------

//__Volume Hierarchy of a Volume________________________________________________________________
// TODO: finish
// template<class Geometry>
// const std::vector<std::string> volume_hierarchy(const std::string& name);
//----------------------------------------------------------------------------------------------

//__Volume Containment Check____________________________________________________________________
template<class Geometry=void>
bool is_inside_volume(const r3_point& point,
                      const structure_value& name) {
  return Geometry::is_inside_volume(point, name);
}
template<>
inline bool is_inside_volume<void>(const r3_point& point,
                                   const structure_value& name) {
  return geometry::is_inside_volume(point, name);
}
template<class Geometry=void>
bool is_inside_volume(const r4_point& point,
                      const structure_value& name) {
  return Geometry::is_inside_volume(point, name);
}
template<>
inline bool is_inside_volume<void>(const r4_point& point,
                                   const structure_value& name) {
  return geometry::is_inside_volume(point, name);
}
//----------------------------------------------------------------------------------------------

//__Volume Hierarchy Search_____________________________________________________________________
template<class Geometry=void>
const structure_vector volume_hierarchy(const r3_point& point) {
  return Geometry::volume_hierarchy(point);
}
template<>
inline const structure_vector volume_hierarchy<void>(const r3_point& point) {
  return geometry::volume_hierarchy(point);
}
template<class Geometry=void>
const structure_vector volume_hierarchy(const r4_point& point) {
  return Geometry::volume_hierarchy(point);
}
template<>
inline const structure_vector volume_hierarchy<void>(const r4_point& point) {
  return geometry::volume_hierarchy(point);
}
//----------------------------------------------------------------------------------------------

//__Volume Search_______________________________________________________________________________
template<class Geometry=void>
const structure_value volume(const r3_point& point) {
  return Geometry::volume(point);
}
template<>
inline const structure_value volume<void>(const r3_point& point) {
  return geometry::volume(point);
}
template<class Geometry=void>
const structure_value volume(const r4_point& point) {
  return Geometry::volume(point);
}
template<>
inline const structure_value volume<void>(const r4_point& point) {
  return geometry::volume(point);
}
//----------------------------------------------------------------------------------------------

//__Limit Box of a Volume_______________________________________________________________________
template<class Geometry=void>
const box_volume limits_of(const structure_value& name) {
  return Geometry::limits_of(name);
}
template<>
inline const box_volume limits_of<void>(const structure_value& name) {
  return geometry::limits_of(name);
}
//----------------------------------------------------------------------------------------------

//__Time Error on a Detector Component__________________________________________________________
template<class Geometry=void>
real time_resolution_of(const structure_value& name) {
  return Geometry::time_resolution_of(name);
}
template<>
inline real time_resolution_of<void>(const structure_value& name) {
  return geometry::time_resolution_of(name);
}
//----------------------------------------------------------------------------------------------

//__Compute Intersection of Box Volumes_________________________________________________________
template<class Geometry=void>
const box_volume coordinatewise_intersection(const box_volume& first,
                                             const box_volume& second) {
  return Geometry::coordinatewise_intersection(first, second);
}
template<>
inline const box_volume coordinatewise_intersection<void>(const box_volume& first,
                                                          const box_volume& second) {
  return geometry::coordinatewise_intersection(first, second);
}
//----------------------------------------------------------------------------------------------

//__Compute Union of Box Volumes________________________________________________________________
template<class Geometry=void>
const box_volume coordinatewise_union(const box_volume& first,
                                      const box_volume& second) {
  return Geometry::coordinatewise_intersection(first, second);
}
template<>
inline const box_volume coordinatewise_union<void>(const box_volume& first,
                                                   const box_volume& second) {
  return geometry::coordinatewise_intersection(first, second);
}
//----------------------------------------------------------------------------------------------

//__Limit Box of the Volume with Point__________________________________________________________
template<class Geometry=void>
const box_volume limits_of_volume(const r3_point& point) {
  return Geometry::limits_of_volume(point);
}
template<>
inline const box_volume limits_of_volume<void>(const r3_point& point) {
  return geometry::limits_of_volume(point);
}
template<class Geometry=void>
const box_volume limits_of_volume(const r4_point& point) {
  return Geometry::limits_of_volume(point);
}
template<>
inline const box_volume limits_of_volume<void>(const r4_point& point) {
  return geometry::limits_of_volume(point);
}
//----------------------------------------------------------------------------------------------

//__Time Error on a Detector Component from a Point_____________________________________________
template<class Geometry=void>
real time_resolution_of_volume(const r3_point& point) {
  return Geometry::time_resolution_of_volume(point);
}
template<>
inline real time_resolution_of_volume<void>(const r3_point& point) {
  return geometry::time_resolution_of_volume(point);
}
template<class Geometry=void>
real time_resolution_of_volume(const r4_point& point) {
  return Geometry::time_resolution_of_volume(point);
}
template<>
inline real time_resolution_of_volume<void>(const r4_point& point) {
  return geometry::time_resolution_of_volume(point);
}
//----------------------------------------------------------------------------------------------

//__Box Volume Containment Check________________________________________________________________
template<class Geometry=void>
bool is_inside_volume(const r3_point& point,
                      const box_volume& box) {
  return Geometry::is_inside_volume(point, box);
}
template<>
inline bool is_inside_volume<void>(const r3_point& point,
                                   const box_volume& box) {
  return geometry::is_inside_volume(point, box);
}
template<class Geometry=void>
bool is_inside_volume(const r4_point& point,
                      const box_volume& box) {
  return Geometry::is_inside_volume(point, box);
}
template<>
inline bool is_inside_volume<void>(const r4_point& point,
                                   const box_volume& box) {
  return geometry::is_inside_volume(point, box);
}
//----------------------------------------------------------------------------------------------

//__Find Center of Geometry_____________________________________________________________________
template<class Geometry=void>
const r3_point find_center(const structure_value& name) {
  return Geometry::find_center(name);
}
template<>
inline const r3_point find_center<void>(const structure_value& name) {
  return geometry::find_center(name);
}
template<class Geometry=void>
const r3_point find_center(const r3_point& point) {
  return Geometry::find_center(point);
}
template<>
inline const r3_point find_center<void>(const r3_point& point) {
  return geometry::find_center(point);
}
template<class Geometry=void>
const r4_point find_center(const r4_point& point) {
  return Geometry::find_center(point);
}
template<>
inline const r4_point find_center<void>(const r4_point& point) {
  return geometry::find_center(point);
}
//----------------------------------------------------------------------------------------------

//__Get Volumes from Set of Points______________________________________________________________
template<class Geometry=void>
const structure_vector volumes(const r3_point_vector& points) {
  return Geometry::volumes(points);
}
template<>
inline const structure_vector volumes<void>(const r3_point_vector& points) {
  return geometry::volumes(points);
}
template<class Geometry=void>
const structure_vector volumes(const r4_point_vector& points) {
  return Geometry::volumes(points);
}
template<>
inline const structure_vector volumes<void>(const r4_point_vector& points) {
  return geometry::volumes(points);
}
//----------------------------------------------------------------------------------------------

//__Get Unique Volumes from Set of Points_______________________________________________________
template<class Geometry=void>
const structure_vector unique_volumes(const r3_point_vector& points) {
  return Geometry::unique_volumes(points);
}
template<>
inline const structure_vector unique_volumes<void>(const r3_point_vector& points) {
  return geometry::unique_volumes(points);
}
template<class Geometry=void>
const structure_vector unique_volumes(const r4_point_vector& points) {
  return Geometry::unique_volumes(points);
}
template<>
inline const structure_vector unique_volumes<void>(const r4_point_vector& points) {
  return geometry::unique_volumes(points);
}
//----------------------------------------------------------------------------------------------

//__Add Volume to Geometry Using Local Coordinates______________________________________________
template<class Geometry=void>
std::size_t add_to_volume_local(const structure_value& parent,
                                const structure_value& name,
                                const box_volume& box) {
  return Geometry::add_to_volume_local(parent, name, box);
}
template<>
inline std::size_t add_to_volume_local<void>(const structure_value& parent,
                                             const structure_value& name,
                                             const box_volume& box) {
  return geometry::add_to_volume_local(parent, name, box);
}
//----------------------------------------------------------------------------------------------

//__Add Volumes to Geometry Using Local Coordinates_____________________________________________
template<class Geometry=void>
std::size_t add_to_volume_local(const structure_value& parent,
                                const structure_vector& names,
                                const box_volume_vector& boxes) {
  return Geometry::add_to_volume_local(parent, names, boxes);
}
template<>
inline std::size_t add_to_volume_local<void>(const structure_value& parent,
                                             const structure_vector& names,
                                             const box_volume_vector& boxes) {
  return geometry::add_to_volume_local(parent, names, boxes);
}
//----------------------------------------------------------------------------------------------

//__Add Volume to Geometry Using Global Coordinates_____________________________________________
template<class Geometry=void>
std::size_t add_to_volume_global(const structure_value& parent,
                                 const structure_value& name,
                                 const box_volume& box) {
  return Geometry::add_to_volume_global(parent, name, box);
}
template<>
inline std::size_t add_to_volume_global<void>(const structure_value& parent,
                                              const structure_value& name,
                                              const box_volume& box) {
  return geometry::add_to_volume_global(parent, name, box);
}
//----------------------------------------------------------------------------------------------

//__Add Volumes to Geometry Using Global Coordinates____________________________________________
template<class Geometry=void>
std::size_t add_to_volume_global(const structure_value& parent,
                                 const structure_vector& names,
                                 const box_volume_vector& boxes) {
  return Geometry::add_to_volume_global(parent, names, boxes);
}
template<>
inline std::size_t add_to_volume_global<void>(const structure_value& parent,
                                              const structure_vector& names,
                                              const box_volume_vector& boxes) {
  return geometry::add_to_volume_global(parent, names, boxes);
}
//----------------------------------------------------------------------------------------------

} /* namespace custom */ ///////////////////////////////////////////////////////////////////////

} /* namespace geometry */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__GEOMETRY_HH */

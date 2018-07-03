/*
 * src/tracker/geometry.cc
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

#include <tracker/geometry.hh>

#include <unordered_map>

#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VUserDetectorConstruction.hh>
#include <Geant4/G4RunManager.hh>
#include <Geant4/G4VoxelLimits.hh>
#include <Geant4/G4AffineTransform.hh>
// TODO: use #include <Geant4/G4IntersectionSolid.hh>
// TODO: use #include <Geant4/G4UnionSolid.hh>

#include <tracker/core/units.hh>

#include <tracker/util/algorithm.hh>

namespace MATHUSLA { namespace TRACKER {

namespace geometry { ///////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Geometric Type______________________________________________________________________________
struct _geometric_volume { G4VPhysicalVolume* volume; G4AffineTransform transform; };
//----------------------------------------------------------------------------------------------

//__Hash of R3 Point for Precomputed Name Map___________________________________________________
struct r3_point_hash {
  size_t operator()(const r3_point& point) const {
    auto out = 17U;
    out = std::fma(out, 31U, std::hash<real>()(point.x));
    out = std::fma(out, 31U, std::hash<real>()(point.y));
    out = std::fma(out, 31U, std::hash<real>()(point.z));
    return out;
  }
};
//----------------------------------------------------------------------------------------------

//__Static Geometry Variables___________________________________________________________________
thread_local std::string _path;
thread_local G4RunManager* _manager;
thread_local G4VPhysicalVolume* _world;
thread_local std::unordered_map<structure_value, _geometric_volume> _geometry;
thread_local std::vector<structure_value> _geometry_insertion_order;
thread_local real _default_time_resolution = 2 * units::time;
thread_local time_resolution_map _time_resolution_map;
thread_local std::unordered_map<r3_point, structure_value, r3_point_hash> _name_cache;
// TODO: use -> thread_local std::unordered_map<std::string, box_volume> _box_cache;
const G4VoxelLimits _blank_voxels;
const G4AffineTransform _blank_transform;
//----------------------------------------------------------------------------------------------

//__Convert From G4ThreeVector to R3_Point______________________________________________________
const r3_point _to_r3_point(const G4ThreeVector& vector) {
  return { vector.x(), vector.y(), vector.z() };
}
//----------------------------------------------------------------------------------------------

//__Convert From R3_Point to G4ThreeVector______________________________________________________
const G4ThreeVector _to_G4ThreeVector(const r3_point& point) {
  return G4ThreeVector(point.x, point.y, point.z);
}
//----------------------------------------------------------------------------------------------

//__Safely Get Name of Volume___________________________________________________________________
const structure_value _get_name(const _geometric_volume& gvolume) {
  const auto& volume = gvolume.volume;
  return volume ? volume->GetName() : "";
}
//----------------------------------------------------------------------------------------------

//__Transform to Local Coordinates______________________________________________________________
const G4ThreeVector _to_local_transform(const G4ThreeVector& vector,
                                        const G4ThreeVector& translation,
                                        const G4RotationMatrix& rotation) {
  return (rotation.inverse() * vector) - translation;
}
//----------------------------------------------------------------------------------------------

//__Transform Vector by Translation then Rotation_______________________________________________
const G4ThreeVector _to_global_transform(const G4ThreeVector& vector,
                                         const G4ThreeVector& translation,
                                         const G4RotationMatrix& rotation) {
  return rotation * (vector + translation);
}
//----------------------------------------------------------------------------------------------

//__Load Geometry from GDML File________________________________________________________________
G4VPhysicalVolume* _load_geometry(const std::string& path) {
  static G4ThreadLocal G4GDMLParser _parser;
  _parser.Clear();
  _parser.Read(path, false);
  return _parser.GetWorldVolume();
}
//----------------------------------------------------------------------------------------------

//__Trivial DetectorConstruction for Geant4_____________________________________________________
class _empty_construction : public G4VUserDetectorConstruction {
  G4VPhysicalVolume* Construct() { return _world; }
};
//----------------------------------------------------------------------------------------------

//__Build Geometry Tree and Insertion-Order List________________________________________________
void _setup_geometry(const _geometric_volume& top) {
  const auto& volume = top.volume->GetLogicalVolume();
  const auto size = volume->GetNoDaughters();
  for (auto i = 0; i < size; ++i) {
    const auto& daughter = volume->GetDaughter(i);
    const auto& name = daughter->GetName();
    _geometric_volume daughter_geometry{
      daughter,
      G4AffineTransform().InverseProduct(
        top.transform,
        G4AffineTransform(daughter->GetFrameRotation(),
                          daughter->GetFrameTranslation()))
    };
    _geometry[name] = daughter_geometry;  // FIXME: overwrite if duplicate name! is this wanted?
    _geometry_insertion_order.push_back(name);
    _setup_geometry(daughter_geometry);
  }
}
//----------------------------------------------------------------------------------------------

//__Geometry Traversal Helper Function__________________________________________________________
template<class BinaryFunction>
BinaryFunction _traverse_geometry(BinaryFunction f) {
  for (const auto& name : _geometry_insertion_order) {
    const auto& search = _geometry.find(name);
    if (search != _geometry.cend())
      if (!f(search->first, search->second)) break;
  }
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Volume Containment Check____________________________________________________________________
bool _in_volume(const r3_point& point,
                const _geometric_volume& gvolume) {
  const auto& volume = gvolume.volume;
  const auto& transform = gvolume.transform;
  const auto& translation = transform.NetTranslation();
  const auto& rotation = transform.NetRotation();
  const auto& position = _to_local_transform(_to_G4ThreeVector(point), translation, rotation);
  // (rotation.inverse() * _to_G4ThreeVector(point)) - translation;
  return volume->GetLogicalVolume()->GetSolid()->Inside(position);
}
//----------------------------------------------------------------------------------------------

//__Find Volume Hierarchy_______________________________________________________________________
std::vector<_geometric_volume> _get_volume_hierarchy(const r3_point& point) {
  std::vector<_geometric_volume> out{};
  _traverse_geometry([&](const auto&, const auto& gvolume) {
    if (_in_volume(point, gvolume))
      out.push_back(gvolume);
    return true;
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Find Volume_________________________________________________________________________________
_geometric_volume _get_volume(const r3_point& point) {
  const auto& search = _name_cache.find(point);
  if (search != _name_cache.cend()) {
    return _geometry[search->second];
  } else {
    _geometric_volume out{};
    _traverse_geometry([&](const auto&, const auto& gvolume) {
      if (_in_volume(point, gvolume))
        out = gvolume;
      return true;
    });
    _name_cache[point] = _get_name(out);
    return out;
  }
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Initialize Geant4 Geometry Manager__________________________________________________________
void open(const std::string& path,
          const real default_time_error,
          const time_resolution_map& map) {
  close();
  _path = path;
  static const std::string _bar(61, '-');
  std::cout << "\nInitialize Geant4 Geometry Manager:\n" << _bar << "\n\n";
  _world = _load_geometry(path);
  _manager = new G4RunManager;
  _manager->SetVerboseLevel(0);
  _manager->SetUserInitialization(new _empty_construction);
  _manager->InitializeGeometry();
  _setup_geometry({_world, G4AffineTransform()});
  _default_time_resolution = default_time_error;
  _time_resolution_map = map;
  std::cout << _bar << "\n\n";
}
//----------------------------------------------------------------------------------------------

//__Close Geant4 Geometry Manager_______________________________________________________________
void close() {
  _path.clear();
  _geometry_insertion_order.clear();
  _geometry.clear();
  _time_resolution_map.clear();
  _name_cache.clear();
  // _box_cache.clear();
  delete _world;
  delete _manager;
}
//----------------------------------------------------------------------------------------------

//__Current Geometry File_______________________________________________________________________
const std::string& current_geometry_path() {
  return _path;
}
//----------------------------------------------------------------------------------------------

//__List of all Geometry Volumes________________________________________________________________
const structure_vector full_structure() {
  structure_vector out;
  out.reserve(_geometry.size());
  std::transform(_geometry.cbegin(), _geometry.cend(), std::back_inserter(out),
    [](const auto& element) { return element.first; });
  return out;
}
//----------------------------------------------------------------------------------------------

//__List of all Geometry Volumes Except those in the List_______________________________________
const structure_vector full_structure_except(const std::vector<std::string>& names) {
  structure_vector out;
  out.reserve(_geometry.size());
  const auto names_begin = names.cbegin();
  const auto names_end = names.cend();
  std::for_each(_geometry.cbegin(), _geometry.cend(),
    [&](const auto& element) {
      const auto& name = element.first;
      if (std::find(names_begin, names_end, name) == names_end)
        out.push_back(name);
    });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Get Default Time Resolution for Detector Volumes____________________________________________
real default_time_resolution() {
  return _default_time_resolution;
}
//----------------------------------------------------------------------------------------------

//__Volume Hierarchy of a Volume________________________________________________________________
const structure_vector volume_hierarchy(const structure_value& name) {
  // TODO: finish
  structure_vector out{};
  _traverse_geometry([&](const auto&, const auto&) {
    out.push_back(name);
    return true;
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Check Volume Containment____________________________________________________________________
bool is_inside_volume(const r3_point& point,
                      const structure_value& name) {
  const auto& search = _geometry.find(name);
  return search != _geometry.cend() && _in_volume(point, search->second);
}
bool is_inside_volume(const r4_point& point,
                      const structure_value& name) {
  return is_inside_volume(reduce_to_r3(point), name);
}
//----------------------------------------------------------------------------------------------

//__Get Volume Hierarchy________________________________________________________________________
const structure_vector volume_hierarchy(const r3_point& point) {
  const auto& hierarchy = _get_volume_hierarchy(point);
  structure_vector out;
  std::transform(hierarchy.cbegin(), hierarchy.cend(), std::back_inserter(out), _get_name);
  return out;
}
const structure_vector volume_hierarchy(const r4_point& point) {
  return volume_hierarchy(reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

//__Get Volume Name_____________________________________________________________________________
const structure_value volume(const r3_point& point) {
  return _get_name(_get_volume(point));
}
const structure_value volume(const r4_point& point) {
  return volume(reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

//__Box Volume Stream Overload__________________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const box_volume& box) {
  return os << "[" << box.min << "..." << box.center << "..." << box.max << "]";
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Get Solid from Geometry Volume______________________________________________________________
const G4VSolid* _get_solid(const _geometric_volume& gvolume) {
  const auto volume = gvolume.volume;
  return volume ? volume->GetLogicalVolume()->GetSolid() : nullptr;
}
//----------------------------------------------------------------------------------------------

//__Calculate Local Bounding Box________________________________________________________________
void _calculate_local_extent(const G4VSolid* solid,
                             G4ThreeVector& min_vector,
                             G4ThreeVector& max_vector) {
  G4double min, max;
  solid->CalculateExtent(EAxis::kXAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setX(min);
  max_vector.setX(max);
  solid->CalculateExtent(EAxis::kYAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setY(min);
  max_vector.setY(max);
  solid->CalculateExtent(EAxis::kZAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setZ(min);
  max_vector.setZ(max);
}
//----------------------------------------------------------------------------------------------

//__Transform Vector then Convert to R3 Point___________________________________________________
const r3_point _convert_transform(const G4ThreeVector& vector,
                                  const G4ThreeVector& translation,
                                  const G4RotationMatrix& rotation) {
  return _to_r3_point(_to_global_transform(vector, translation, rotation));
}
//----------------------------------------------------------------------------------------------

//__Calculate Global Bounding Box_______________________________________________________________
const box_volume _calculate_global_extent(const G4VSolid* solid,
                                          const G4ThreeVector& translation,
                                          const G4RotationMatrix& rotation) {
  G4ThreeVector min_vector, max_vector;
  _calculate_local_extent(solid, min_vector, max_vector);
  const auto min = _convert_transform(min_vector, translation, rotation);
  const auto max = _convert_transform(max_vector, translation, rotation);
  const auto minmax_x = std::minmax(min.x, max.x);
  const auto minmax_y = std::minmax(min.y, max.y);
  const auto minmax_z = std::minmax(min.z, max.z);
  return {
    _convert_transform(G4ThreeVector(), translation, rotation),
    {minmax_x.first, minmax_y.first, minmax_z.first},
    {minmax_x.second, minmax_y.second, minmax_z.second}
  };
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Check Volume Limits_________________________________________________________________________
const box_volume limits_of(const structure_value& name) {
  const auto& geometry_search = _geometry.find(name);
  if (geometry_search == _geometry.cend())
    return {};

  const auto& gvolume = geometry_search->second;
  const auto solid = _get_solid(gvolume);
  if (!solid)
    return {};

  const auto& transform = gvolume.transform;
  const auto out = _calculate_global_extent(solid, transform.NetTranslation(), transform.NetRotation());

  const auto& name_search = _name_cache.find(out.center);
  if (name_search == _name_cache.cend())
    _name_cache[out.center] = name;

  return out;
}
//----------------------------------------------------------------------------------------------

//__Time Error on a Detector Component__________________________________________________________
real time_resolution_of(const structure_value& name) {
  const auto search = _time_resolution_map.find(name);
  return search != _time_resolution_map.cend() ? search->second : _default_time_resolution;
}
//----------------------------------------------------------------------------------------------

//__Compute Intersection of Box Volumes_________________________________________________________
const box_volume coordinatewise_intersection(const box_volume& first,
                                             const box_volume& second) {
  box_volume out{};

  if (!(first.max.x < second.min.x || second.max.x < first.min.x)) {
    out.min.x = std::max(first.min.x, second.min.x);
    out.max.x = std::min(first.max.x, second.max.x);
  }
  if (!(first.max.y < second.min.y || second.max.y < first.min.y)) {
    out.min.y = std::max(first.min.y, second.min.y);
    out.max.y = std::min(first.max.y, second.max.y);
  }
  if (!(first.max.z < second.min.z || second.max.z < first.min.z)) {
    out.min.z = std::max(first.min.z, second.min.z);
    out.max.z = std::min(first.max.z, second.max.z);
  }
  out.center = 0.5L * (out.min + out.max);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Compute Union of Box Volumes________________________________________________________________
const box_volume coordinatewise_union(const box_volume& first,
                                      const box_volume& second) {
  box_volume out{};
  out.min.x = std::min(first.min.x, second.min.x);
  out.max.x = std::max(first.max.x, second.max.x);
  out.min.y = std::min(first.min.y, second.min.y);
  out.max.y = std::max(first.max.y, second.max.y);
  out.min.z = std::min(first.min.z, second.min.z);
  out.max.z = std::max(first.max.z, second.max.z);
  out.center = 0.5L * (out.min + out.max);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Limit Box of the Volume with Point__________________________________________________________
const box_volume limits_of_volume(const r3_point& point) {
  return limits_of(volume(point));
}
const box_volume limits_of_volume(const r4_point& point) {
  return limits_of(volume(point));
}
//----------------------------------------------------------------------------------------------

//__Time Error on a Detector Component from a Point_____________________________________________
real time_resolution_of_volume(const r3_point& point) {
  return time_resolution_of(volume(point));
}
real time_resolution_of_volume(const r4_point& point) {
  return time_resolution_of_volume(reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

//__Box Volume Containment Check________________________________________________________________
constexpr bool is_inside_volume(const r3_point& point,
                                const box_volume& box) {
  return util::algorithm::between(point.x, box.min.x, box.max.x)
      && util::algorithm::between(point.y, box.min.y, box.max.y)
      && util::algorithm::between(point.z, box.min.z, box.max.z);
}
constexpr bool is_inside_volume(const r4_point& point,
                                const box_volume& box) {
  return is_inside_volume(reduce_to_r3(point), box);
}
//----------------------------------------------------------------------------------------------

//__Find Center of Geometry around Point________________________________________________________
const r3_point find_center(const structure_value& name) {
  return limits_of(name).center;
}
const r3_point find_center(const r3_point& point) {
  return limits_of_volume(point).center;
}
const r4_point find_center(const r4_point& point) {
  const auto& center = limits_of_volume(point).center;
  return {point.t, center.x, center.y, center.z};
}
//----------------------------------------------------------------------------------------------

} /* namespace geometry */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

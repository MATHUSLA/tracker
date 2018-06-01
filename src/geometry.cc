/*
 * src/geometry.cc
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

#include <unordered_map>

#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VUserDetectorConstruction.hh>
#include <Geant4/G4RunManager.hh>
#include <Geant4/G4VoxelLimits.hh>
#include <Geant4/G4AffineTransform.hh>

#include "util/algorithm.hh"

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
thread_local std::unordered_map<std::string, _geometric_volume> _geometry;
thread_local std::vector<std::string> _geometry_insertion_order;
thread_local std::unordered_map<r3_point, std::string, r3_point_hash> _name_cache;
//thread_local std::unordered_map<r3_point, box_volume, r3_point_hash> _box_cache;
//----------------------------------------------------------------------------------------------

//__Convert From G4ThreeVector to R3_Point______________________________________________________
r3_point _to_r3_point(const G4ThreeVector& vector) {
  return { vector.x(), vector.y(), vector.z() };
}
//----------------------------------------------------------------------------------------------

//__Convert From R3_Point to G4ThreeVector______________________________________________________
G4ThreeVector _to_G4ThreeVector(const r3_point& point) {
  return G4ThreeVector(point.x, point.y, point.z);
}
//----------------------------------------------------------------------------------------------

//__Safely Get Name of Volume___________________________________________________________________
const std::string _safe_get_name(const _geometric_volume& gvolume) {
  const auto& volume = gvolume.volume;
  return volume ? volume->GetName() : "";
}
//----------------------------------------------------------------------------------------------

//__Unsafely Get Name of Volume_________________________________________________________________
const std::string _unsafe_get_name(const _geometric_volume& gvolume) {
  return gvolume.volume->GetName();
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
    if (search != _geometry.end())
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
  const auto& position = (rotation.inverse() * _to_G4ThreeVector(point)) - translation;
  return volume->GetLogicalVolume()->GetSolid()->Inside(position);
}
//----------------------------------------------------------------------------------------------

//__Find Volume Hierarchy_______________________________________________________________________
std::vector<_geometric_volume> _get_volume_hierarchy(const r3_point& point) {
  std::vector<_geometric_volume> out{};
  _traverse_geometry([&](const auto& _, const auto& gvolume) {
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
  if (search != _name_cache.end()) {
    return _geometry[search->second];
  } else {
    _geometric_volume out{};
    _traverse_geometry([&](const auto& _, const auto& gvolume) {
      if (_in_volume(point, gvolume))
        out = gvolume;
      return true;
    });
    _name_cache[point] = _unsafe_get_name(out);
    return out;
  }
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Initialize Geant4 Geometry Manager__________________________________________________________
void open(const std::string& path) {
  close();
  _path = path;
  const static std::string&& _bar = std::string(61, '-');
  std::cout << "\nInitialize Geant4 Geometry Manager:\n" << _bar << "\n\n";
  _world = _load_geometry(path);
  _manager = new G4RunManager;
  _manager->SetVerboseLevel(0);
  _manager->SetUserInitialization(new _empty_construction);
  _manager->InitializeGeometry();
  _setup_geometry({_world, G4AffineTransform()});
  std::cout << _bar << "\n\n";
}
//----------------------------------------------------------------------------------------------

//__Close Geant4 Geometry Manager_______________________________________________________________
void close() {
  _path.clear();
  _geometry_insertion_order.clear();
  _geometry.clear();
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
const std::vector<std::string> full_structure() {
  std::vector<std::string> out;
  out.reserve(_geometry.size());
  std::transform(_geometry.cbegin(), _geometry.cend(), std::back_inserter(out),
    [](const auto& element) { return element.first; });
  return out;
}
//----------------------------------------------------------------------------------------------

//__List of all Geometry Volumes Except those in the List_______________________________________
const std::vector<std::string> full_structure_except(const std::vector<std::string>& names) {
  std::vector<std::string> out;
  out.reserve(_geometry.size());
  const auto names_begin = names.cbegin();
  const auto names_end = names.cend();
  std::for_each(_geometry.cbegin(), _geometry.cend(),
    [&](const auto& element) {
      const auto& name = element.first;
      if (std::find(names_begin, names_end, name) == names_end) out.push_back(name);
    });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Volume Hierarchy of a Volume________________________________________________________________
const std::vector<std::string> volume_hierarchy(const std::string& name) {
  std::vector<std::string> out{};
  _traverse_geometry([&](const auto& _, const auto& __) {
    out.push_back(name);
    return true;
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Check Volume Containment____________________________________________________________________
bool is_inside_volume(const r3_point& point,
                      const std::string& name) {
  const auto& search = _geometry.find(name);
  return search != _geometry.end() && _in_volume(point, search->second);
}
bool is_inside_volume(const r4_point& point,
                      const std::string& name) {
  return is_inside_volume(reduce_to_r3(point), name);
}
//----------------------------------------------------------------------------------------------

//__Get Volume Hierarchy________________________________________________________________________
const std::vector<std::string> volume_hierarchy(const r3_point& point) {
  const auto& hierarchy = _get_volume_hierarchy(point);
  std::vector<std::string> out;
  std::transform(hierarchy.cbegin(), hierarchy.cend(), std::back_inserter(out), _safe_get_name);
  return out;
}
const std::vector<std::string> volume_hierarchy(const r4_point& point) {
  return volume_hierarchy(reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

//__Get Volume Name_____________________________________________________________________________
const std::string volume(const r3_point& point) {
  return _safe_get_name(_get_volume(point));
}
const std::string volume(const r4_point& point) {
  return volume(reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

//__Check Volume Limits_________________________________________________________________________
const box_volume limits_of(const std::string& name) {
  static const G4VoxelLimits _blank_voxels;
  static const G4AffineTransform _blank_transform;

  box_volume out{};

  const auto& geometry_search = _geometry.find(name);
  if (geometry_search == _geometry.end())
    return out;

  const auto& gvolume = geometry_search->second;
  const auto& volume = gvolume.volume;
  if (!volume)
    return out;

  const auto& transform = gvolume.transform;
  const auto& rotation = transform.NetRotation();
  const auto& translation = transform.NetTranslation();

  const auto& center = rotation * translation;
  out.center = _to_r3_point(center);

  const auto& name_search = _name_cache.find(out.center);
  if (name_search == _name_cache.end())
    _name_cache[out.center] = name;

  G4double min, max;
  G4ThreeVector min_vector, max_vector;

  const auto& solid = volume->GetLogicalVolume()->GetSolid();
  solid->CalculateExtent(EAxis::kXAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setX(min);
  max_vector.setX(max);
  solid->CalculateExtent(EAxis::kYAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setY(min);
  max_vector.setY(max);
  solid->CalculateExtent(EAxis::kZAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setZ(min);
  max_vector.setZ(max);

  out.min = _to_r3_point(rotation * (min_vector + translation));
  out.max = _to_r3_point(rotation * (max_vector + translation));

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

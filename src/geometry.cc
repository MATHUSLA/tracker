#include "geometry.hh"

#include <unordered_map>

#include "Geant4/G4GDMLParser.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4VUserDetectorConstruction.hh"
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4VoxelLimits.hh"
#include "Geant4/G4AffineTransform.hh"

namespace MATHUSLA { namespace TRACKER {

namespace geometry { ///////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

struct _geometric_volume { G4VPhysicalVolume* volume; G4AffineTransform transform; };

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

//__Static Variables____________________________________________________________________________
thread_local G4RunManager* _manager;
thread_local G4VPhysicalVolume* _world;
thread_local std::unordered_map<std::string, _geometric_volume> _geometry;
thread_local std::vector<std::string> _geometry_insertion_order;
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
  const auto&& size = volume->GetNoDaughters();
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
      f(search->first, search->second);
  }
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Volume Containment Check____________________________________________________________________
bool _in_volume(const r3_point& point, const _geometric_volume& gvolume) {
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
  _traverse_geometry([&](const auto& name, const auto& gvolume) {
    if (_in_volume(point, gvolume)) out.push_back(gvolume);
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Find Volume_________________________________________________________________________________
_geometric_volume _get_volume(const r3_point& point) {
  _geometric_volume out{};
  _traverse_geometry([&](const auto& name, const auto& gvolume) {
    if (_in_volume(point, gvolume)) out = gvolume;
  });
  return out;
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Initialize Geant4 Geometry Manager__________________________________________________________
void open(const std::string& path) {
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
  delete _manager;
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
  std::vector<std::string> out{};
  std::transform(hierarchy.cbegin(), hierarchy.cend(), std::back_inserter(out),
    [](const auto& gvolume) { return gvolume.volume->GetName(); });
  return out;
}
const std::vector<std::string> volume_hierarchy(const r4_point& point) {
  return volume_hierarchy(reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

//__Get Volume Name_____________________________________________________________________________
const std::string volume(const r3_point& point) {
  const auto& gvolume = _get_volume(point);
  const auto& volume = gvolume.volume;
  return volume ? volume->GetName() : "";
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

  const auto& search = _geometry.find(name);
  if (search == _geometry.end())
    return out;

  const auto& gvolume = search->second;
  const auto& volume = gvolume.volume;
  if (!volume)
    return out;

  const auto& transform = gvolume.transform;
  const auto& rotation = transform.NetRotation();
  const auto& translation = transform.NetTranslation();

  const auto& center = rotation * translation;
  out.center = _to_r3_point(center);

  G4double min, max;
  G4ThreeVector min_vector, max_vector;

  const auto& solid = volume->GetLogicalVolume()->GetSolid();
  solid->CalculateExtent(kXAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setX(min);
  max_vector.setX(max);
  solid->CalculateExtent(kYAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setY(min);
  max_vector.setY(max);
  solid->CalculateExtent(kZAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setZ(min);
  max_vector.setZ(max);

  out.min = _to_r3_point(rotation * (min_vector + translation));
  out.max = _to_r3_point(rotation * (max_vector + translation));

  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace geometry */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

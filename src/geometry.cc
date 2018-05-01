#include "geometry.hh"

#include <unordered_map>

#include "Geant4/G4GDMLParser.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4VTouchable.hh"
#include "Geant4/G4VUserDetectorConstruction.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/G4Navigator.hh"
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4VisExecutive.hh"
#include "Geant4/G4VVisManager.hh"
#include "Geant4/G4VMarker.hh"
#include "Geant4/G4Circle.hh"
#include "Geant4/G4VisAttributes.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4UImanager.hh"
#include "Geant4/G4UIExecutive.hh"
#include "Geant4/G4PhysicalVolumeStore.hh"
#include "Geant4/G4VoxelLimits.hh"
#include "Geant4/G4AffineTransform.hh"

namespace MATHUSLA { namespace TRACKER {

namespace { ////////////////////////////////////////////////////////////////////////////////////

struct _geometric_volume { G4VPhysicalVolume* volume; G4AffineTransform transform; };

//__Geant4 Static Variables_____________________________________________________________________
static G4ThreadLocal G4GDMLParser        _parser;
static G4ThreadLocal std::string         _path;
static G4ThreadLocal G4RunManager*       _manager;
static G4ThreadLocal G4VPhysicalVolume*  _world;
static G4ThreadLocal G4Navigator         _navigator;
static G4ThreadLocal G4NavigationHistory _history;
static G4ThreadLocal std::unordered_map<std::string, _geometric_volume> _geometry_tree;
//----------------------------------------------------------------------------------------------

//__Convert Between R3_Point and G4ThreeVector__________________________________________________
type::r3_point _convert(const G4ThreeVector& vector) {
  return { vector.x(), vector.y(), vector.z() };
}
//----------------------------------------------------------------------------------------------

//__Load Geometry from GDML File________________________________________________________________
static G4VPhysicalVolume* _load_geometry(const std::string& path) {
  _parser.Clear();
  _parser.Read(path, false);
  return _parser.GetWorldVolume();
}
//----------------------------------------------------------------------------------------------

static void _setup_geometry(const _geometric_volume& top) {
  const auto& volume = top.volume->GetLogicalVolume();
  const auto&& size = volume->GetNoDaughters();
  for (auto i = 0; i < size; ++i) {
    const auto& daughter = volume->GetDaughter(i);
    _geometric_volume daughter_geometry{
      daughter,
      G4AffineTransform().InverseProduct(
        top.transform,
        G4AffineTransform(daughter->GetFrameRotation(),
                          daughter->GetFrameTranslation()))
    };
    _geometry_tree[daughter->GetName()] = daughter_geometry;  // FIXME: overwrite if duplicate name!
    _setup_geometry(daughter_geometry);
  }
}

static G4VPhysicalVolume* _get_volume(const std::string& name) {
  const auto& search = _geometry_tree.find(name);
  return search != _geometry_tree.end() ? search->second.volume : nullptr;
}

static const G4AffineTransform _get_transform(G4VPhysicalVolume* volume) {
  const auto& search = _geometry_tree.find(volume->GetName());
  return search != _geometry_tree.end() ? search->second.transform : G4AffineTransform();
}

//__Trivial DetectorConstruction for Geant4_____________________________________________________
class _empty_construction : public G4VUserDetectorConstruction {
  G4VPhysicalVolume* Construct() {
    _world = _load_geometry(_path);
    _navigator.SetWorldVolume(_world);
    return _world;
  }
};
//----------------------------------------------------------------------------------------------

//__Get Volume from World_______________________________________________________________________
static G4VPhysicalVolume* _get_volume(const long double x,
                                      const long double y,
                                      const long double z) {
  static G4ThreadLocal G4ThreeVector _vector;
  _vector.set(x * Units::Length, y * Units::Length, z * Units::Length);
  return _navigator.LocateGlobalPointAndSetup(_vector, nullptr, false, false);
}
//----------------------------------------------------------------------------------------------

//__Get Touchable from World____________________________________________________________________
static G4VTouchable* _get_touchable(const long double x,
                                    const long double y,
                                    const long double z) {
  static G4ThreadLocal G4ThreeVector _vector;
  _vector.set(x * Units::Length, y * Units::Length, z * Units::Length);
  auto out = new G4TouchableHistory();
  _navigator.LocateGlobalPointAndUpdateTouchable(_vector, G4ThreeVector(), out);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Get Name of Physical Volume_________________________________________________________________
static const std::string _volume_name(const long double x,
                                      const long double y,
                                      const long double z) {
  const auto volume = _get_volume(x, y, z);
  return volume ? volume->GetName() : "";
}
//----------------------------------------------------------------------------------------------

//__Get Hierarchy of Volumes____________________________________________________________________
static const std::vector<std::string> _volume_hierarchy(const long double x,
                                                        const long double y,
                                                        const long double z) {
  std::vector<std::string> out{};
  const auto touchable = _get_touchable(x, y, z);
  const auto depth = touchable->GetHistoryDepth();
  for (auto i = 0; i < depth; ++i)
    out.push_back(touchable->GetVolume(i)->GetName());
  return out;
}
//----------------------------------------------------------------------------------------------

//__Check Volume Containment____________________________________________________________________
static inline bool _is_inside_volume(const long double x,
                                     const long double y,
                                     const long double z,
                                     const std::string& name) {
  const auto touchable = _get_touchable(x, y, z);
  const auto depth = touchable->GetHistoryDepth();
  for (auto i = 0; i < depth; ++i)
    if (name == touchable->GetVolume(i)->GetName())
      return true;
  return false;
}
//----------------------------------------------------------------------------------------------

static inline bool _in_volume(const long double x,
                              const long double y,
                              const long double z,
                              const std::string& name) {
  const auto& search = _geometry_tree.find(name);
  if (search != _geometry_tree.end()) {
    const auto& gvolume = search->second;
    const auto& volume = gvolume.volume;
    const auto& transform = gvolume.transform;
    const auto& translation = transform.NetTranslation();
    const auto& rotation = transform.NetRotation();
    const auto& position = (rotation.inverse() * G4ThreeVector(x, y, z)) - translation;

    std::cout << position << "\n";

    return volume->GetLogicalVolume()->GetSolid()->Inside(position);
  }
  return false;
}

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

namespace geometry { ///////////////////////////////////////////////////////////////////////////

//__Initialize Geant4 Geometry Manager__________________________________________________________
void open(const std::string& path) {
  const static std::string&& _bar = std::string(61, '-');
  std::cout << "\nInitialize Geant4 Geometry Manager:\n\n" << _bar;
  _path = path;
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

const detector_limits limits_of(const std::string name) {
  static G4VoxelLimits _blank_voxels;
  static G4AffineTransform _blank_transform;

  detector_limits out{};
  const auto detector = _get_volume(name);
  if (!detector)
    return out;

  const auto transform = _get_transform(detector);
  const auto rotation = transform.NetRotation();
  const auto translation = transform.NetTranslation();
  std::cout << transform << "\n";

  const auto center = rotation * translation;
  out.center = _convert(center);

  auto min_vector = G4ThreeVector(0, 0, 0);
  auto max_vector = G4ThreeVector(0, 0, 0);

  real min, max;
  const auto solid = detector->GetLogicalVolume()->GetSolid();
  solid->CalculateExtent(kXAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setX(min);
  max_vector.setX(max);
  solid->CalculateExtent(kYAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setY(min);
  max_vector.setY(max);
  solid->CalculateExtent(kZAxis, _blank_voxels, _blank_transform, min, max);
  min_vector.setZ(min);
  max_vector.setZ(max);

  out.min = _convert(rotation * (min_vector + translation));
  out.max = _convert(rotation * (max_vector + translation));

  return out;
}

//__Check Volume Containment____________________________________________________________________
bool is_inside_volume(const r3_point& point,
                      const std::string& name) {
  return _in_volume(point.x, point.y, point.z, name);
}
bool is_inside_volume(const r4_point& point,
                      const std::string& name) {
  return is_inside_volume(reduce_to_r3(point), name);
}
//----------------------------------------------------------------------------------------------

//__Get Volume Name_____________________________________________________________________________
const std::string volume(const r3_point& point) {
  return _volume_name(point.x, point.y, point.z);
}
const std::string volume(const r4_point& point) {
  return volume(reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

//__Get Volume Hierarchy________________________________________________________________________
const std::vector<std::string> volume_hierarchy(const r3_point& point) {
  return _volume_hierarchy(point.x, point.y, point.z);
}
const std::vector<std::string> volume_hierarchy(const r4_point& point) {
  return volume_hierarchy(reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

} /* namespace geometry */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

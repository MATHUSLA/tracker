#include "geometry.hh"

#include "Geant4/G4GDMLParser.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4VTouchable.hh"
#include "Geant4/G4VUserDetectorConstruction.hh"
#include "Geant4/G4Navigator.hh"
#include "Geant4/G4RunManager.hh"

namespace MATHUSLA { namespace TRACKER {

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Geant4 Static Variables_____________________________________________________________________
static G4ThreadLocal G4GDMLParser       _parser;
static G4ThreadLocal std::string        _path;
static G4ThreadLocal G4VPhysicalVolume* _world;
static G4ThreadLocal G4Navigator        _navigator;
//----------------------------------------------------------------------------------------------

//__Load Geometry from GDML File________________________________________________________________
static G4VPhysicalVolume* _load_geometry(const std::string& path) {
  _parser.Clear();
  _parser.Read(path, false);
  return _parser.GetWorldVolume();
}
//----------------------------------------------------------------------------------------------

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
static G4VPhysicalVolume* _get_volume(const double x,
                                      const double y,
                                      const double z) {
  static G4ThreadLocal G4ThreeVector _vector;
  _vector.set(x * Units::Length, y * Units::Length, z * Units::Length);
  return _navigator.LocateGlobalPointAndSetup(_vector);
}
//----------------------------------------------------------------------------------------------

//__Get Touchable from World____________________________________________________________________
static G4VTouchable* _get_touchable(const double x,
                                    const double y,
                                    const double z) {
  static G4ThreadLocal G4ThreeVector _vector;
  _vector.set(x * Units::Length, y * Units::Length, z * Units::Length);
  auto out = new G4TouchableHistory();
  _navigator.LocateGlobalPointAndUpdateTouchable(_vector, G4ThreeVector(), out);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Get Name of Physical Volume_________________________________________________________________
static const std::string _volume_name(const double x,
                                      const double y,
                                      const double z) {
  const auto volume = _get_volume(x, y, z);
  return volume ? volume->GetName() : "";
}
//----------------------------------------------------------------------------------------------

//__Get Hierarchy of Volumes____________________________________________________________________
static const std::vector<std::string> _volume_hierarchy(const double x,
                                                        const double y,
                                                        const double z) {
  std::vector<std::string> out{};
  const auto touchable = _get_touchable(x, y, z);
  const auto depth = touchable->GetHistoryDepth();
  for (auto i = 0; i < depth; ++i)
    out.push_back(touchable->GetVolume(i)->GetName());
  return out;
}
//----------------------------------------------------------------------------------------------

//__Check Volume Containment____________________________________________________________________
static inline bool _within_volume(const double x,
                                  const double y,
                                  const double z,
                                  const std::string& name) {
  const auto touchable = _get_touchable(x, y, z);
  const auto depth = touchable->GetHistoryDepth();
  for (auto i = 0; i < depth; ++i)
    if (name == touchable->GetVolume(i)->GetName())
      return true;
  return false;
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Initialize Geant4 Geometry Manager__________________________________________________________
void Geometry::Initialize(const std::string& path) {
  const static std::string&& _bar = std::string(61, '-');
  std::cout << "\nInitialize Geant4 Geometry Manager:\n\n" << _bar;
  _path = path;
  auto manager = new G4RunManager;
  manager->SetVerboseLevel(0);
  manager->SetUserInitialization(new _empty_construction);
  manager->InitializeGeometry();
  std::cout << _bar << "\n\n";
}
//----------------------------------------------------------------------------------------------

//__Get Volume Name_____________________________________________________________________________
const std::string Geometry::Volume(const r3_point& point) {
  return _volume_name(point.x, point.y, point.z);
}
const std::string Geometry::Volume(const r4_point& point) {
  return _volume_name(point.x, point.y, point.z);
}
//----------------------------------------------------------------------------------------------

//__Get Volume Hierarchy________________________________________________________________________
const std::vector<std::string> Geometry::VolumeHierarchy(const r3_point& point) {
  return _volume_hierarchy(point.x, point.y, point.z);
}
const std::vector<std::string> Geometry::VolumeHierarchy(const r4_point& point) {
  return _volume_hierarchy(point.x, point.y, point.z);
}
//----------------------------------------------------------------------------------------------

//__Check Volume Containment____________________________________________________________________
bool Geometry::WithinVolume(const r3_point& point,
                            const std::string& name) {
  return _within_volume(point.x, point.y, point.z, name);
}
bool Geometry::WithinVolume(const r4_point& point,
                            const std::string& name) {
  return _within_volume(point.x, point.y, point.z, name);
}
//----------------------------------------------------------------------------------------------

} } /* namespace MATHUSLA::TRACKER */

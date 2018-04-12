#include "geometry.hh"

#include "Geant4/G4GeometryManager.hh"
#include "Geant4/G4PhysicalVolumeStore.hh"
#include "Geant4/G4LogicalVolumeStore.hh"
#include "Geant4/G4SolidStore.hh"
#include "Geant4/G4GDMLParser.hh"
#include "Geant4/G4Navigator.hh"

namespace MATHUSLA { namespace TRACKER {

Geometry::Geometry(const std::string& path)
    : G4VUserDetectorConstruction(), _geometry(path) {}

namespace {

static G4ThreadLocal G4GDMLParser _parser;
static G4ThreadLocal G4VPhysicalVolume* _world;
static G4ThreadLocal G4Navigator _navigator;

static inline G4VPhysicalVolume* _load_geometry(const std::string& path) {
  _parser.Clear();
  _parser.Read(path, false);
  return _parser.GetWorldVolume();
}

} /* anonymous namespace */

G4VPhysicalVolume* Geometry::Construct() {
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  _world = _load_geometry(_geometry);
  _navigator.SetWorldVolume(_world);
  return _world;
}

const std::string Geometry::DetectorName(const V3& point) {
  static G4ThreadLocal G4ThreeVector _vector;
  _vector.set(point.x, point.y, point.z);
  auto volume = _navigator.LocateGlobalPointAndSetup(_vector);
  return volume ? volume->GetName() : "";
}

} } /* namespace MATHUSLA::TRACKER */

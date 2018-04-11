#include "geometry/Construction.hh"

#include "Geant4/G4GeometryManager.hh"
#include "Geant4/G4PhysicalVolumeStore.hh"
#include "Geant4/G4LogicalVolumeStore.hh"
#include "Geant4/G4SolidStore.hh"

#include "geometry/Reader.hh"

namespace MATHUSLA { namespace TRACKER {

Construction::Construction(const std::string& geometry)
    : G4VUserDetectorConstruction(), _geometry(geometry) {}

G4VPhysicalVolume* Construction::Construct() {
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  return Reader::Load(_geometry);
}

} } /* namespace MATHUSLA::TRACKER */

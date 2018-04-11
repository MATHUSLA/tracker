#ifndef GEOMETRY_CONSTRUCTION_HH
#define GEOMETRY_CONSTRUCTION_HH
#pragma once

#include "Geant4/G4VUserDetectorConstruction.hh"
#include "Geant4/G4VPhysicalVolume.hh"

namespace MATHUSLA { namespace TRACKER {

class Construction : public G4VUserDetectorConstruction {
public:
  Construction(const std::string& _geometry);
  G4VPhysicalVolume* Construct();

private:
  std::string _geometry;
};

} } /* namespace MATHUSLA::TRACKER */

#endif /* GEOMETRY_CONSTRUCTION_HH */

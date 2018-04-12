#ifndef TRACKER__GEOMETRY_HH
#define TRACKER__GEOMETRY_HH
#pragma once

#include "Geant4/G4VUserDetectorConstruction.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4SystemOfUnits.hh"

namespace MATHUSLA { namespace TRACKER {

class Geometry : public G4VUserDetectorConstruction {
public:
  Geometry(const std::string& path);
  G4VPhysicalVolume* Construct();

  struct V3 { double    x, y, z; };
  struct V4 { double t, x, y, z; };

  static const std::string DetectorName(const V3& point);

private:
  std::string _geometry;
};

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__GEOMETRY_HH */

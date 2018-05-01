#ifndef TRACKER__UNITS_HH
#define TRACKER__UNITS_HH
#pragma once

#include "Geant4/G4UnitsTable.hh"
#include "Geant4/G4SystemOfUnits.hh"
#include "Geant4/G4PhysicalConstants.hh"

namespace MATHUSLA {

//__Extention Momentum Units____________________________________________________________________
constexpr auto GeVperC = GeV;
constexpr auto MeVperC = MeV;
constexpr auto keVperC = keV;
constexpr auto  eVperC =  eV;
//----------------------------------------------------------------------------------------------

namespace Units { //////////////////////////////////////////////////////////////////////////////

//__Standard Units______________________________________________________________________________
constexpr auto Length   = cm;
constexpr auto Time     = ns;
constexpr auto Energy   = MeV;
constexpr auto Momentum = MeVperC;
//----------------------------------------------------------------------------------------------

//__Install Momentum Units into Geant4__________________________________________________________
inline void Define() {
  new G4UnitDefinition("GeV/c", "GeV/c", "Momentum", GeVperC);
  new G4UnitDefinition("MeV/c", "MeV/c", "Momentum", MeVperC);
  new G4UnitDefinition("keV/c", "keV/c", "Momentum", keVperC);
  new G4UnitDefinition( "eV/c",  "eV/c", "Momentum",  eVperC);
}
//----------------------------------------------------------------------------------------------

} /* namespace Units */ ////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* TRACKER__UNITS_HH */

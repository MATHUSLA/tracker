/*
 * include/units.hh
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

namespace units { //////////////////////////////////////////////////////////////////////////////

//__Standard Units______________________________________________________________________________
constexpr auto length   = cm;
constexpr auto time     = ns;
constexpr auto energy   = MeV;
constexpr auto momentum = MeVperC;
constexpr auto velocity = length / time;
//----------------------------------------------------------------------------------------------

//__Install Momentum Units into Geant4__________________________________________________________
inline void define() {
  new G4UnitDefinition("GeV/c", "GeV/c", "Momentum", GeVperC);
  new G4UnitDefinition("MeV/c", "MeV/c", "Momentum", MeVperC);
  new G4UnitDefinition("keV/c", "keV/c", "Momentum", keVperC);
  new G4UnitDefinition( "eV/c",  "eV/c", "Momentum",  eVperC);
}
//----------------------------------------------------------------------------------------------

//__Alias Units and Constants___________________________________________________________________
constexpr auto speed_of_light = c_light;
//----------------------------------------------------------------------------------------------

} /* namespace units */ ////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* TRACKER__UNITS_HH */

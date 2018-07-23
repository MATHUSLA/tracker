/*
 * include/tracker/core/units.hh
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

#ifndef TRACKER__CORE__UNITS_HH
#define TRACKER__CORE__UNITS_HH
#pragma once

#include <tracker/core/type.hh>

namespace MATHUSLA {

namespace units { //////////////////////////////////////////////////////////////////////////////

//__Shielded Include of Units and Constants_____________________________________________________
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4PhysicalConstants.hh>
//----------------------------------------------------------------------------------------------

//__Extention Momentum Units____________________________________________________________________
constexpr auto GeVperC = GeV;
constexpr auto MeVperC = MeV;
constexpr auto keVperC = keV;
constexpr auto  eVperC =  eV;
//----------------------------------------------------------------------------------------------

//__Standard Units______________________________________________________________________________
constexpr auto length   = cm;
constexpr auto time     = ns;
constexpr auto energy   = MeV;
constexpr auto momentum = MeVperC;
constexpr auto velocity = length / time;
constexpr auto angle    = rad;
//----------------------------------------------------------------------------------------------

//__Standard Unit String________________________________________________________________________
static const std::string& length_string   = "cm";
static const std::string& time_string     = "ns";
static const std::string& energy_string   = "MeV";
static const std::string& momentum_string = "MeVperC";
static const std::string& velocity_string = "cm/ns";
static const std::string& angle_string    = "rad";
//----------------------------------------------------------------------------------------------

//__Alias Units and Constants___________________________________________________________________
constexpr auto speed_of_light = c_light;
//----------------------------------------------------------------------------------------------

//__Scale R3 by Units___________________________________________________________________________
inline type::r3_point& scale_r3_length(type::r3_point& point) {
  return point /= length;
}
//----------------------------------------------------------------------------------------------

//__Scale R3 by Units___________________________________________________________________________
inline const type::r3_point scale_r3_length(const type::r3_point& point) {
  return point / length;
}
//----------------------------------------------------------------------------------------------

//__Scale R4 by Units___________________________________________________________________________
inline type::r4_point& scale_r4_length(type::r4_point& point) {
  point.t /= time;
  point.x /= length;
  point.y /= length;
  point.z /= length;
  return point;
}
//----------------------------------------------------------------------------------------------

//__Scale R4 by Units___________________________________________________________________________
inline const type::r4_point scale_r4_length(const type::r4_point& point) {
  return type::r4_point{point.t / time, point.x / length, point.y / length, point.z / length};
}
//----------------------------------------------------------------------------------------------

//__Scale R3 by Units___________________________________________________________________________
inline type::r3_point& scale_r3_momentum(type::r3_point& point) {
  return point /= momentum;
}
//----------------------------------------------------------------------------------------------

//__Scale R3 by Units___________________________________________________________________________
inline const type::r3_point scale_r3_momentum(const type::r3_point& point) {
  return point / momentum;
}
//----------------------------------------------------------------------------------------------

//__Scale R4 by Units___________________________________________________________________________
inline type::r4_point& scale_r4_momentum(type::r4_point& point) {
  point.t /= energy;
  point.x /= momentum;
  point.y /= momentum;
  point.z /= momentum;
  return point;
}
//----------------------------------------------------------------------------------------------

//__Scale R4 by Units___________________________________________________________________________
inline const type::r4_point scale_r4_momentum(const type::r4_point& point) {
  return type::r4_point{point.t / energy, point.x / momentum, point.y / momentum, point.z / momentum};
}
//----------------------------------------------------------------------------------------------

} /* namespace units */ ////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* TRACKER__CORE__UNITS_HH */

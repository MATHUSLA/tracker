#ifndef TRACKER__GEOMETRY_HH
#define TRACKER__GEOMETRY_HH
#pragma once

#include "types.hh"
#include "units.hh"

namespace MATHUSLA { namespace TRACKER {

namespace Geometry { ///////////////////////////////////////////////////////////////////////////

//__Initialize Geometry Navigation System_______________________________________________________
void Initialize(const std::string& path);
//----------------------------------------------------------------------------------------------

//__Volume Search_______________________________________________________________________________
const std::string Volume(const r3_point& point);
const std::string Volume(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Volume Hierarchy Search_____________________________________________________________________
const std::vector<std::string> VolumeHierarchy(const r3_point& point);
const std::vector<std::string> VolumeHierarchy(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Within Volume Check_________________________________________________________________________
bool WithinVolume(const r3_point& point, const std::string& name);
bool WithinVolume(const r4_point& point, const std::string& name);
//----------------------------------------------------------------------------------------------

} /* namespace Geometry */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__GEOMETRY_HH */

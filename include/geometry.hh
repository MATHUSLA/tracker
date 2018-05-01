#ifndef TRACKER__GEOMETRY_HH
#define TRACKER__GEOMETRY_HH
#pragma once

#include "point.hh"
#include "units.hh"

namespace MATHUSLA { namespace TRACKER {

namespace geometry { ///////////////////////////////////////////////////////////////////////////

using namespace type;

//__Geometry Navigation System__________________________________________________________________
void open(const std::string& path);
void close();
//----------------------------------------------------------------------------------------------

struct detector_limits { r3_point center, min, max; };
const detector_limits limits_of(const std::string name);

//__Volume Containment Check____________________________________________________________________
bool is_inside_volume(const r3_point& point, const std::string& name);
bool is_inside_volume(const r4_point& point, const std::string& name);
//----------------------------------------------------------------------------------------------

//__Volume Search_______________________________________________________________________________
const std::string volume(const r3_point& point);
const std::string volume(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Volume Hierarchy Search_____________________________________________________________________
const std::vector<std::string> volume_hierarchy(const r3_point& point);
const std::vector<std::string> volume_hierarchy(const r4_point& point);
//----------------------------------------------------------------------------------------------

} /* namespace geometry */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__GEOMETRY_HH */

#ifndef TRACKER__READER_HH
#define TRACKER__READER_HH
#pragma once

#include <unordered_map>

#include "analysis.hh"

namespace MATHUSLA { namespace TRACKER {

namespace reader { /////////////////////////////////////////////////////////////////////////////

using namespace type;

namespace root { ///////////////////////////////////////////////////////////////////////////////

//__ROOT Directory Search_______________________________________________________________________
std::vector<std::string> search_directory(const std::string& path);
//----------------------------------------------------------------------------------------------

//__ROOT Import Types___________________________________________________________________________
using point_keys    = std::array<std::string, 4>;
using detector_keys = std::array<std::string, 2>;
using detector_map  = std::unordered_map<type::integer, std::string>;
//----------------------------------------------------------------------------------------------

//__ROOT Detector Map Import____________________________________________________________________
detector_map import_detector_map(const std::string& path);
//----------------------------------------------------------------------------------------------

//__ROOT Event Import___________________________________________________________________________
analysis::event_vector import_events(const std::string& path,
                                     const point_keys& keys);
analysis::event_vector import_events(const std::string& path,
                                     const detector_keys& keys,
                                     const detector_map& map);
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

namespace script { /////////////////////////////////////////////////////////////////////////////

//__Tracking Script Options Map_________________________________________________________________
using tracking_options = std::unordered_map<std::string, std::string>;
//----------------------------------------------------------------------------------------------

//__Allowed Keys for Tracking Script Options Map________________________________________________
static const std::array<std::string, 8> allowed_keys{{
  "geometry-file",
  "geometry-map",
  "root-data",
  "root-keys",
  "collapse-size",
  "layer-depth",
  "line-width",
  "seed-size"}};
//----------------------------------------------------------------------------------------------

//__Tracking Script Key Check___________________________________________________________________
bool is_key_allowed(const std::string& key);
//----------------------------------------------------------------------------------------------

//__Tracking Script Options Parser______________________________________________________________
const tracking_options read(const std::string& path);
//----------------------------------------------------------------------------------------------

} /* namespace script */ ///////////////////////////////////////////////////////////////////////

} /* namespace reader */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__READER_HH */

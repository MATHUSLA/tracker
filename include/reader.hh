/*
 * include/reader.hh
 *
 * Copyright 2018 Brandon Gomes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

//__ROOT Detector Map___________________________________________________________________________
using detector_map = std::unordered_map<integer, std::string>;
//----------------------------------------------------------------------------------------------

//__ROOT Detector Map Import____________________________________________________________________
detector_map import_detector_map(const std::string& path);
//----------------------------------------------------------------------------------------------

//__ROOT Event Import___________________________________________________________________________
analysis::event_vector import_events(const std::string& path,
                                     const std::string& time_key,
                                     const std::string& x_key,
                                     const std::string& y_key,
                                     const std::string& z_key);
analysis::event_vector import_events(const std::string& path,
                                     const std::string& time_key,
                                     const std::string& detector_key,
                                     const detector_map& map);
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

namespace script { /////////////////////////////////////////////////////////////////////////////

//__Tracking Script Options_____________________________________________________________________
struct tracking_options {
  std::string geometry_file     = "";
  std::string geometry_map_file = "";
  std::string root_directory    = "";
  std::string root_time_key     = "Time";
  std::string root_x_key        = "X";
  std::string root_y_key        = "Y";
  std::string root_z_key        = "Z";
  std::string root_detector_key = "Detector";
  r4_point    collapse_size     = {0, 0, 0, 0};
  real        layer_depth       = 500;
  real        line_width        = 1;
  integer     seed_size         = 3;
};
//----------------------------------------------------------------------------------------------

//__Tracking Script Options Parser______________________________________________________________
const tracking_options read(const std::string& path);
//----------------------------------------------------------------------------------------------

} /* namespace script */ ///////////////////////////////////////////////////////////////////////

} /* namespace reader */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__READER_HH */

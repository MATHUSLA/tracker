/*
 * include/tracker/reader.hh
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

#ifndef TRACKER__READER_HH
#define TRACKER__READER_HH
#pragma once

#include <tracker/analysis.hh>
#include <tracker/geometry.hh>
#include <tracker/units.hh>

namespace MATHUSLA { namespace TRACKER {

namespace reader { /////////////////////////////////////////////////////////////////////////////

using namespace type;

//__Data File Collection Mode___________________________________________________________________
enum class CollectionMode { Position, Detector };
//----------------------------------------------------------------------------------------------

//__Tracking Options Structure__________________________________________________________________
struct tracking_options {
  std::string    geometry_file      = "";
  std::string    geometry_map_file  = "";
  std::string    geometry_time_file = "";
  std::string    data_directory     = "";
  std::string    data_t_key         = "T";
  std::string    data_x_key         = "X";
  std::string    data_y_key         = "Y";
  std::string    data_z_key         = "Z";
  std::string    data_dt_key        = "dT";
  std::string    data_dx_key        = "dX";
  std::string    data_dy_key        = "dY";
  std::string    data_dz_key        = "dZ";
  std::string    data_detector_key  = "Detector";
  CollectionMode mode               = CollectionMode::Detector;
  real           default_time_error = 2 * units::time;
  r4_point       compression_size   = {0, 0, 0, 0};
  Coordinate     layer_axis         = Coordinate::Z;
  real           layer_depth        = 50 * units::length;
  real           line_width         = 1;
  size_t         seed_size          = 3;
};
//----------------------------------------------------------------------------------------------

//__Detector Map Import_________________________________________________________________________
const geometry::detector_map import_detector_map(const std::string& path);
//----------------------------------------------------------------------------------------------

//__Detector Time Resolution Map Import_________________________________________________________
const geometry::time_resolution_map import_time_resolution_map(const std::string& path);
//----------------------------------------------------------------------------------------------

namespace root { ///////////////////////////////////////////////////////////////////////////////

//__ROOT Directory Search_______________________________________________________________________
const std::vector<std::string> search_directory(const std::string& path,
                                                const std::string& ext="root");
//----------------------------------------------------------------------------------------------

//__ROOT Event Import___________________________________________________________________________
const analysis::event_vector import_events(const std::string& path,
                                           const std::string& t_key,
                                           const std::string& x_key,
                                           const std::string& y_key,
                                           const std::string& z_key);
const analysis::event_vector import_events(const std::string& path,
                                           const std::string& t_key,
                                           const std::string& detector_key,
                                           const geometry::detector_map& map);
const analysis::event_vector import_events(const std::string& path,
                                           const tracking_options& options,
                                           const geometry::detector_map& map);
const analysis::event_vector import_events(const std::string& path,
                                           const tracking_options& options);
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

namespace script { /////////////////////////////////////////////////////////////////////////////

//__Tracking Script Options Parser______________________________________________________________
const tracking_options read(const std::string& path);
//----------------------------------------------------------------------------------------------

} /* namespace script */ ///////////////////////////////////////////////////////////////////////

//__Parse Command Line Arguments________________________________________________________________
const tracking_options parse_input(int& argc,
                                   char* argv[]);
//----------------------------------------------------------------------------------------------

} /* namespace reader */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__READER_HH */

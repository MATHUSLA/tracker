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

#include <tracker/core/type.hh>
#include <tracker/core/units.hh>
#include <tracker/analysis/type.hh>
#include <tracker/analysis/monte_carlo.hh>
#include <tracker/geometry.hh>

namespace MATHUSLA { namespace TRACKER {

namespace reader { /////////////////////////////////////////////////////////////////////////////

using namespace type;

//__Tracking Options Structure__________________________________________________________________
struct tracking_options {
  std::string  geometry_file              = "";
  std::string  geometry_map_file          = "";
  std::string  geometry_time_file         = "";
  real         default_time_error         = 2 * units::time;

  std::string  data_directory             = "";
  std::string  data_file_extension        = "root";
  std::string  data_t_key                 = "";
  std::string  data_x_key                 = "";
  std::string  data_y_key                 = "";
  std::string  data_z_key                 = "";
  std::string  data_dt_key                = "";
  std::string  data_dx_key                = "";
  std::string  data_dy_key                = "";
  std::string  data_dz_key                = "";
  std::string  data_detector_key          = "Detector";
  std::string  data_track_id_key          = "Track";
  std::string  data_parent_id_key         = "Parent";
  std::string  data_e_key                 = "";
  std::string  data_px_key                = "";
  std::string  data_py_key                = "";
  std::string  data_pz_key                = "";

  std::string  statistics_directory       = "";
  std::string  statistics_file_prefix     = "statistics";
  std::string  statistics_file_extension  = "root";

  bool         time_smearing              = true;
  real         simulated_efficiency       = 1;
  real         simulated_noise_rate       = 0;
  real_range   event_time_window          = {0, 0};
  Coordinate   layer_axis                 = Coordinate::Z;
  real         layer_depth                = 0;
  real         line_width                 = 1;
  size_t       seed_size                  = 3;
  real         event_density_limit        = 1;
  real         event_overload_limit       = 2;
  real         track_density_limit        = 1;

  bool         verbose_output             = false;
  bool         draw_events                = false;
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

//__Import Mode Type____________________________________________________________________________
enum class ImportMode { Detector, Widths };
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
                                           const tracking_options& options,
                                           const ImportMode mode);
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Import______________________________________________________________________
const analysis::full_event_vector import_full_events(const std::string& path,
                                                     const std::string& t_key,
                                                     const std::string& x_key,
                                                     const std::string& y_key,
                                                     const std::string& z_key,
                                                     const std::string& dt_key,
                                                     const std::string& dx_key,
                                                     const std::string& dy_key,
                                                     const std::string& dz_key);
const analysis::full_event_vector import_full_events(const std::string& path,
                                                     const std::string& t_key,
                                                     const std::string& detector_key,
                                                     const geometry::detector_map& map);
const analysis::full_event_vector import_full_events(const std::string& path,
                                                     const tracking_options& options,
                                                     const geometry::detector_map& map);
const analysis::full_event_vector import_full_events(const std::string& path,
                                                     const tracking_options& options,
                                                     const ImportMode mode);
//----------------------------------------------------------------------------------------------

//__ROOT Event Import with Monte-Carlo Tracks___________________________________________________
const analysis::mc::event_vector_bundle import_event_mc_bundle(const std::string& path,
                                                               const std::string& track_key,
                                                               const std::string& t_key,
                                                               const std::string& x_key,
                                                               const std::string& y_key,
                                                               const std::string& z_key);
const analysis::mc::event_vector_bundle import_event_mc_bundle(const std::string& path,
                                                               const std::string& track_key,
                                                               const std::string& t_key,
                                                               const std::string& x_key,
                                                               const std::string& y_key,
                                                               const std::string& z_key,
                                                               const std::string& detector_key,
                                                               const geometry::detector_map& map);
const analysis::mc::event_vector_bundle import_event_mc_bundle(const std::string& path,
                                                               const tracking_options& options,
                                                               const geometry::detector_map& map);
const analysis::mc::event_vector_bundle import_event_mc_bundle(const std::string& path,
                                                               const tracking_options& options,
                                                               const ImportMode mode);
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Import with Monte-Carlo Tracks______________________________________________
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const std::string& path,
                                                                         const std::string& track_key,
                                                                         const std::string& t_key,
                                                                         const std::string& x_key,
                                                                         const std::string& y_key,
                                                                         const std::string& z_key,
                                                                         const std::string& dt_key,
                                                                         const std::string& dx_key,
                                                                         const std::string& dy_key,
                                                                         const std::string& dz_key);
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const std::string& path,
                                                                         const std::string& track_key,
                                                                         const std::string& t_key,
                                                                         const std::string& x_key,
                                                                         const std::string& y_key,
                                                                         const std::string& z_key,
                                                                         const std::string& dt_key,
                                                                         const std::string& detector_key,
                                                                         const geometry::detector_map& map);
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const std::string& path,
                                                                         const tracking_options& options,
                                                                         const geometry::detector_map& map);
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const std::string& path,
                                                                         const tracking_options& options,
                                                                         const ImportMode mode);
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

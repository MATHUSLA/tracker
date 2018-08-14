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

#include <fstream>

#include <tracker/core/type.hh>
#include <tracker/core/units.hh>
#include <tracker/analysis/type.hh>
#include <tracker/analysis/monte_carlo.hh>
#include <tracker/geometry.hh>

#include <tracker/util/command_line_parser.hh>
#include <tracker/util/error.hh>
#include <tracker/util/io.hh>
#include <tracker/util/string.hh>

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

namespace reserved { ///////////////////////////////////////////////////////////////////////////
//__Reserved Symbols____________________________________________________________________________
static const char comment_character           = '#';
static const char space_character             = ' ';
static const char key_value_separator         = ':';
static const std::string& continuation_string = "...";
//----------------------------------------------------------------------------------------------
} /* namespace reserved */ /////////////////////////////////////////////////////////////////////

//__Parse Line from Tracking Script_____________________________________________________________
void parse_line(const std::string& line,
                std::string& key,
                std::string& value);
//----------------------------------------------------------------------------------------------

//__Parse Key Value Pair into File Path_________________________________________________________
void parse_file_path(const std::string& key,
                     const std::string& value,
                     std::string& out,
                     bool exit_on_error=true);
//----------------------------------------------------------------------------------------------

//__Parse Data Key Value Pair___________________________________________________________________
void parse_data_keys(const std::string& key,
                     const std::string& value,
                     std::string& t_key,
                     std::string& x_key,
                     std::string& y_key,
                     std::string& z_key,
                     bool exit_on_error=true);
//----------------------------------------------------------------------------------------------

//__Parse Boolean Key Value_____________________________________________________________________
void parse_boolean(const std::string& key,
                   const std::string& value,
                   bool& out,
                   bool exit_on_error=true);
//----------------------------------------------------------------------------------------------

//__Parse Real Key Value________________________________________________________________________
void parse_real(const std::string& key,
                const std::string& value,
                real& out,
                bool exit_on_error=true);
//----------------------------------------------------------------------------------------------

//__Parse Positive Real Key Value_______________________________________________________________
void parse_positive_real(const std::string& key,
                         const std::string& value,
                         real& out,
                         bool exit_on_error=true);
//----------------------------------------------------------------------------------------------

//__Parse Real Range Key Value__________________________________________________________________
void parse_real_range(const std::string& key,
                      const std::string& value,
                      real_range& out,
                      bool exit_on_error=true);
//----------------------------------------------------------------------------------------------

//__Parse Integer Key Value_____________________________________________________________________
void parse_integer(const std::string& key,
                   const std::string& value,
                   integer& out,
                   bool exit_on_error=true);
//----------------------------------------------------------------------------------------------

//__Parse Size Type Key Value___________________________________________________________________
void parse_size_type(const std::string& key,
                     const std::string& value,
                     std::size_t& out,
                     bool exit_on_error=true);
//----------------------------------------------------------------------------------------------

//__Parse R4 Key Value__________________________________________________________________________
void parse_r4_point(const std::string& key,
                    const std::string& value,
                    r4_point& out,
                    bool use_units=true,
                    bool exit_on_error=true);
//----------------------------------------------------------------------------------------------

//__Parse R3 Coordinate Key Value_______________________________________________________________
void parse_r3_coordinate(const std::string& key,
                         const std::string& value,
                         Coordinate& coordinate,
                         bool exit_on_error=true);
//----------------------------------------------------------------------------------------------

//__Parse R4 Coordinate Key Value_______________________________________________________________
void parse_r4_coordinate(const std::string& key,
                         const std::string& value,
                         Coordinate& coordinate,
                         bool exit_on_error=true);
//----------------------------------------------------------------------------------------------

//__Read and Parse Lines from Tracking Script___________________________________________________
template<class TernaryFunction>
void parse_lines(std::ifstream& file,
                 TernaryFunction f) {
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) continue;
    std::string key, value;
    f(line, key, value);
  }
}
template<class TernaryFunction>
void parse_lines(const std::string& path,
                 TernaryFunction f) {
  std::ifstream file{path};
  parse_lines(file, f);
}
//----------------------------------------------------------------------------------------------

//__Read Lines from Tracking Script_____________________________________________________________
template<class BinaryFunction>
void read_lines(std::ifstream& file,
                BinaryFunction f) {
  parse_lines(file, [&](const auto& line, auto& key, auto& value) {
    parse_line(line, key, value);
    f(key, value);
  });
}
template<class BinaryFunction>
void read_lines(const std::string& path,
                BinaryFunction f) {
  std::ifstream file{path};
  read_lines(file, f);
}
//----------------------------------------------------------------------------------------------

//__Default Tracking Script Options Extension Parser____________________________________________
inline void default_extension_parser(const std::string& key,
                                     const std::string&,
                                     tracking_options&) {
  util::error::exit("[FATAL ERROR] Invalid Key in Tracking Script: \"", key, "\".\n");
}
//----------------------------------------------------------------------------------------------

//__Tracking Script Options Parser______________________________________________________________
template<class ExtensionParser=decltype(default_extension_parser)>
const tracking_options read(const std::string& path,
                            ExtensionParser& parser=default_extension_parser) {
  tracking_options out{};
  read_lines(path, [&](const auto& key, const auto& value) {
    if (key.empty()) return;
    if (value.empty()) {
      util::error::exit("[FATAL ERROR] Missing Value For Key: \"", key, "\".\n");
    } else {
      if (key == "geometry-file") {
        parse_file_path(key, value, out.geometry_file);
      } else if (key == "geometry-map-file") {
        parse_file_path(key, value, out.geometry_map_file);
      } else if (key == "geometry-map") {
        util::error::exit_when(value != reserved::continuation_string,
          "[FATAL ERROR] \"Geometry Map\" Key Requires Continuation String \"",
          reserved::continuation_string, "\" Before Map Entries.\n");
        // TODO: implement ...
      } else if (key == "data-directory") {
        parse_file_path(key, value, out.data_directory);
      } else if (key == "data-file-extension") {
        out.data_file_extension = value;
      } else if (key == "data-position-keys") {
        parse_data_keys(key, value,
          out.data_t_key, out.data_x_key, out.data_y_key, out.data_z_key);
      } else if (key == "data-position-error-keys") {
        parse_data_keys(key, value,
          out.data_dt_key, out.data_dx_key, out.data_dy_key, out.data_dz_key);
      } else if (key == "data-detector-key") {
        out.data_detector_key = value;
      } else if (key == "data-track-id-key") {
        out.data_track_id_key = value;
      } else if (key == "data-parent-id-key") {
        out.data_parent_id_key = value;
      } else if (key == "data-momentum-keys") {
        parse_data_keys(key, value,
          out.data_e_key, out.data_px_key, out.data_py_key, out.data_pz_key);
      } else if (key == "geometry-default-time-error") {
        parse_positive_real(key, value, out.default_time_error);
        out.default_time_error *= units::time;
      } else if (key == "time-smearing") {
        parse_boolean(key, value, out.time_smearing);
      } else if (key == "simulated-efficiency") {
        parse_positive_real(key, value, out.simulated_efficiency);
      } else if (key == "simulated-noise-rate") {
        parse_positive_real(key, value, out.simulated_noise_rate);
      } else if (key == "event-time-window") {
        parse_real_range(key, value, out.event_time_window);
      } else if (key == "layer-axis") {
        parse_r3_coordinate(key, value, out.layer_axis);
      } else if (key == "layer-depth") {
        parse_positive_real(key, value, out.layer_depth);
        out.layer_depth *= units::length;
      } else if (key == "line-width") {
        parse_positive_real(key, value, out.line_width);
        out.line_width *= units::length;
      } else if (key == "seed-size") {
        parse_size_type(key, value, out.seed_size);
      } else if (key == "event-density-limit") {
        parse_positive_real(key, value, out.event_density_limit);
      } else if (key == "event-overload-limit") {
        parse_positive_real(key, value, out.event_overload_limit);
      } else if (key == "track-density-limit") {
        parse_positive_real(key, value, out.track_density_limit);
      } else if (key == "statistics-directory") {
        parse_file_path(key, value, out.statistics_directory);
      } else if (key == "statistics-file-prefix") {
        out.statistics_file_prefix = value;
      } else if (key == "statistics-file-extension") {
        out.statistics_file_extension = value;
      } else if (key == "verbose-output") {
        parse_boolean(key, value, out.verbose_output);
      } else if (key == "draw-events") {
        parse_boolean(key, value, out.draw_events);
      } else {
        parser(key, value, out);
      }
    }
  });
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace script */ ///////////////////////////////////////////////////////////////////////

//__Parse Command Line Arguments________________________________________________________________
template<class ExtensionParser=decltype(script::default_extension_parser)>
const tracking_options parse_input(int& argc,
                                   char* argv[],
                                   ExtensionParser& parser=script::default_extension_parser) {
  using util::cli::option;
  option help_opt    ('h', "help",        "MATHUSLA Tracking Algorithm", option::no_arguments);
  option verbose_opt ('v', "verbose",     "Verbose Output",              option::no_arguments);
  option quiet_opt   ('q', "quiet",       "Quiet Output",                option::no_arguments);
  option event_opt   ( 0 , "draw-events", "Draw Events",                 option::no_arguments);
  option script_opt  ('s', "script",      "Tracking Script",             option::required_arguments);

  util::cli::parse(argv, {&help_opt, &verbose_opt, &quiet_opt, &event_opt, &script_opt});
  util::error::exit_when(!script_opt.count || argc == 1,
    "[FATAL ERROR] Insufficient Arguments:\n",
    "              Must include arguments for a tracking script (-s).\n");
  util::error::exit_when(!util::io::path_exists(script_opt.argument),
    "[FATAL ERROR] Tracking Script Missing: The file \"", script_opt.argument, "\" cannot be found.\n");

  auto out = script::read(script_opt.argument, parser);

  if (!quiet_opt.count) {
    out.verbose_output |= verbose_opt.count;
  } else {
    out.verbose_output = false;
  }

  if (event_opt.count)
    out.draw_events |= event_opt.count;

  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace reader */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__READER_HH */

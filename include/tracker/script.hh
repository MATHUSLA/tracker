/*
 * include/tracker/script.hh
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

#ifndef TRACKER__SCRIPT_HH
#define TRACKER__SCRIPT_HH
#pragma once

#include <fstream>

#include <tracker/core/type.hh>
#include <tracker/core/units.hh>

#include <tracker/util/command_line_parser.hh>
#include <tracker/util/error.hh>
#include <tracker/util/io.hh>
#include <tracker/util/string.hh>

namespace MATHUSLA { namespace TRACKER {

namespace script { /////////////////////////////////////////////////////////////////////////////

using namespace type;

//__Path Types__________________________________________________________________________________
using path_type = std::string;
using path_vector = std::vector<path_type>;
//----------------------------------------------------------------------------------------------

//__Tracking Options Structure__________________________________________________________________
struct tracking_options {
  path_type    geometry_file              = "";
  path_type    geometry_map_file          = "";
  path_type    geometry_time_file         = "";
  real         default_time_error         = 2 * units::time;

  path_vector  data_directories           = {""};
  real_vector  data_timing_offsets        = {0};
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

  path_type    statistics_directory       = "";
  path_type    statistics_file_prefix     = "statistics";
  std::string  statistics_file_extension  = "root";
  bool         merge_input                = false;

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

namespace reserved { ///////////////////////////////////////////////////////////////////////////
//__Reserved Symbols____________________________________________________________________________
static const char comment_character           = '#';
static const char space_character             = ' ';
static const char key_value_separator         = ':';
static const std::string& continuation_string = "...";
static const char continuation_line_character = '-';
//----------------------------------------------------------------------------------------------
} /* namespace reserved */ /////////////////////////////////////////////////////////////////////

//__Parse Line from Tracking Script_____________________________________________________________
void parse_line(const std::string& line,
                std::string& key,
                std::string& value);
//----------------------------------------------------------------------------------------------

//__Check that a Line has a Continuation String_________________________________________________
bool is_continuation_header(const std::string& key,
                            const std::string& value);
//----------------------------------------------------------------------------------------------

//__Check that a Line has a Continuation String_________________________________________________
bool is_continuation_body(const std::string& key,
                          const std::string& value);
//----------------------------------------------------------------------------------------------

//__Parse Line from a Continuation______________________________________________________________
void parse_continuation_line(const std::string& line,
                             std::string& entry);
//----------------------------------------------------------------------------------------------

//__Parse Key Value Pair into File Path_________________________________________________________
void parse_file_path(const std::string& key,
                     const std::string& value,
                     path_type& out,
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
template<class LineKeyValueParser>
void parse_lines(std::ifstream& file,
                 LineKeyValueParser f) {
  std::string line, key, value;
  bool continuation = false;
  while (std::getline(file, line)) {
    if (line.empty()) continue;
    continuation = f(line, continuation, key, value);
    if (!continuation)
      key.clear();
    value.clear();
  }
}
template<class LineKeyValueParser>
void parse_lines(const path_type& path,
                 LineKeyValueParser f) {
  std::ifstream file{path};
  parse_lines(file, f);
}
//----------------------------------------------------------------------------------------------

//__Read Lines from Tracking Script_____________________________________________________________
template<class KeyValueParser>
void read_lines(std::ifstream& file,
                KeyValueParser f) {
  parse_lines(file, [&](const auto& line, const auto continuation, auto& key, auto& value) {
    if (continuation) {
      std::string test_key, test_value;
      parse_line(line, test_key, test_value);
      if (!is_continuation_body(test_key, test_value)) {
        return f(test_key, test_value);
      } else {
        parse_continuation_line(line, value);
        return f(key, value);
      }
    } else {
      parse_line(line, key, value);
      return f(key, value);
    }
  });
}
template<class KeyValueParser>
void read_lines(const path_type& path,
                KeyValueParser f) {
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

//__Default Tracking Script Options Extension Parser Type_______________________________________
using default_extension_parser_t = decltype(default_extension_parser);
//----------------------------------------------------------------------------------------------

//__Tracking Script Options Parser______________________________________________________________
template<class ExtensionParser=default_extension_parser_t>
const tracking_options read(const path_type& path,
                            ExtensionParser& parser=default_extension_parser) {
  tracking_options out{};
  read_lines(path, [&](const auto& key, const auto& value) {
    if (key.empty()) return false;
    if (value.empty()) {
      util::error::exit("[FATAL ERROR] Missing Value For Key: \"", key, "\".\n");
    } else {
      if (key == "geometry-file") {
        parse_file_path(key, value, out.geometry_file);
      } else if (key == "geometry-map-file") {
        parse_file_path(key, value, out.geometry_map_file);
      } else if (key == "geometry-map") {
        if (is_continuation_header(key, value))
          return true;
        // TODO: implement
      } else if (key == "data-directory") {
        out.data_directories.clear();
        out.data_directories.emplace_back();
        parse_file_path(key, value, out.data_directories.back());
      } else if (key == "data-directories") {
        if (is_continuation_header(key, value)) {
          out.data_directories.clear();
          out.data_timing_offsets.clear();
          return true;
        }
        util::string_vector tokens;
        util::string::split(value, tokens, ",");
        util::error::exit_when(tokens.size() != 2UL,
          "[FATAL ERROR] Parallel Import Formatting for Data Directories is Invalid. \n",
          "              Expected \"path, offset\" but received \"", value, "\".\n");
        out.data_directories.emplace_back();
        auto& directories_back = out.data_directories.back();
        parse_file_path(key, tokens[0], directories_back);
        if (directories_back.empty()) {
          out.data_directories.pop_back();
        } else {
          out.data_timing_offsets.emplace_back();
          auto& offsets_back = out.data_timing_offsets.back();
          parse_real(tokens[0], tokens[1], offsets_back);
        }
        return true;
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
      } else if (key == "merge-input") {
        parse_boolean(key, value, out.merge_input);
      } else if (key == "verbose-output") {
        parse_boolean(key, value, out.verbose_output);
      } else if (key == "draw-events") {
        parse_boolean(key, value, out.draw_events);
      } else {
        parser(key, value, out);
      }
    }
    return false;
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Parse Command Line Arguments________________________________________________________________
template<class ExtensionParser=default_extension_parser_t>
const tracking_options parse_command_line(int& argc,
                                          char* argv[],
                                          ExtensionParser& parser=default_extension_parser) {
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

} /* namespace script */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__SCRIPT_HH */

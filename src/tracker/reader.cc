/*
 * src/tracker/reader.cc
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

#include <tracker/reader.hh>

#include <array>
#include <fstream>

#include <ROOT/TFile.h>
#include <ROOT/TKey.h>
#include <ROOT/TSystemDirectory.h>
#include <ROOT/TTree.h>

#include <tracker/geometry.hh>
#include <tracker/units.hh>

#include <tracker/util/command_line_parser.hh>
#include <tracker/util/error.hh>
#include <tracker/util/io.hh>
#include <tracker/util/string.hh>

namespace MATHUSLA { namespace TRACKER {

namespace reader { /////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Parse Line for Detector Map_________________________________________________________________
void _parse_detector_map_entry(const uint_fast64_t line_count,
                               std::string& line,
                               geometry::detector_map& out) {
  const auto colon = util::string::strip(line).find(":");
  if (colon == std::string::npos || line[0] == '#') {
    return;
  } else {
    try {
      auto i = std::stoll(line.substr(1 + colon));
      auto s = line.substr(0, colon);
      // std::cout << i << " " << s << "\n";
      out.insert({i, s});
    } catch (...) {
      util::error::exit("[FATAL ERROR] Detector Map has an invalid ID.\n",
                        "              See Line ", line_count, ".\n");
    }
  }
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Import Detector Map from File_______________________________________________________________
const geometry::detector_map import_detector_map(const std::string& path) {
  std::ifstream file(path);
  geometry::detector_map out{};
  std::string line;
  uint_fast64_t line_counter{};
  while (std::getline(file, line))
    _parse_detector_map_entry(++line_counter, line, out);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Detector Time Resolution Map Import_________________________________________________________
const geometry::time_resolution_map import_time_resolution_map(const std::string& path) {
  return {};
}
//----------------------------------------------------------------------------------------------

namespace root { ///////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Recursive TSystemFile Traversal_____________________________________________________________
void _collect_paths(TSystemDirectory* dir,
                    std::vector<std::string>& paths,
                    const std::string& ext) {
  if (!dir || !dir->GetListOfFiles()) return;
  for (const auto&& obj : *dir->GetListOfFiles()) {
    const auto file = static_cast<TSystemFile*>(obj);
    const auto name = std::string(file->GetName());
    if (!file->IsDirectory()) {
      if (ext == "" || name.substr(1 + name.find_last_of(".")) == ext)
        paths.push_back(std::string(file->GetTitle()) + "/" + name);
    } else if (!(name == "." || name == "..")) {
      _collect_paths(static_cast<TSystemDirectory*>(file), paths, ext);
    }
  }
}
//----------------------------------------------------------------------------------------------

//__ROOT File Key Traversal_____________________________________________________________________
template<class BinaryFunction>
BinaryFunction _traverse_file(const std::string& path,
                              BinaryFunction f) {
  auto file = TFile::Open(path.c_str());
  if (file && !file->IsZombie()) {
    TIter next(file->GetListOfKeys());
    TKey* key = nullptr;
    while ((key = static_cast<TKey*>(next())))
      f(file, key);
    file->Close();
  }
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Get TKey Classname__________________________________________________________________________
const std::string _get_TKey_classname(const TKey* key) {
  return key->GetClassName();
}
//----------------------------------------------------------------------------------------------

//__Get TTree Object From File Without Checking Type____________________________________________
TTree* _unchecked_get_TTree(TFile* file,
                            const TKey* key) {
  return static_cast<TTree*>(file->Get(key->GetName()));
}
//----------------------------------------------------------------------------------------------

//__Set TTree Branch Base Function______________________________________________________________
void _set_TTree_branches(TTree* tree,
                         const std::string& name,
                         Double_t* value) {
  tree->SetBranchAddress(name.c_str(), value);
}
//----------------------------------------------------------------------------------------------

//__Recursively Set TTree Branches______________________________________________________________
template<class... Args>
void _set_TTree_branches(TTree* tree,
                         const std::string& name,
                         Double_t* value,
                         Args ...args) {
  _set_TTree_branches(tree, name, value);
  _set_TTree_branches(tree, args...);
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Search and Collect ROOT File Paths__________________________________________________________
const std::vector<std::string> search_directory(const std::string& path,
                                                const std::string& ext) {
  std::vector<std::string> paths{};
  _collect_paths(new TSystemDirectory("search", path.c_str()), paths, ext);
  return paths;
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
const analysis::event_vector import_events(const std::string& path,
                                           const std::string& t_key,
                                           const std::string& x_key,
                                           const std::string& y_key,
                                           const std::string& z_key) {
  analysis::event_vector out{};
  TTree* tree = nullptr;
  _traverse_file(path, [&](const auto& file, const auto& key) {
    if (_get_TKey_classname(key) == "TTree") {
      Double_t t, x, y, z;
      tree = _unchecked_get_TTree(file, key);
      _set_TTree_branches(tree, t_key, &t, x_key, &x, y_key, &y, z_key, &z);
      const auto size = tree->GetEntries();
      analysis::event points;
      points.reserve(size);
      for (auto i = 0; i < size; ++i) {
        tree->GetEntry(i);
        points.push_back({t * units::time, x * units::length, y * units::length, z * units::length});
      }
      out.push_back(points);
    }
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
const analysis::event_vector import_events(const std::string& path,
                                           const std::string& t_key,
                                           const std::string& detector_key,
                                           const geometry::detector_map& map) {
  analysis::event_vector out{};
  TTree* tree = nullptr;
  _traverse_file(path, [&](const auto& file, const auto& key) {
    if (_get_TKey_classname(key) == "TTree") {
      Double_t t, detector;
      tree = _unchecked_get_TTree(file, key);
      _set_TTree_branches(tree, t_key, &t, detector_key, &detector);
      const auto size = tree->GetEntries();
      analysis::event points;
      points.reserve(size);
      for (auto i = 0; i < size; ++i) {
        tree->GetEntry(i);
        const auto& search = map.find(detector);
        const auto& name = search != map.end() ? search->second : std::to_string(std::llround(std::trunc(detector)));
        const auto center = geometry::limits_of(name).center;
        points.push_back({t * units::time, center.x, center.y, center.z});
      }
      out.push_back(points);
    }
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
const analysis::event_vector import_events(const std::string& path,
                                           const tracking_options& options,
                                           const geometry::detector_map& map) {
  return options.mode == reader::CollectionMode::Detector
    ? import_events(path, options.data_t_key, options.data_detector_key, map)
    : import_events(path, options.data_t_key, options.data_x_key, options.data_y_key, options.data_z_key);
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
const analysis::event_vector import_events(const std::string& path,
                                           const tracking_options& options) {
  return options.mode == reader::CollectionMode::Detector
    ? import_events(path, options.data_t_key, options.data_detector_key,
        import_detector_map(options.geometry_map_file))
    : import_events(path, options.data_t_key, options.data_x_key, options.data_y_key, options.data_z_key);
}
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

namespace script { /////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Allowed Keys for Tracking Script Options Map________________________________________________
const std::array<std::string, 11> _allowed_keys{{
  "geometry-file",
  "geometry-map-file",
  "geometry-map",
  "root-data",
  "root-position-keys",
  "root-detector-key",
  "compression-size",
  "layer-axis",
  "layer-depth",
  "line-width",
  "seed-size"}};
//----------------------------------------------------------------------------------------------

//__Reserved Symbols____________________________________________________________________________
const char _comment_character           = '#';
const char _space_character             = ' ';
const char _key_value_separator         = ':';
const std::string& _continuation_string = "...";
//----------------------------------------------------------------------------------------------

//__Parse Line from Script Into Key-Value Pair__________________________________________________
void _parse_script_line(const std::string& line,
                        std::string& key,
                        std::string& value) {
  bool on_key = true;
  for (const auto& ch : line) {
    if (ch == _comment_character) break;
    if (ch == _space_character) continue;
    if (ch == _key_value_separator) {
      on_key = false;
      continue;
    }
    if (on_key) key.push_back(ch);
    else value.push_back(ch);
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Key Value Pair into File Path_________________________________________________________
void _parse_key_value_file_path(const std::string& key,
                                const std::string& value,
                                std::string& out) {
  util::error::exit_when(!util::io::path_exists(value),
    "[FATAL ERROR] Missing File Path \"", value, "\" for Key \"", key, "\" in Tracking Script.\n");
  out = value;
}
//----------------------------------------------------------------------------------------------

//__Parse Data Key Value Pair___________________________________________________________________
void _parse_key_value_data_keys(const std::string& key,
                                const std::string& value,
                                CollectionMode& mode,
                                std::string& t_key,
                                std::string& x_key,
                                std::string& y_key,
                                std::string& z_key) {
  std::vector<std::string> data_keys;
  util::string::split(value, data_keys, ",;");
  const auto size = data_keys.size();
  const auto& end = data_keys.cend();
  util::error::exit_when(std::find(data_keys.cbegin(), end, "") != end,
    "[FATAL ERROR] Invalid Data Argument for \"", key, "\" in Tracking Script.\n"
    "              Expected 1 argument \"T\" or 4 arguments \"T, X, Y, Z\".\n");
  if (size == 1) {
    t_key = data_keys[0];
    mode = CollectionMode::Detector;
  } else if (size == 4) {
    t_key = data_keys[0];
    x_key = data_keys[1];
    y_key = data_keys[2];
    z_key = data_keys[3];
    mode = CollectionMode::Position;
  } else {
    util::error::exit(
      "[FATAL ERROR] Invalid Data Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected 1 argument \"T\" or 4 arguments \"T, X, Y, Z\" "
                    "but received ", size, ".\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Real Key Value________________________________________________________________________
void _parse_key_value_real(const std::string& key,
                           const std::string& value,
                           real& out) {
  try {
    out = static_cast<real>(std::stold(value));
  } catch (...) {
    util::error::exit(
      "[FATAL ERROR] Invalid Real Number Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible to floating point value.\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Positive Real Key Value________________________________________________________________________
void _parse_key_value_positive_real(const std::string& key,
                                    const std::string& value,
                                    real& out) {
  try {
    out = static_cast<real>(std::stold(value));
    if (out < 0) {
      out = 0;
      throw std::invalid_argument{"Real Number Must be Greater than or Equal to Zero."};
    }
  } catch (...) {
    util::error::exit(
      "[FATAL ERROR] Invalid Real Number Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible to positive floating point value.\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Integer Key Value_____________________________________________________________________
void _parse_key_value_integer(const std::string& key,
                              const std::string& value,
                              integer& out) {
  try {
    out = static_cast<integer>(std::stoll(value));
  } catch (...) {
    util::error::exit(
      "[FATAL ERROR] Invalid Integer Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible to integral value.\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Size Type Key Value___________________________________________________________________
void _parse_key_value_size_type(const std::string& key,
                                const std::string& value,
                                std::size_t& out) {
  try {
    out = static_cast<std::size_t>(std::stoull(value));
  } catch (...) {
    util::error::exit(
      "[FATAL ERROR] Invalid Size Type Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible to positive integral value.\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse R4 Key Value__________________________________________________________________________
void _parse_key_value_r4_point(const std::string& key,
                               const std::string& value,
                               r4_point& out) {
  std::vector<std::string> point;
  util::string::split(value, point, ",;");
  const auto size = point.size();
  const auto& end = point.cend();
  util::error::exit_when(size != 4 || (std::find(point.cbegin(), end, "") != end),
    "[FATAL ERROR] Invalid R4 Point Argument for \"", key, "\" in Tracking Script.\n"
    "              Expected 4 arguments \"T, X, Y, Z\" but received ", size, ".\n");
  try {
    out = {
      static_cast<real>(std::stold(point[0])),
      static_cast<real>(std::stold(point[1])),
      static_cast<real>(std::stold(point[2])),
      static_cast<real>(std::stold(point[3])) };
  } catch (...) {
    util::error::exit(
      "[FATAL ERROR] Invalid R4 Point Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected 4 arguments \"T, X, Y, Z\" convertible to floating point values.\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse R3 Coordinate Key Value_______________________________________________________________
void _parse_key_value_r3_coordinate(const std::string& key,
                                    const std::string& value,
                                    Coordinate& coordinate) {
  if (value == "x" || value == "X") {
    coordinate = Coordinate::X;
  } else if (value == "y" || value == "Y") {
    coordinate = Coordinate::Y;
  } else if (value == "z" || value == "Z") {
    coordinate = Coordinate::Z;
  } else {
    util::error::exit("[FATAL ERROR] Invalid R3 Coordinate Argument for \"", key, "\" in Tracking Script.\n",
                      "              Expected X/x, Y/y, or Z/z but received \"", value, "\".\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse R4 Coordinate Key Value_______________________________________________________________
void _parse_key_value_r4_coordinate(const std::string& key,
                                    const std::string& value,
                                    Coordinate& coordinate) {
  if (value == "t" || value == "T") {
    coordinate = Coordinate::T;
  } else if (value == "x" || value == "X") {
    coordinate = Coordinate::X;
  } else if (value == "y" || value == "Y") {
    coordinate = Coordinate::Y;
  } else if (value == "z" || value == "Z") {
    coordinate = Coordinate::Z;
  } else {
    util::error::exit("[FATAL ERROR] Invalid R4 Coordinate Argument for \"", key, "\" in Tracking Script.\n",
                      "              Expected T/t, X/x, Y/y, or Z/z but received \"", value, "\".\n");
  }
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Tracking Script Options Parser______________________________________________________________
const tracking_options read(const std::string& path) {
  std::ifstream file(path);
  tracking_options out;

  std::string line;
  // TODO: use for in-place map reader -> `uint_fast64_t line_counter{};`
  while (std::getline(file, line)) {
    if (line.empty()) continue;

    std::string key, value;
    _parse_script_line(line, key, value);

    if (key.empty()) continue;

    if (value.empty()) {
      util::error::exit("[FATAL ERROR] Missing Value For Key: \"", key, "\".\n");
    } else {
      if (key == "geometry-file") {
        _parse_key_value_file_path(key, value, out.geometry_file);
      } else if (key == "geometry-map-file") {
        _parse_key_value_file_path(key, value, out.geometry_map_file);
      } else if (key == "geometry-map") {
        util::error::exit_when(value != _continuation_string,
          "[FATAL ERROR] \"Geometry Map\" Key Requires Continuation String \"",
          _continuation_string, "\" Before Map Entries.\n");

      } else if (key == "root-data") {
        _parse_key_value_file_path(key, value, out.data_directory);
      } else if (key == "root-position-keys") {
        _parse_key_value_data_keys(key, value,
          out.mode, out.data_t_key, out.data_x_key, out.data_y_key, out.data_z_key);
      } else if (key == "root-position-error-keys") {
        _parse_key_value_data_keys(key, value,
          out.mode, out.data_dt_key, out.data_dx_key, out.data_dy_key, out.data_dz_key);
      } else if (key == "root-detector-key") {
        out.data_detector_key = value;
      } else if (key == "geometry-default-time-error") {
        _parse_key_value_positive_real(key, value, out.default_time_error);
      } else if (key == "compression-size") {
        _parse_key_value_r4_point(key, value, out.compression_size);
      } else if (key == "layer-axis") {
        _parse_key_value_r3_coordinate(key, value, out.layer_axis);
      } else if (key == "layer-depth") {
        _parse_key_value_positive_real(key, value, out.layer_depth);
      } else if (key == "line-width") {
        _parse_key_value_positive_real(key, value, out.line_width);
      } else if (key == "seed-size") {
        _parse_key_value_size_type(key, value, out.seed_size);
      } else {
        util::error::exit("[FATAL ERROR] Invalid Key in Tracking Script: \"", key, "\".\n");
      }
    }
  }
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace script */ ///////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Missing Path Exit Command___________________________________________________________________
void _exit_on_missing_path(const std::string& path,
                           const std::string& name) {
  util::error::exit_when(!util::io::path_exists(path),
    "[FATAL ERROR] ", name, " Missing: The file \"", path, "\" cannot be found.\n");
}
//----------------------------------------------------------------------------------------------
} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Parse Command Line Arguments________________________________________________________________
const tracking_options parse_input(int& argc,
                                   char* argv[]) {
  using util::cli::option;
  option help_opt    ('h', "help",     "MATHUSLA Tracking Algorithm", option::no_arguments);
  option verbose_opt ('v', "",         "Verbosity",                   option::required_arguments);
  option geo_opt     ('g', "geometry", "Geometry Import",             option::required_arguments);
  option map_opt     ('m', "map",      "Detector Map",                option::required_arguments);
  option data_opt    ('d', "data",     "ROOT Data Directory",         option::required_arguments);
  option script_opt  ('s', "script",   "Tracking Script",             option::required_arguments);

  util::cli::parse(argv, {&help_opt, &verbose_opt, &geo_opt, &data_opt, &map_opt, &script_opt});

  util::error::exit_when((geo_opt.count && !data_opt.count)
                      || (data_opt.count && !geo_opt.count)
                      || !(script_opt.count || geo_opt.count || data_opt.count)
                      || argc == 1,
    "[FATAL ERROR] Insufficient Arguments:\n",
    "              Must include arguments for geometry ",
                  "and ROOT directory (-gd) or arguments for a tracking script (-s).\n");

  if (script_opt.count) _exit_on_missing_path(script_opt.argument, "Tracking Script");

  reader::tracking_options out;

  if (script_opt.count) {
    out = script::read(script_opt.argument);
    geo_opt.count += !out.geometry_file.empty();
    map_opt.count += !out.geometry_map_file.empty();
    data_opt.count += !out.data_directory.empty();
  } else {
    out.geometry_file = geo_opt.count ? geo_opt.argument : "";
    out.geometry_map_file = map_opt.count ? map_opt.argument : "";
    out.data_directory = data_opt.count ? data_opt.argument : "";
  }

  if (geo_opt.count) _exit_on_missing_path(out.geometry_file, "Geometry File");
  if (map_opt.count) _exit_on_missing_path(out.geometry_map_file, "Geometry Map");
  if (data_opt.count) _exit_on_missing_path(out.data_directory, "ROOT Directory");

  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace reader */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

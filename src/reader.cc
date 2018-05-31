/*
 * src/reader.cc
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

#include "reader.hh"

#include <array>
#include <fstream>

#include <ROOT/TFile.h>
#include <ROOT/TKey.h>
#include <ROOT/TSystemDirectory.h>
#include <ROOT/TTree.h>

#include "geometry.hh"
#include "units.hh"
#include "util/command_line_parser.hh"
#include "util/error.hh"
#include "util/io.hh"
#include "util/string.hh"

namespace MATHUSLA { namespace TRACKER {

namespace reader { /////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Recursive TSystemFile Traversal_____________________________________________________________
void _collect_paths(TSystemDirectory* dir,
                    std::vector<std::string>& paths,
                    const std::string& ext) {
  if (!dir || !dir->GetListOfFiles()) return;
  for (const auto&& obj : *dir->GetListOfFiles()) {
    const auto file = static_cast<TSystemFile*>(obj);
    const auto&& name = std::string(file->GetName());
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
  TIter next(file->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = static_cast<TKey*>(next())))
    f(file, key);
  delete key;
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

namespace root { ///////////////////////////////////////////////////////////////////////////////

//__Search and Collect ROOT File Paths__________________________________________________________
std::vector<std::string> search_directory(const std::string& path) {
  std::vector<std::string> paths{};
  _collect_paths(new TSystemDirectory("data", path.c_str()), paths, "root");
  return paths;
}
//----------------------------------------------------------------------------------------------

//__Import Detector Map from File_______________________________________________________________
detector_map import_detector_map(const std::string& path,
                                 bool silence_errors) {
  std::ifstream file(path);
  detector_map out{};
  std::string line;
  uint_fast64_t line_counter;
  while (std::getline(file, line)) {
    ++line_counter;
    const auto& space = line.find(" ");
    if (space == std::string::npos || line[0] == '#') {
      continue;
    } else {
      try {
        out.insert({std::stoll(line.substr(1 + space)), line.substr(0, space)});
      } catch (...) {
        util::error::exit_when(!silence_errors,
          "[FATAL ERROR] Detector Map has an invalid ID.\n",
          "              See Line ", line_counter, ".\n");
      }
    }
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
analysis::event_vector import_events(const std::string& path,
                                     const std::string& time_key,
                                     const std::string& x_key,
                                     const std::string& y_key,
                                     const std::string& z_key) {
  analysis::event_vector out{};
  TTree* tree = nullptr;
  _traverse_file(path, [&](const auto& file, const auto& key) {
    if (std::string(key->GetClassName()) == "TTree") {
      Double_t t, x, y, z;
      tree = static_cast<TTree*>(file->Get(key->GetName()));
      tree->SetBranchAddress(time_key.c_str(), &t);
      tree->SetBranchAddress(x_key.c_str(), &x);
      tree->SetBranchAddress(y_key.c_str(), &y);
      tree->SetBranchAddress(z_key.c_str(), &z);
      const auto size = tree->GetEntries();
      analysis::event_points points;
      points.reserve(size);
      for (auto i = 0; i < size; ++i) {
        tree->GetEntry(i);
        points.push_back({
          t * units::time,
          x * units::length,
          y * units::length,
          z * units::length});
      }
      out.push_back(points);
    }
  });
  delete tree;
  return out;
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
analysis::event_vector import_events(const std::string& path,
                                     const std::string& time_key,
                                     const std::string& detector_key,
                                     const detector_map& map) {
  analysis::event_vector out{};
  TTree* tree = nullptr;
  _traverse_file(path, [&](const auto& file, const auto& key) {
    if (std::string(key->GetClassName()) == "TTree") {
      Double_t t, detector;
      tree = static_cast<TTree*>(file->Get(key->GetName()));
      tree->SetBranchAddress(time_key.c_str(), &t);
      tree->SetBranchAddress(detector_key.c_str(), &detector);
      const auto size = tree->GetEntries();
      analysis::event_points points;
      points.reserve(size);
      for (auto i = 0; i < size; ++i) {
        tree->GetEntry(i);
        const auto& search = map.find(detector);
        const auto& name = search != map.end() ? search->second : std::to_string(std::llround(std::trunc(detector)));
        const auto center = geometry::limits_of(name).center;
        points.push_back({t, center.x, center.y, center.z});
      }
      out.push_back(points);
    }
  });
  delete tree;
  return out;
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
analysis::event_vector import_events(const std::string& path,
                                     const tracking_options& options,
                                     const detector_map& map) {
 return options.mode == reader::CollectionMode::Detector
   ? import_events(path, options.root_time_key, options.root_detector_key, map)
   : import_events(path, options.root_time_key, options.root_x_key, options.root_y_key, options.root_z_key);
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
analysis::event_vector import_events(const std::string& path,
                                     const tracking_options& options) {
 const auto map = import_detector_map(options.geometry_map_file);
 return import_events(path, options, map);
}
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

namespace script { /////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Allowed Keys for Tracking Script Options Map________________________________________________
const std::array<std::string, 10> _allowed_keys{{
  "geometry-file",
  "geometry-map",
  "root-data",
  "root-position-keys",
  "root-detector-key",
  "collapse-size",
  "layer-axis",
  "layer-depth",
  "line-width",
  "seed-size"}};
//----------------------------------------------------------------------------------------------

//__Tracking Script Allowed Keys End Iterator___________________________________________________
const auto _allowed_keys_end = _allowed_keys.cend();
//----------------------------------------------------------------------------------------------

//__Tracking Script Key Check___________________________________________________________________
bool _is_key_allowed(const std::string& key) {
  return _allowed_keys_end != std::find(_allowed_keys.cbegin(), _allowed_keys_end, key);
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Tracking Script Options Parser______________________________________________________________
const tracking_options read(const std::string& path) {
  std::ifstream file(path);
  tracking_options out;

  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) continue;

    std::string key, value;
    bool on_key = true;
    for (const auto& ch : line) {
      if (ch == '#') break;
      if (ch == ' ') continue;
      if (ch == ':') {
        on_key = false;
        continue;
      }
      if (on_key) key.push_back(ch);
      else value.push_back(ch);
    }

    if (!value.empty() && _is_key_allowed(key)) {
      if (key == "geometry-file") {
        out.geometry_file = value;
      } else if (key == "geometry-map") {
        out.geometry_map_file = value;
      } else if (key == "root-data") {
        out.root_directory = value;

      } else if (key == "root-position-keys") {
        std::vector<std::string> root_keys;
        util::string::split(value, root_keys, ",;");
        const auto size = root_keys.size();
        const auto& end = root_keys.cend();
        util::error::exit_when(std::find(root_keys.cbegin(), end, "") != end,
          "[FATAL ERROR] \"ROOT Keys\" argument in Tracking Script is invalid.\n"
          "              Expected 1 arguments \"T\" or 4 arguments \"T, X, Y, Z\".\n");

        if (size == 1) {
          out.root_time_key = root_keys[0];
          out.mode = CollectionMode::Detector;
        } else if (size == 4) {
          out.root_time_key = root_keys[0];
          out.root_x_key = root_keys[1];
          out.root_y_key = root_keys[2];
          out.root_z_key = root_keys[3];
          out.mode = CollectionMode::Position;
        } else {
          util::error::exit(
            "[FATAL ERROR] \"ROOT Keys\" argument in Tracking Script is invalid.\n"
            "              Expected 1 arguments \"T\" or 4 arguments \"T, X, Y, Z\" "
                          "but received ", size, ".\n");
        }
      } else if (key == "root-detector-key") {
        out.root_detector_key = value;
      } else if (key == "collapse-size") {

        std::vector<std::string> point;
        util::string::split(value, point, ",;");
        const auto size = point.size();
        const auto& end = point.cend();

        util::error::exit_when(size != 4,
          "[FATAL ERROR] \"Collapse Size\" argument in Tracking Script is invalid.\n"
          "              Expected 4 arguments \"dt, dx, dy, dz\" but received ", size, ".\n");

        util::error::exit_when(std::find(point.cbegin(), end, "") != end,
          "[FATAL ERROR] \"Collapse Size\" argument in Tracking Script is invalid.\n"
          "              Expected 4 arguments \"dt, dx, dy, dz\".\n");

        out.collapse_size = {
          static_cast<real>(std::stold(point[0])),
          static_cast<real>(std::stold(point[1])),
          static_cast<real>(std::stold(point[2])),
          static_cast<real>(std::stold(point[3])) };
      } else if (key == "layer-axis") {
        if (value == "x" || value == "X") {
          out.layer_axis = Coordinate::X;
        } else if (value == "y" || value == "Y") {
          out.layer_axis = Coordinate::Y;
        } else if (value == "z" || value == "Z") {
          out.layer_axis = Coordinate::Z;
        } else {
          util::error::exit("[FATAL ERROR] \"Layer Axis\" argument in Tracking Script is invalid.\n",
                            "              Expected X/x, Y/y, or Z/z but received ", value, ".\n");
        }
      } else if (key == "layer-depth") {
        out.layer_depth = static_cast<real>(std::stold(value));
      } else if (key == "line-width") {
        out.line_width = static_cast<real>(std::stold(value));
      } else if (key == "seed-size") {
        out.seed_size = static_cast<size_t>(std::stoull(value));
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
    "[FATAL ERROR] ", name, " Missing: The file ", path, " cannot be found.\n");
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
  option root_opt    ('d', "data",     "ROOT Data Directory",         option::required_arguments);
  option script_opt  ('s', "script",   "Tracking Script",             option::required_arguments);

  util::cli::parse(argv, {&help_opt, &verbose_opt, &geo_opt, &root_opt, &map_opt, &script_opt});

  util::error::exit_when((geo_opt.count && !root_opt.count)
                      || (root_opt.count && !geo_opt.count)
                      || !(script_opt.count || geo_opt.count || root_opt.count)
                      || argc == 1,
    "[FATAL ERROR] Insufficient Arguments:\n",
    "              Must include arguments for geometry ",
    "and ROOT directory \"(-gd)\" or arguments for a tracking script \"(-s)\".\n");

  if (script_opt.count) _exit_on_missing_path(script_opt.argument, "Tracking Script");

  reader::tracking_options out;

  if (script_opt.count) {
    out = script::read(script_opt.argument);
    geo_opt.count += !out.geometry_file.empty();
    map_opt.count += !out.geometry_map_file.empty();
    root_opt.count += !out.root_directory.empty();
  } else {
    out.geometry_file = geo_opt.count ? geo_opt.argument : "";
    out.geometry_map_file = map_opt.count ? map_opt.argument : "";
    out.root_directory = root_opt.count ? root_opt.argument : "";
  }

  if (geo_opt.count) _exit_on_missing_path(out.geometry_file, "Geometry File");
  if (map_opt.count) _exit_on_missing_path(out.geometry_map_file, "Geometry Map");
  if (root_opt.count) _exit_on_missing_path(out.root_directory, "ROOT Directory");

  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace reader */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

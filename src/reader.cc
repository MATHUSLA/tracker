/*
 * src/reader.cc
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

#include "reader.hh"

#include <array>
#include <fstream>

#include "ROOT/TFile.h"
#include "ROOT/TKey.h"
#include "ROOT/TSystemDirectory.h"
#include "ROOT/TTree.h"

#include "geometry.hh"
#include "units.hh"
#include "util/error.hh"
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
BinaryFunction _traverse_file(const std::string& path, BinaryFunction f) {
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
detector_map import_detector_map(const std::string& path) {
  std::ifstream file(path);
  detector_map out{};
  std::string line;
  while (std::getline(file, line)) {
    // TODO: how about std::isspace?
    const auto& space = line.find(" ");
    if (space == std::string::npos || line[0] == '#') {
      continue;
    } else {
      try {
        const auto& name = line.substr(0, space);
        const auto& id = std::stoll(line.substr(1 + space));
        out.insert({id, name}); // overwrites duplicate names
      } catch (...) {}
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
      analysis::event_points points{};
      const auto&& size = tree->GetEntries();
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
      analysis::event_points points{};
      const auto&& size = tree->GetEntries();
      for (auto i = 0; i < size; ++i) {
        tree->GetEntry(i);
        const auto& search = map.find(detector);
        const auto& name = search != map.end() ? search->second : std::to_string(detector);
        const auto& center = geometry::limits_of(name).center;
        points.push_back({t, center.x, center.y, center.z});
      }
      out.push_back(points);
    }
  });
  delete tree;
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

namespace script { /////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Allowed Keys for Tracking Script Options Map________________________________________________
const std::array<std::string, 8> _allowed_keys{{
  "geometry-file",
  "geometry-map",
  "root-data",
  "root-keys",
  "collapse-size",
  "layer-depth",
  "line-width",
  "seed-size"}};
//----------------------------------------------------------------------------------------------

//__Tracking Script Allowed Keys End Iterator___________________________________________________
const auto& _allowed_keys_end = _allowed_keys.cend();
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
      } else if (key == "root-keys") {

        std::vector<std::string> root_keys;
        util::string::split(value, root_keys, ",;");
        const auto& size = root_keys.size();
        const auto& end = root_keys.cend();

        util::error::exit_when(std::find(root_keys.cbegin(), end, "") != end,
          "[FATAL ERROR] \"ROOT Keys\" argument in Tracking Script is invalid.\n"
          "              Expected 2 arguments \"Time, Detector\" or 4 arguments \"Time, X, Y, Z\".\n");

        if (size == 2) {
          out.root_time_key = root_keys[0];
          out.root_detector_key = root_keys[1];
        } else if (size == 4) {
          out.root_time_key = root_keys[0];
          out.root_x_key = root_keys[1];
          out.root_y_key = root_keys[2];
          out.root_z_key = root_keys[3];
        } else {
          util::error::exit(
            "[FATAL ERROR] \"ROOT Keys\" argument in Tracking Script is invalid.\n"
            "              Expected 2 arguments \"Time, Detector\" or 4 arguments \"Time, X, Y, Z\" "
                          "but received ", size, ".\n");
        }

      } else if (key == "collapse-size") {

        std::vector<std::string> point;
        util::string::split(value, point, ",;");
        const auto& size = point.size();
        const auto& end = point.cend();

        util::error::exit_when(size != 4,
          "[FATAL ERROR] \"Collapse Size\" argument in Tracking Script is invalid.\n"
          "              Expected 4 arguments \"dt, dx, dy, dz\" but received ", size, ".\n");

        util::error::exit_when(std::find(point.cbegin(), end, "") != end,
          "[FATAL ERROR] \"Collapse Size\" argument in Tracking Script is invalid.\n"
          "              Expected 4 arguments \"dt, dx, dy, dz\".\n");

        out.collapse_size = {
          static_cast<type::real>(std::stold(point[0])),
          static_cast<type::real>(std::stold(point[1])),
          static_cast<type::real>(std::stold(point[2])),
          static_cast<type::real>(std::stold(point[3])) };

      } else if (key == "layer-depth") {
        out.layer_depth = static_cast<type::real>(std::stold(value));
      } else if (key == "line-width") {
        out.line_width = static_cast<type::real>(std::stold(value));
      } else if (key == "seed-size") {
        out.seed_size = static_cast<type::integer>(std::stoull(value));
      }
    }
  }
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace script */ ///////////////////////////////////////////////////////////////////////

} /* namespace reader */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

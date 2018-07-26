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

#include <tracker/core/units.hh>

#include <tracker/util/command_line_parser.hh>
#include <tracker/util/error.hh>
#include <tracker/util/io.hh>
#include <tracker/util/string.hh>

#include "helper/root.hh"

namespace MATHUSLA { namespace TRACKER {

namespace reader { /////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Parse Line for Detector Map_________________________________________________________________
void _parse_detector_map_entry(const uint_fast64_t line_count,
                               std::string& line,
                               geometry::detector_map& out) {
  const auto colon = util::string::strip(line).find(':');
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

//__Parse Line for Time Resolution Map__________________________________________________________
void _parse_time_resolution_map_entry(const uint_fast64_t line_count,
                                      std::string& line,
                                      geometry::time_resolution_map& out) {
  const auto colon = util::string::strip(line).find(':');
  if (colon == std::string::npos || line[0] == '#') {
    return;
  } else {
    try {
      auto i = std::stold(line.substr(1 + colon));
      auto s = line.substr(0, colon);
      // std::cout << i << " " << s << "\n";
      out.insert({s, i});
    } catch (...) {
      util::error::exit("[FATAL ERROR] Time Resolution Map has an invalid Time.\n",
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
  std::ifstream file(path);
  geometry::time_resolution_map out{};
  std::string line;
  uint_fast64_t line_counter{};
  while (std::getline(file, line))
    _parse_time_resolution_map_entry(++line_counter, line, out);
  return out;
}
//----------------------------------------------------------------------------------------------

namespace root { ///////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Helper Namespace Import_____________________________________________________________________
using namespace ::MATHUSLA::TRACKER::root;
//----------------------------------------------------------------------------------------------

//__Search for Name in Detector Map or Use ID as Name___________________________________________
const std::string _detector_name(const Double_t detector,
                                 const geometry::detector_map& map) {
  const auto& search = map.find(detector);
  return search != map.end() ? search->second : std::to_string(std::llround(std::trunc(detector)));
}
//----------------------------------------------------------------------------------------------

//__Convert Track Data to Unsigned ID___________________________________________________________
std::size_t _track_id(const Double_t track) {
  return static_cast<std::size_t>(std::llround(track));
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Search and Collect ROOT File Paths__________________________________________________________
const std::vector<std::string> search_directory(const std::string& path,
                                                const std::string& ext) {
  helper::init();
  return helper::search_directory(path, ext);
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
const analysis::event_vector import_events(const std::string& path,
                                           const std::string& t_key,
                                           const std::string& x_key,
                                           const std::string& y_key,
                                           const std::string& z_key) {
  analysis::event_vector out{};
  helper::traverse_keys(path, "READ", "TTree", [&](const auto& file, const auto& key) {
    auto tree = helper::tree::get_tree(file, key);
    if (!tree)
      return;

    helper::tree::vector_data_type* t = nullptr;
    helper::tree::vector_data_type* x = nullptr;
    helper::tree::vector_data_type* y = nullptr;
    helper::tree::vector_data_type* z = nullptr;
    helper::tree::set_branches(tree, t_key, &t, x_key, &x, y_key, &y, z_key, &z);

    helper::tree::traverse_entries(tree,
      [&](const auto entries) {
        out.reserve(entries);
      }, [&](){
        const auto size = t->size();
        if (size == 0)
          return;

        analysis::event points;
        points.reserve(size);
        for (std::size_t i = 0; i < size; ++i) {
          points.push_back({
            (*t)[i] * units::time, (*x)[i] * units::length, (*y)[i] * units::length, (*z)[i] * units::length});
        }
        out.push_back(points);
      });
    return; // FIXME: find other way to get only one TTree
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
  helper::traverse_keys(path, "READ", "TTree", [&](const auto& file, const auto& key) {
    auto tree = helper::tree::get_tree(file, key);
    if (!tree)
      return;

    helper::tree::vector_data_type* detector = nullptr;
    helper::tree::vector_data_type* t        = nullptr;
    helper::tree::set_branches(tree, t_key, &t, detector_key, &detector);

    helper::tree::traverse_entries(tree,
      [&](const auto entries) {
        out.reserve(entries);
      }, [&](){
        const auto size = detector->size();
        if (size == 0)
          return;

        analysis::event points;
        points.reserve(size);

        for (std::size_t i = 0; i < size; ++i) {
          const auto center = geometry::limits_of(_detector_name((*detector)[i], map)).center;
          points.push_back({(*t)[i] * units::time, center.x, center.y, center.z});
        }
        out.push_back(points);
      });
    return; // FIXME: find other way to get only one TTree
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
const analysis::event_vector import_events(const std::string& path,
                                           const tracking_options& options,
                                           const geometry::detector_map& map) {
  return import_events(path, options.data_t_key, options.data_detector_key, map);
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
const analysis::event_vector import_events(const std::string& path,
                                           const tracking_options& options,
                                           const ImportMode mode) {
  return mode == ImportMode::Detector
    ? import_events(path, options, import_detector_map(options.geometry_map_file))
    : import_events(path, options.data_t_key, options.data_x_key, options.data_y_key, options.data_z_key);
}
//----------------------------------------------------------------------------------------------

//__Import Full Events from ROOT File___________________________________________________________
const analysis::full_event_vector import_full_events(const std::string& path,
                                                     const std::string& t_key,
                                                     const std::string& x_key,
                                                     const std::string& y_key,
                                                     const std::string& z_key,
                                                     const std::string& dt_key,
                                                     const std::string& dx_key,
                                                     const std::string& dy_key,
                                                     const std::string& dz_key) {
  analysis::full_event_vector out{};
  helper::traverse_keys(path, "READ", "TTree", [&](const auto& file, const auto& key) {
    auto tree = helper::tree::get_tree(file, key);
    if (!tree)
      return;

    helper::tree::vector_data_type* t = nullptr;
    helper::tree::vector_data_type* x = nullptr;
    helper::tree::vector_data_type* y = nullptr;
    helper::tree::vector_data_type* z = nullptr;
    helper::tree::vector_data_type* dt = nullptr;
    helper::tree::vector_data_type* dx = nullptr;
    helper::tree::vector_data_type* dy = nullptr;
    helper::tree::vector_data_type* dz = nullptr;
    helper::tree::set_branches(tree,
      t_key, &t, x_key, &x, y_key, &y, z_key, &z,
      dt_key, &dt, dx_key, &dx, dy_key, &dy, dz_key, &dz);

    helper::tree::traverse_entries(tree,
      [&](const auto entries) {
        out.reserve(entries);
      }, [&](){
        const auto size = t->size();
        if (size == 0)
          return;

        analysis::full_event points;
        points.reserve(size);

        for (std::size_t i = 0; i < size; ++i) {
          points.push_back({
            (*t)[i] * units::time, (*x)[i] * units::length, (*y)[i] * units::length, (*z)[i] * units::length,
            {(*dt)[i] * units::time, (*dx)[i] * units::length, (*dy)[i] * units::length, (*dz)[i] * units::length}});
        }
        out.push_back(points);
      });
    return; // FIXME: find other way to get only one TTree
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Import Full Events from ROOT File___________________________________________________________
const analysis::full_event_vector import_full_events(const std::string& path,
                                                     const std::string& t_key,
                                                     const std::string& dt_key,
                                                     const std::string& detector_key,
                                                     const geometry::detector_map& map) {
  analysis::full_event_vector out{};
  helper::traverse_keys(path, "READ", "TTree", [&](const auto& file, const auto& key) {
    auto tree = helper::tree::get_tree(file, key);
    if (!tree)
      return;

    helper::tree::vector_data_type* detector = nullptr;
    helper::tree::vector_data_type* t        = nullptr;
    helper::tree::vector_data_type* dt       = nullptr;
    helper::tree::set_branches(tree, t_key, &t, dt_key, &dt, detector_key, &detector);

    helper::tree::traverse_entries(tree,
      [&](const auto entries) {
        out.reserve(entries);
      }, [&](){
        const auto size = detector->size();
        if (size == 0)
          return;

        analysis::full_event points;
        points.reserve(size);

        for (std::size_t i = 0; i < size; ++i) {
          const auto limits = geometry::limits_of(_detector_name((*detector)[i], map));
          const auto& center = limits.center;
          const auto& min = limits.min;
          const auto& max = limits.max;
          points.push_back({
            (*t)[i] * units::time, center.x, center.y, center.z,
            {(*dt)[i] * units::time, max.x - min.x, max.y - min.y, max.z - min.z}});
        }
        out.push_back(points);
      });
    return; // FIXME: find other way to get only one TTree
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Import Full Events from ROOT File___________________________________________________________
const analysis::full_event_vector import_full_events(const std::string& path,
                                                     const tracking_options& options,
                                                     const geometry::detector_map& map) {
  return import_full_events(path, options.data_t_key, options.data_dt_key, options.data_detector_key, map);
}
//----------------------------------------------------------------------------------------------

//__Import Full Events from ROOT File___________________________________________________________
const analysis::full_event_vector import_full_events(const std::string& path,
                                                     const tracking_options& options,
                                                     const ImportMode mode) {
  return mode == ImportMode::Detector
    ? import_full_events(path, options, import_detector_map(options.geometry_map_file))
    : import_full_events(path,
        options.data_t_key,
        options.data_x_key,
        options.data_y_key,
        options.data_z_key,
        options.data_dt_key,
        options.data_dx_key,
        options.data_dy_key,
        options.data_dz_key);
}
//----------------------------------------------------------------------------------------------

//__ROOT Event Import with Monte-Carlo Tracks___________________________________________________
const analysis::mc::event_vector_bundle import_event_mc_bundle(const std::string& path,
                                                               const std::string& track_key,
                                                               const std::string& t_key,
                                                               const std::string& x_key,
                                                               const std::string& y_key,
                                                               const std::string& z_key) {
  analysis::mc::event_vector_bundle out{{}, {}};
  helper::traverse_keys(path, "READ", "TTree", [&](const auto& file, const auto& key) {
    auto tree = helper::tree::get_tree(file, key);
    if (!tree)
      return;

    helper::tree::vector_data_type* track = nullptr;
    helper::tree::vector_data_type* t     = nullptr;
    helper::tree::vector_data_type* x     = nullptr;
    helper::tree::vector_data_type* y     = nullptr;
    helper::tree::vector_data_type* z     = nullptr;
    helper::tree::set_branches(tree,
      track_key, &track, t_key, &t, x_key, &x, y_key, &y, z_key, &z);

    helper::tree::traverse_entries(tree,
      [&](const auto entries) {
        out.events.reserve(entries);
        out.true_events.reserve(entries);
      }, [&](){
        const auto size = track->size();
        if (size == 0)
          return;

        analysis::event points;
        analysis::mc::event true_points;
        points.reserve(size);
        true_points.reserve(size);

        for (std::size_t i = 0; i < size; ++i) {
          true_points.push_back({
            _track_id((*track)[i]),
            (*t)[i] * units::time, (*x)[i] * units::length, (*y)[i] * units::length, (*z)[i] * units::length});
          points.push_back(reduce_to_r4(true_points.back()));
        }
        out.events.push_back(points);
        out.true_events.push_back(true_points);
      });
    return; // FIXME: find other way to get only one TTree
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__ROOT Event Import with Monte-Carlo Tracks___________________________________________________
const analysis::mc::event_vector_bundle import_event_mc_bundle(const std::string& path,
                                                               const std::string& track_key,
                                                               const std::string& t_key,
                                                               const std::string& x_key,
                                                               const std::string& y_key,
                                                               const std::string& z_key,
                                                               const std::string& detector_key,
                                                               const geometry::detector_map& map) {
  analysis::mc::event_vector_bundle out{{}, {}};
  helper::traverse_keys(path, "READ", "TTree", [&](const auto& file, const auto& key) {
    auto tree = helper::tree::get_tree(file, key);
    if (!tree)
      return;

    helper::tree::vector_data_type* track     = nullptr;
    helper::tree::vector_data_type* detector  = nullptr;
    helper::tree::vector_data_type* t         = nullptr;
    helper::tree::vector_data_type* x         = nullptr;
    helper::tree::vector_data_type* y         = nullptr;
    helper::tree::vector_data_type* z         = nullptr;
    helper::tree::set_branches(tree,
      track_key, &track, detector_key, &detector, t_key, &t, x_key, &x, y_key, &y, z_key, &z);

    helper::tree::traverse_entries(tree,
      [&](const auto entries) {
        out.events.reserve(entries);
        out.true_events.reserve(entries);
      }, [&](){
        const auto size = track->size();
        if (size == 0)
          return;

        analysis::event points;
        analysis::mc::event true_points;
        points.reserve(size);
        true_points.reserve(size);

        for (std::size_t i = 0; i < size; ++i) {
          const auto center = geometry::limits_of(_detector_name((*detector)[i], map)).center;
          points.push_back({(*t)[i] * units::time, center.x, center.y, center.z});
          true_points.push_back({
            _track_id((*track)[i]),
            (*t)[i] * units::time, (*x)[i] * units::length, (*y)[i] * units::length, (*z)[i] * units::length});
        }
        out.events.push_back(points);
        out.true_events.push_back(true_points);
      });
    return; // FIXME: find other way to get only one TTree
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__ROOT Event Import with Monte-Carlo Tracks___________________________________________________
const analysis::mc::event_vector_bundle import_event_mc_bundle(const std::string& path,
                                                               const tracking_options& options,
                                                               const geometry::detector_map& map) {
  return import_event_mc_bundle(path,
    options.data_track_id_key,
    options.data_t_key,
    options.data_x_key,
    options.data_y_key,
    options.data_z_key,
    options.data_detector_key,
    map);
}
//----------------------------------------------------------------------------------------------

//__ROOT Event Import with Monte-Carlo Tracks___________________________________________________
const analysis::mc::event_vector_bundle import_event_mc_bundle(const std::string& path,
                                                               const tracking_options& options,
                                                               const ImportMode mode) {
  return mode == ImportMode::Detector
    ? import_event_mc_bundle(path, options, import_detector_map(options.geometry_map_file))
    : import_event_mc_bundle(path,
        options.data_track_id_key,
        options.data_t_key,
        options.data_x_key,
        options.data_y_key,
        options.data_z_key);
}
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
                                                                         const std::string& dz_key) {
  analysis::mc::full_event_vector_bundle out{{}, {}};
  helper::traverse_keys(path, "READ", "TTree", [&](const auto& file, const auto& key) {
    auto tree = helper::tree::get_tree(file, key);
    if (!tree)
      return;

    helper::tree::vector_data_type* track = nullptr;
    helper::tree::vector_data_type* t = nullptr;
    helper::tree::vector_data_type* x = nullptr;
    helper::tree::vector_data_type* y = nullptr;
    helper::tree::vector_data_type* z = nullptr;
    helper::tree::vector_data_type* dt = nullptr;
    helper::tree::vector_data_type* dx = nullptr;
    helper::tree::vector_data_type* dy = nullptr;
    helper::tree::vector_data_type* dz = nullptr;
    helper::tree::set_branches(tree, track_key, &track,
      t_key, &t, x_key, &x, y_key, &y, z_key, &z,
      dt_key, &dt, dx_key, &dx, dy_key, &dy, dz_key, &dz);

    helper::tree::traverse_entries(tree,
      [&](const auto entries) {
        out.events.reserve(entries);
        out.true_events.reserve(entries);
      }, [&](){
        const auto size = track->size();
        if (size == 0)
          return;

        analysis::full_event points;
        analysis::mc::event true_points;
        points.reserve(size);
        true_points.reserve(size);

        for (std::size_t i = 0; i < size; ++i) {
          points.push_back({
            (*t)[i] * units::time, (*x)[i] * units::length, (*y)[i] * units::length, (*z)[i] * units::length,
            {(*dt)[i] * units::time, (*dx)[i] * units::length, (*dy)[i] * units::length, (*dz)[i] * units::length}});
          true_points.push_back({
            _track_id((*track)[i]),
            (*t)[i] * units::time, (*x)[i] * units::length, (*y)[i] * units::length, (*z)[i] * units::length});
        }
        out.events.push_back(points);
        out.true_events.push_back(true_points);
      });
    return; // FIXME: find other way to get only one TTree
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Import with Monte-Carlo Tracks______________________________________________
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const std::string& path,
                                                                         const std::string& track_key,
                                                                         const std::string& t_key,
                                                                         const std::string& x_key,
                                                                         const std::string& y_key,
                                                                         const std::string& z_key,
                                                                         const std::string& dt_key,
                                                                         const std::string& detector_key,
                                                                         const geometry::detector_map& map) {
  analysis::mc::full_event_vector_bundle out{{}, {}};
  helper::traverse_keys(path, "READ", "TTree", [&](const auto& file, const auto& key) {
    auto tree = helper::tree::get_tree(file, key);
    if (!tree)
      return;

    helper::tree::vector_data_type* track = nullptr;
    helper::tree::vector_data_type* detector = nullptr;
    helper::tree::vector_data_type* t = nullptr;
    helper::tree::vector_data_type* x = nullptr;
    helper::tree::vector_data_type* y = nullptr;
    helper::tree::vector_data_type* z = nullptr;
    helper::tree::vector_data_type* dt = nullptr;
    helper::tree::set_branches(tree,
      track_key, &track, detector_key, &detector,
      dt_key, &dt, t_key, &t, x_key, &x, y_key, &y, z_key, &z);

    helper::tree::traverse_entries(tree,
      [&](const auto entries) {
        out.events.reserve(entries);
        out.true_events.reserve(entries);
      }, [&](){
        const auto size = track->size();
        if (size == 0)
          return;

        analysis::full_event points;
        analysis::mc::event true_points;
        points.reserve(size);
        true_points.reserve(size);

        for (std::size_t i = 0; i < size; ++i) {
          const auto limits = geometry::limits_of(_detector_name((*detector)[i], map));
          const auto& center = limits.center;
          const auto& min = limits.min;
          const auto& max = limits.max;
          points.push_back({
            (*t)[i] * units::time, center.x, center.y, center.z,
            {(*dt)[i] * units::time, max.x - min.x, max.y - min.y, max.z - min.z}});
          true_points.push_back({
            _track_id((*track)[i]),
            (*t)[i] * units::time, (*x)[i] * units::length, (*y)[i] * units::length, (*z)[i] * units::length});
        }
        out.events.push_back(points);
        out.true_events.push_back(true_points);
      });
    return; // FIXME: find other way to get only one TTree
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Import with Monte-Carlo Tracks______________________________________________
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const std::string& path,
                                                                         const tracking_options& options,
                                                                         const geometry::detector_map& map) {
  return import_full_event_mc_bundle(path,
    options.data_track_id_key,
    options.data_t_key,
    options.data_x_key,
    options.data_y_key,
    options.data_z_key,
    options.data_dt_key,
    options.data_detector_key,
    map);
}
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Import with Monte-Carlo Tracks______________________________________________
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const std::string& path,
                                                             const tracking_options& options,
                                                             const ImportMode mode) {
  return mode == ImportMode::Detector
    ? import_full_event_mc_bundle(path, options, import_detector_map(options.geometry_map_file))
    : import_full_event_mc_bundle(path,
        options.data_track_id_key,
        options.data_t_key,
        options.data_x_key,
        options.data_y_key,
        options.data_z_key,
        options.data_dt_key,
        options.data_dx_key,
        options.data_dy_key,
        options.data_dz_key);
}
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

namespace script { /////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

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
  } else if (size == 4) {
    t_key = data_keys[0];
    x_key = data_keys[1];
    y_key = data_keys[2];
    z_key = data_keys[3];
  } else {
    util::error::exit(
      "[FATAL ERROR] Invalid Data Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected 1 argument \"T\" or 4 arguments \"T, X, Y, Z\" "
                    "but received ", size, ".\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Boolean Key Value_____________________________________________________________________
void _parse_key_value_boolean(const std::string& key,
                              const std::string& value,
                              bool& out) {
  const auto value_lowercase = util::string::tolower(value);
  if (value_lowercase == "true" || value_lowercase == "t" || value_lowercase == "1") {
    out = true;
  } else if (value_lowercase == "false" || value_lowercase == "f" || value_lowercase == "0") {
    out = false;
  } else {
    util::error::exit(
      "[FATAL ERROR] Invalid Boolean Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible boolean value.\n");
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

//__Parse Positive Real Key Value_______________________________________________________________
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

//__Parse Real Range Key Value__________________________________________________________________
void _parse_key_value_real_range(const std::string& key,
                                 const std::string& value,
                                 real_range& out) {
  try {
    std::vector<std::string> values;
    util::string::split(value, values, ",;");
    if (values.size() != 2)
      throw std::invalid_argument{"Size Mismatch."};

    out = real_range{static_cast<real>(std::stold(values[0])),
                     static_cast<real>(std::stold(values[1]))};

  } catch (...) {
    util::error::exit(
      "[FATAL ERROR] Invalid Real Number Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible to Real range (min, max).\n");
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
                               r4_point& out,
                               bool use_units) {
  std::vector<std::string> point;
  util::string::split(value, point, ",;");
  const auto size = point.size();
  const auto& end = point.cend();
  util::error::exit_when(size != 4 || (std::find(point.cbegin(), end, "") != end),
    "[FATAL ERROR] Invalid R4 Point Argument for \"", key, "\" in Tracking Script.\n"
    "              Expected 4 arguments \"T, X, Y, Z\" but received ", size, ".\n");
  try {
    out = r4_point{
      (use_units ? units::time   : 1.0L) * static_cast<real>(std::stold(point[0])),
      (use_units ? units::length : 1.0L) * static_cast<real>(std::stold(point[1])),
      (use_units ? units::length : 1.0L) * static_cast<real>(std::stold(point[2])),
      (use_units ? units::length : 1.0L) * static_cast<real>(std::stold(point[3])) };
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
  tracking_options out{};

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
        // TODO: implement ...
      } else if (key == "data-directory") {
        _parse_key_value_file_path(key, value, out.data_directory);
      } else if (key == "data-file-extension") {
        out.data_file_extension = value;
      } else if (key == "data-position-keys") {
        _parse_key_value_data_keys(key, value,
          out.data_t_key, out.data_x_key, out.data_y_key, out.data_z_key);
      } else if (key == "data-position-error-keys") {
        _parse_key_value_data_keys(key, value,
          out.data_dt_key, out.data_dx_key, out.data_dy_key, out.data_dz_key);
      } else if (key == "data-detector-key") {
        out.data_detector_key = value;
      } else if (key == "data-track-id-key") {
        out.data_track_id_key = value;
      } else if (key == "data-parent-id-key") {
        out.data_parent_id_key = value;
      } else if (key == "data-momentum-keys") {
        _parse_key_value_data_keys(key, value,
          out.data_e_key, out.data_px_key, out.data_py_key, out.data_pz_key);
      } else if (key == "geometry-default-time-error") {
        _parse_key_value_positive_real(key, value, out.default_time_error);
        out.default_time_error *= units::time;
      } else if (key == "time-smearing") {
        _parse_key_value_boolean(key, value, out.time_smearing);
      } else if (key == "simulated-efficiency") {
        _parse_key_value_positive_real(key, value, out.simulated_efficiency);
      } else if (key == "simulated-noise-rate") {
        _parse_key_value_positive_real(key, value, out.simulated_noise_rate);
      } else if (key == "event-time-window") {
        _parse_key_value_real_range(key, value, out.event_time_window);
      } else if (key == "layer-axis") {
        _parse_key_value_r3_coordinate(key, value, out.layer_axis);
      } else if (key == "layer-depth") {
        _parse_key_value_positive_real(key, value, out.layer_depth);
        out.layer_depth *= units::length;
      } else if (key == "line-width") {
        _parse_key_value_positive_real(key, value, out.line_width);
        out.line_width *= units::length;
      } else if (key == "seed-size") {
        _parse_key_value_size_type(key, value, out.seed_size);
      } else if (key == "event-density-limit") {
        _parse_key_value_positive_real(key, value, out.event_density_limit);
      } else if (key == "event-overload-limit") {
        _parse_key_value_positive_real(key, value, out.event_overload_limit);
      } else if (key == "track-density-limit") {
        _parse_key_value_positive_real(key, value, out.track_density_limit);
      } else if (key == "statistics-directory") {
        _parse_key_value_file_path(key, value, out.statistics_directory);
      } else if (key == "statistics-file-prefix") {
        // FIXME: add checking parser
        out.statistics_file_prefix = value;
      } else if (key == "statistics-file-extension") {
        // FIXME: add checking parser
        out.statistics_file_extension = value;
      } else if (key == "verbose-output") {
        _parse_key_value_boolean(key, value, out.verbose_output);
      } else if (key == "draw-events") {
        _parse_key_value_boolean(key, value, out.draw_events);
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
  option help_opt    ('h', "help",      "MATHUSLA Tracking Algorithm", option::no_arguments);
  option verbose_opt ('v', "",          "Verbose Output",              option::no_arguments);
  option quiet_opt   ('q', "",          "Quiet Output",                option::no_arguments);
  option event_opt   (0, "draw-events", "Draw Events",                 option::no_arguments);

  // TODO: remove
  option geo_opt     ('g', "geometry",  "Geometry Import",             option::required_arguments);
  option map_opt     ('m', "map",       "Detector Map",                option::required_arguments);
  option data_opt    ('d', "data",      "ROOT Data Directory",         option::required_arguments);

  option script_opt  ('s', "script",    "Tracking Script",             option::required_arguments);

  util::cli::parse(argv, {&help_opt, &verbose_opt, &quiet_opt, &event_opt, &geo_opt, &data_opt, &map_opt, &script_opt});

   // TODO: remove
  util::error::exit_when((geo_opt.count && !data_opt.count)
                      || (data_opt.count && !geo_opt.count)
                      || !(script_opt.count || geo_opt.count || data_opt.count)
                      || argc == 1,
    "[FATAL ERROR] Insufficient Arguments:\n",
    "              Must include arguments for geometry ",
                  "and ROOT directory (-gd) or arguments for a tracking script (-s).\n");

  if (script_opt.count)
    _exit_on_missing_path(script_opt.argument, "Tracking Script");

  reader::tracking_options out;

  if (script_opt.count) {
    out = script::read(script_opt.argument);
    geo_opt.count += !out.geometry_file.empty();
    map_opt.count += !out.geometry_map_file.empty();
    data_opt.count += !out.data_directory.empty();
  } else {
     // TODO: remove
    out.geometry_file = geo_opt.count ? geo_opt.argument : "";
    out.geometry_map_file = map_opt.count ? map_opt.argument : "";
    out.data_directory = data_opt.count ? data_opt.argument : "";
  }

   // TODO: remove
  if (geo_opt.count) _exit_on_missing_path(out.geometry_file, "Geometry File");
  if (map_opt.count) _exit_on_missing_path(out.geometry_map_file, "Geometry Map");
  if (data_opt.count) _exit_on_missing_path(out.data_directory, "ROOT Directory");

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

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

#include <unordered_map>

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
      // FIXME: clean up and add range check
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
const geometry::detector_map import_detector_map(const script::path_type& path) {
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
const geometry::time_resolution_map import_time_resolution_map(const script::path_type& path) {
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

//__Search and Collect File Paths_______________________________________________________________
const script::path_vector search_directory(const script::path_type& path,
                                           const std::string& ext) {
  helper::init();
  return helper::search_directory(path, ext);
}
//----------------------------------------------------------------------------------------------

//__Search and Collect File Paths_______________________________________________________________
const std::vector<script::path_vector> search_directories(const script::path_vector& paths,
                                                          const std::string& ext) {
  std::vector<script::path_vector> out;
  out.reserve(paths.size());
  util::algorithm::back_insert_transform(paths, out,
    [&](const auto& path) { return search_directory(path, ext); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Transpose Search Directories________________________________________________________________
const std::vector<script::path_vector> transpose_search_directories(const script::path_vector& paths,
                                                                    const std::string& ext) {
  const auto directories = search_directories(paths, ext);
  std::vector<script::path_vector> transposed;
  std::size_t path_counter{};
  bool added_events = false;
  do {
    added_events = false;
    transposed.emplace_back();
    for (const auto& inner_paths : directories) {
      if (inner_paths.size() <= path_counter) {
        transposed.back().push_back("");
      } else {
        transposed.back().push_back(inner_paths[path_counter]);
        added_events = true;
      }
    }
    ++path_counter;
  } while (added_events);
  transposed.pop_back();
  transposed.shrink_to_fit();
  return transposed;
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Fix Function Overload Selection Using Lambda________________________________________________
#define _LAMBDA(f) [&](auto&&... args) { return f(std::forward<decltype(args)>(args)...); }
//----------------------------------------------------------------------------------------------

//__Append Events to EventVector________________________________________________________________
template<class EventVector>
void _merge_events(EventVector& front,
                   const EventVector& back) {
  const auto front_size = front.size();
  const auto back_size = back.size();
  const auto min_size = std::min(front_size, back_size);

  std::size_t index{};
  for (; index < min_size; ++index) {
    auto& event = front[index];
    auto& add_in = back[index];
    event.reserve(event.size() + add_in.size());
    event.insert(event.cend(), add_in.cbegin(), add_in.cend());
  }
  if (front_size < back_size) {
    for (; index < back_size; ++index)
      front.push_back(back[index]);
  }
}
//----------------------------------------------------------------------------------------------

//__Time Shift Events By Offset_________________________________________________________________
template<class EventVector>
void _time_shift_events(EventVector& events,
                        const real timing_offset) {
  for (auto& event : events)
    for (auto& hit : event)
      hit.t += timing_offset;
}
//----------------------------------------------------------------------------------------------

//__Track ID Helper Types_______________________________________________________________________
using _track_id_vector = std::vector<std::size_t>;
using _track_id_map = std::unordered_map<std::size_t, std::size_t>;
//----------------------------------------------------------------------------------------------

//__Find Track IDs in Map without Matching Path Index___________________________________________
const _track_id_vector _find_ids_without_matching_path_index(const _track_id_map& track_ids,
                                                             const std::size_t path_index) {
  _track_id_vector out{};
  for (const auto& entry : track_ids) {
    if (entry.second != path_index)
      out.push_back(entry.first);
  }
  out.erase(std::unique(out.begin(), out.end()), out.cend());
  return util::algorithm::sort_range(out);
}
//----------------------------------------------------------------------------------------------

//__Reassign Track ID for Conflicting Track in Event from Parent Information____________________
bool _change_track_id_from_parent(std::size_t& hit_id,
                                  _track_id_vector& changed_ids,
                                  _track_id_vector& resultant_ids) {
  const auto search = util::algorithm::range_binary_find_first(changed_ids, hit_id);
  if (search != changed_ids.cend()) {
    hit_id = resultant_ids[search - changed_ids.cbegin()];
    return true;
  }
  return false;
}
//----------------------------------------------------------------------------------------------

//__Reassign Track IDs for Conflicting Tracks in Event__________________________________________
template<class MCEvent>
void _change_track_ids(const std::size_t path_index,
                       MCEvent& true_event,
                       _track_id_map& track_ids,
                       _track_id_vector& changed_ids,
                       _track_id_vector& resultant_ids) {
  const auto size = true_event.size();
  const auto reserved_ids = _find_ids_without_matching_path_index(track_ids, path_index);

  for (std::size_t i{}; i < size; ++i) {
    auto& hit_id = true_event[i].track_id;
    if (_change_track_id_from_parent(hit_id, changed_ids, resultant_ids)) {
      continue;
    } else if (!util::algorithm::binary_search_range(reserved_ids, hit_id)) {
      changed_ids.push_back(hit_id);
      resultant_ids.push_back(hit_id);
    }

    auto attempt = hit_id + 1UL;
    while (util::algorithm::binary_search_range(reserved_ids, attempt)
        || util::algorithm::binary_search_range(changed_ids, attempt)
        || util::algorithm::binary_search_range(resultant_ids, attempt)) {
      ++attempt;
    }
    changed_ids.push_back(hit_id);
    resultant_ids.push_back(attempt);
    hit_id = attempt;
    track_ids.emplace(hit_id, path_index);
  }
}
//----------------------------------------------------------------------------------------------

//__Reassign Track IDs for Conflicting Tracks___________________________________________________
template<class MCEventVector>
void _reassign_track_ids(const std::size_t path_index,
                         MCEventVector& true_events,
                         _track_id_map& track_ids) {
  _track_id_vector changed_ids, resultant_ids;
  for (auto& event : true_events) {
    for (auto& hit : event) {
      const auto search = track_ids.find(hit.track_id);
      if (search != track_ids.cend()) {
        if (search->second != path_index) {
          if (changed_ids.empty() || !_change_track_id_from_parent(hit.track_id, changed_ids, resultant_ids)) {
            _change_track_ids(path_index, event, track_ids, changed_ids, resultant_ids);
            break;
          }
        }
        continue;
      }
      track_ids.emplace(hit.track_id, path_index);
    }
  }
}
//----------------------------------------------------------------------------------------------

//__Parallel Import of Events___________________________________________________________________
template<class EventVector, class Importer, class ...Args>
const EventVector _merge_import_events(const script::path_vector& paths,
                                       const real_vector& timing_offsets,
                                       Importer import,
                                       Args&& ...args) {
  EventVector out{};
  const auto path_count = paths.size();
  if (path_count == 0UL || path_count != timing_offsets.size())
    return out;

  for (std::size_t i{}; i < path_count; ++i) {
    if (paths[i].empty())
      continue;
    auto events = import(paths[i], std::forward<Args>(args)...);
    _time_shift_events(events, timing_offsets[i]);
    _merge_events(out, events);
  }

  return out;
}
//----------------------------------------------------------------------------------------------

//__Parallel Import of Events and Monte Carlo___________________________________________________
template<class EventVectorBundle, class Importer, class ...Args>
const EventVectorBundle _merge_import_events_and_mc(const script::path_vector& paths,
                                                    const real_vector& timing_offsets,
                                                    Importer import,
                                                    Args&& ...args) {
  EventVectorBundle out{{}, {}};
  const auto path_count = paths.size();
  if (path_count == 0UL || path_count != timing_offsets.size())
    return out;

  _track_id_map track_ids;
  for (std::size_t i{}; i < path_count; ++i) {
    if (paths[i].empty())
      continue;
    auto bundle = import(paths[i], std::forward<Args>(args)...);
    auto& events = bundle.events;
    _time_shift_events(events, timing_offsets[i]);
    _merge_events(out.events, events);
    auto& true_events = bundle.true_events;
    _time_shift_events(true_events, timing_offsets[i]);
    _reassign_track_ids(i, true_events, track_ids);
    _merge_events(out.true_events, true_events);
  }

  return out;
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Import Events from ROOT File________________________________________________________________
const analysis::event_vector import_events(const script::path_type& path,
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
const analysis::event_vector import_events(const script::path_type& path,
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
const analysis::event_vector import_events(const script::path_type& path,
                                           const script::tracking_options& options,
                                           const geometry::detector_map& map) {
  return import_events(path, options.data_t_key, options.data_detector_key, map);
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
const analysis::event_vector import_events(const script::path_type& path,
                                           const script::tracking_options& options,
                                           const ImportMode mode) {
  return mode == ImportMode::Detector
    ? import_events(path, options, import_detector_map(options.geometry_map_file))
    : import_events(path, options.data_t_key, options.data_x_key, options.data_y_key, options.data_z_key);
}
//----------------------------------------------------------------------------------------------

//__ROOT Event Parallel Import__________________________________________________________________
const analysis::event_vector import_events(const script::path_vector& paths,
                                           const real_vector& timing_offsets,
                                           const std::string& t_key,
                                           const std::string& x_key,
                                           const std::string& y_key,
                                           const std::string& z_key) {
  return _merge_import_events<analysis::event_vector>(paths, timing_offsets,
    _LAMBDA(import_events), t_key, x_key, y_key, z_key);
}
//----------------------------------------------------------------------------------------------

//__ROOT Event Parallel Import__________________________________________________________________
const analysis::event_vector import_events(const script::path_vector& paths,
                                           const real_vector& timing_offsets,
                                           const std::string& t_key,
                                           const std::string& detector_key,
                                           const geometry::detector_map& map) {
  return _merge_import_events<analysis::event_vector>(paths, timing_offsets,
    _LAMBDA(import_events), t_key, detector_key, map);
}
//----------------------------------------------------------------------------------------------

//__ROOT Event Parallel Import__________________________________________________________________
const analysis::event_vector import_events(const script::path_vector& paths,
                                           const script::tracking_options& options,
                                           const geometry::detector_map& map) {
  return import_events(paths, options.data_timing_offsets, options.data_t_key, options.data_detector_key, map);
}
//----------------------------------------------------------------------------------------------

//__ROOT Event Parallel Import__________________________________________________________________
const analysis::event_vector import_events(const script::path_vector& paths,
                                           const script::tracking_options& options,
                                           const ImportMode mode) {
  return mode == ImportMode::Detector
    ? import_events(paths, options, import_detector_map(options.geometry_map_file))
    : import_events(paths, options.data_timing_offsets, options.data_t_key, options.data_x_key, options.data_y_key, options.data_z_key);
}
//----------------------------------------------------------------------------------------------

//__Import Full Events from ROOT File___________________________________________________________
const analysis::full_event_vector import_full_events(const script::path_type& path,
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
const analysis::full_event_vector import_full_events(const script::path_type& path,
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
const analysis::full_event_vector import_full_events(const script::path_type& path,
                                                     const script::tracking_options& options,
                                                     const geometry::detector_map& map) {
  return import_full_events(path, options.data_t_key, options.data_dt_key, options.data_detector_key, map);
}
//----------------------------------------------------------------------------------------------

//__Import Full Events from ROOT File___________________________________________________________
const analysis::full_event_vector import_full_events(const script::path_type& path,
                                                     const script::tracking_options& options,
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

//__ROOT Full Event Parallel Import_____________________________________________________________
const analysis::full_event_vector import_full_events(const script::path_vector& paths,
                                                     const real_vector& timing_offsets,
                                                     const std::string& t_key,
                                                     const std::string& x_key,
                                                     const std::string& y_key,
                                                     const std::string& z_key,
                                                     const std::string& dt_key,
                                                     const std::string& dx_key,
                                                     const std::string& dy_key,
                                                     const std::string& dz_key) {
  return _merge_import_events<analysis::full_event_vector>(paths, timing_offsets,
    _LAMBDA(import_full_events), t_key, x_key, y_key, z_key, dt_key, dx_key, dy_key, dz_key);
}
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Parallel Import_____________________________________________________________
const analysis::full_event_vector import_full_events(const script::path_vector& paths,
                                                     const real_vector& timing_offsets,
                                                     const std::string& t_key,
                                                     const std::string& dt_key,
                                                     const std::string& detector_key,
                                                     const geometry::detector_map& map) {
  return _merge_import_events<analysis::full_event_vector>(paths, timing_offsets,
    _LAMBDA(import_full_events), t_key, dt_key, detector_key, map);
}
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Parallel Import_____________________________________________________________
const analysis::full_event_vector import_full_events(const script::path_vector& paths,
                                                     const script::tracking_options& options,
                                                     const geometry::detector_map& map) {
  return import_full_events(paths, options.data_timing_offsets, options.data_t_key, options.data_dt_key, options.data_detector_key, map);
}
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Parallel Import_____________________________________________________________
const analysis::full_event_vector import_full_events(const script::path_vector& paths,
                                                     const script::tracking_options& options,
                                                     const ImportMode mode) {
  return mode == ImportMode::Detector
    ? import_full_events(paths, options, import_detector_map(options.geometry_map_file))
    : import_full_events(paths,
        options.data_timing_offsets,
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
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_type& path,
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
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_type& path,
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
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_type& path,
                                                               const script::tracking_options& options,
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
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_type& path,
                                                               const script::tracking_options& options,
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

//__ROOT Event Parallel Import with Monte-Carlo Tracks__________________________________________
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_vector& paths,
                                                               const real_vector& timing_offsets,
                                                               const std::string& track_key,
                                                               const std::string& t_key,
                                                               const std::string& x_key,
                                                               const std::string& y_key,
                                                               const std::string& z_key) {
  return _merge_import_events_and_mc<analysis::mc::event_vector_bundle>(paths, timing_offsets,
    _LAMBDA(import_event_mc_bundle), track_key, t_key, x_key, y_key, z_key);
}
//----------------------------------------------------------------------------------------------

//__ROOT Event Parallel Import with Monte-Carlo Tracks__________________________________________
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_vector& paths,
                                                               const real_vector& timing_offsets,
                                                               const std::string& track_key,
                                                               const std::string& t_key,
                                                               const std::string& x_key,
                                                               const std::string& y_key,
                                                               const std::string& z_key,
                                                               const std::string& detector_key,
                                                               const geometry::detector_map& map) {
  return _merge_import_events_and_mc<analysis::mc::event_vector_bundle>(paths, timing_offsets,
    _LAMBDA(import_event_mc_bundle), track_key, t_key, x_key, y_key, z_key, detector_key, map);
}
//----------------------------------------------------------------------------------------------

//__ROOT Event Parallel Import with Monte-Carlo Tracks__________________________________________
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_vector& paths,
                                                               const script::tracking_options& options,
                                                               const geometry::detector_map& map) {
  return import_event_mc_bundle(paths,
    options.data_timing_offsets,
    options.data_track_id_key,
    options.data_t_key,
    options.data_x_key,
    options.data_y_key,
    options.data_z_key,
    options.data_detector_key,
    map);
}
//----------------------------------------------------------------------------------------------

//__ROOT Event Parallel Import with Monte-Carlo Tracks__________________________________________
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_vector& paths,
                                                               const script::tracking_options& options,
                                                               const ImportMode mode) {
  return mode == ImportMode::Detector
    ? import_event_mc_bundle(paths, options, import_detector_map(options.geometry_map_file))
    : import_event_mc_bundle(paths,
        options.data_timing_offsets,
        options.data_track_id_key,
        options.data_t_key,
        options.data_x_key,
        options.data_y_key,
        options.data_z_key);
}
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Import with Monte-Carlo Tracks______________________________________________
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_type& path,
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
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_type& path,
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
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_type& path,
                                                                         const script::tracking_options& options,
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
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_type& path,
                                                                         const script::tracking_options& options,
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

//__ROOT Full Event Parallel Import with Monte-Carlo Tracks_____________________________________
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_vector& paths,
                                                                         const real_vector& timing_offsets,
                                                                         const std::string& track_key,
                                                                         const std::string& t_key,
                                                                         const std::string& x_key,
                                                                         const std::string& y_key,
                                                                         const std::string& z_key,
                                                                         const std::string& dt_key,
                                                                         const std::string& dx_key,
                                                                         const std::string& dy_key,
                                                                         const std::string& dz_key) {
  return _merge_import_events_and_mc<analysis::mc::full_event_vector_bundle>(paths, timing_offsets,
    _LAMBDA(import_full_event_mc_bundle), track_key, t_key, x_key, y_key, z_key, dt_key, dx_key, dy_key, dz_key);
}
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Parallel Import with Monte-Carlo Tracks_____________________________________
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_vector& paths,
                                                                         const real_vector& timing_offsets,
                                                                         const std::string& track_key,
                                                                         const std::string& t_key,
                                                                         const std::string& x_key,
                                                                         const std::string& y_key,
                                                                         const std::string& z_key,
                                                                         const std::string& dt_key,
                                                                         const std::string& detector_key,
                                                                         const geometry::detector_map& map) {
  return _merge_import_events_and_mc<analysis::mc::full_event_vector_bundle>(paths, timing_offsets,
    _LAMBDA(import_full_event_mc_bundle), track_key, t_key, x_key, y_key, z_key, dt_key, detector_key, map);
}
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Parallel Import with Monte-Carlo Tracks_____________________________________
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_vector& paths,
                                                                         const script::tracking_options& options,
                                                                         const geometry::detector_map& map) {
  return import_full_event_mc_bundle(paths,
    options.data_timing_offsets,
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

//__ROOT Full Event Parallel Import with Monte-Carlo Tracks_____________________________________
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_vector& paths,
                                                                         const script::tracking_options& options,
                                                                         const ImportMode mode) {
  return mode == ImportMode::Detector
    ? import_full_event_mc_bundle(paths, options, import_detector_map(options.geometry_map_file))
    : import_full_event_mc_bundle(paths,
        options.data_timing_offsets,
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

//__Merge Save Files____________________________________________________________________________
void merge_save(const script::path_type& home_path,
                const script::path_vector& paths,
                const std::vector<std::string>& prefixes) {
  const auto path_count = paths.size();
  if (path_count == 0UL || (path_count != prefixes.size() && !prefixes.empty()))
    return;

  auto prefixes_copy = prefixes;
  if (prefixes_copy.empty()) {
    prefixes_copy.reserve(path_count);
    prefixes_copy.resize(path_count, "");
  }

  helper::while_open(home_path, "UPDATE", [&](auto file) {
    for (std::size_t i{}; i < path_count; ++i) {
      helper::traverse_keys(paths[i], "READ", [&](auto inner_file, auto key) {
        auto clazz = helper::get_key_class(key);
        if (!clazz) return;
        if (clazz->InheritsFrom("TDirectory")) {
          /* TODO: implement
          source->cd(key->GetName());
          TDirectory *subdir = gDirectory;
          adir->cd();
          CopyDir(subdir);
          adir->cd();
          */
        } else if (clazz->InheritsFrom("TTree")) {
          inner_file->cd();
          auto tree = helper::tree::unchecked_get_tree(inner_file, key);
          file->cd();
          auto new_tree = tree->CloneTree();
          new_tree->SetName((prefixes_copy[i] + new_tree->GetName()).c_str());
          new_tree->Write();
        } else {
          inner_file->cd();
          auto obj = key->ReadObj();
          if (clazz->InheritsFrom("TNamed")) {
            auto named = dynamic_cast<TNamed*>(obj);
            named->SetName((prefixes_copy[i] + named->GetName()).c_str());
            file->cd();
            named->Write();
          } else {
            file->cd();
            obj->Write();
          }
          delete obj;
        }
      });
    }
  });
}
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

} /* namespace reader */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

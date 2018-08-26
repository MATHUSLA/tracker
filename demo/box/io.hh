/*
 * demo/box/io.hh
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

#ifndef TRACKER__BOX__IO_HH
#define TRACKER__BOX__IO_HH
#pragma once

#include <iostream>

#include <tracker/analysis/monte_carlo.hh>
#include <tracker/analysis/track.hh>
#include <tracker/analysis/vertex.hh>
#include <tracker/plot.hh>
#include <tracker/script.hh>

#include <tracker/util/time.hh>

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace plot     = MATHUSLA::TRACKER::plot;
namespace script   = MATHUSLA::TRACKER::script;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

namespace box { namespace io { /////////////////////////////////////////////////////////////////

//__Extension Parser for Tracking Script________________________________________________________
struct extension_parser {
  std::size_t layer_count;
  type::real scintillator_x_width;
  type::real scintillator_y_width;
  type::real scintillator_height;
  type::real layer_spacing;
  type::real x_displacement;
  type::real y_displacement;
  type::real x_edge_length;
  type::real y_edge_length;

  extension_parser();
  void operator()(const std::string& key,
                  const std::string& value,
                  script::tracking_options& options);
};
//----------------------------------------------------------------------------------------------

//__Draw Main Detector To Canvas________________________________________________________________
void draw_detector(plot::canvas& canvas,
                   std::size_t layer_count);
//----------------------------------------------------------------------------------------------

//__Add Track and Intersecting Geometry to Canvas_______________________________________________
void draw_track(plot::canvas& canvas,
                const analysis::track& track);
//----------------------------------------------------------------------------------------------

//__Draw MC Tracks to Canvas____________________________________________________________________
void draw_mc_tracks(plot::canvas& canvas,
                    const analysis::mc::track_vector& tracks);
//----------------------------------------------------------------------------------------------

//__Add Track and Intersecting Geometry to Canvas_______________________________________________
void draw_vertex_and_guess(plot::canvas& canvas,
                           const analysis::vertex& vertex);
//----------------------------------------------------------------------------------------------

//__Show and Add Tracks to Statistics___________________________________________________________
void save_tracks(const analysis::track_vector& tracks,
                 plot::canvas& canvas,
                 analysis::track::tree& tree,
                 const script::tracking_options& options);
//----------------------------------------------------------------------------------------------

//__Show and Add Vertices to Statistics_________________________________________________________
void save_vertices(const analysis::vertex_vector& vertices,
                   plot::canvas& canvas,
                   analysis::vertex::tree& tree,
                   const script::tracking_options& options);
//----------------------------------------------------------------------------------------------

//__Get Statistics File Path and Add Directories________________________________________________
const std::string add_statistics_path(const script::tracking_options& options);
//----------------------------------------------------------------------------------------------

//__Calculate Value Tags for Paths______________________________________________________________
const plot::value_tag_vector data_paths_value_tags(const script::path_vector& paths,
                                                   const type::real_vector& timing_offsets,
                                                   const std::size_t starting_index=0UL);
//----------------------------------------------------------------------------------------------

//__Save and Merge Simulation and Tracking Files________________________________________________
void save_files(const script::path_type& save_path,
                analysis::track::tree& track_tree,
                analysis::vertex::tree& vertex_tree,
                const script::path_vector& paths,
                const bool merge);
//----------------------------------------------------------------------------------------------

//__Print Bar___________________________________________________________________________________
inline void print_bar(const size_t count=99,
                      const char b='=') {
  std::cout << "\n" << std::string(count, b) << "\n\n";
}
//----------------------------------------------------------------------------------------------

//__Print Directories For Tracking______________________________________________________________
inline void print_tracking_directories(const script::path_vector& directories) {
  if (directories.size() <= 1UL) {
    std::cout << "Begin Tracking in " << directories.front() << ":\n\n";
  } else {
    std::cout << "Begin Parallel Tracking in: \n";
    for (const auto& path : directories)
      std::cout << "  - " << path << "\n";
    std::cout << "\n";
  }
}
//----------------------------------------------------------------------------------------------

//__Print Paths For Tracking____________________________________________________________________
inline void print_tracking_paths(const script::path_vector& paths) {
  if (paths.size() <= 1UL) {
    std::cout << "Read Path: " << paths.front() << "\n";
  } else {
    std::cout << "Read Paths: \n";
    for (const auto& path : paths)
      std::cout << "  - " << path << "\n";
  }
}
//----------------------------------------------------------------------------------------------

//__Print Event Summary_________________________________________________________________________
inline void print_event_summary(const std::size_t event_counter,
                                const std::size_t event_size,
                                const type::real compression_size,
                                const type::real event_density) {
  std::cout << "Event " << event_counter << " with " << event_size << " hits.\n"
            << "  Compression Size: " << compression_size << "\n"
            << "  Event Density: "    << event_density * 100.0L << " %\n";
}
//----------------------------------------------------------------------------------------------

//__Print Tracking Summary______________________________________________________________________
inline void print_tracking_summary(const analysis::full_event& event,
                                   const analysis::track_vector& tracks) {
  const auto size = std::accumulate(tracks.cbegin(), tracks.cend(), 0UL,
    [](const auto sum, const auto& track) { return sum + track.fit_converged(); });
  std::cout << "  Track Count: "   << size << "\n"
            << "  Track Density: " << size / static_cast<type::real>(event.size())
                                        * 100.0L << " %\n";
}
//----------------------------------------------------------------------------------------------

} } /* namespace box::io */ ////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* TRACKER__BOX__IO_HH */

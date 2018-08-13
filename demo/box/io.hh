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
#include <tracker/reader.hh>

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace plot     = MATHUSLA::TRACKER::plot;
namespace reader   = MATHUSLA::TRACKER::reader;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

namespace box { ////////////////////////////////////////////////////////////////////////////////

//__Extension Parser for Tracking Script________________________________________________________
struct extension_parser {
  std::size_t layer_count;
  type::real scintillator_x_width, scintillator_y_width;

  extension_parser();
  void operator()(const std::string& key,
                  const std::string& value,
                  reader::tracking_options& options);
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
                 const reader::tracking_options& options);
//----------------------------------------------------------------------------------------------

//__Show and Add Vertices to Statistics_________________________________________________________
void save_vertices(const analysis::vertex_vector& vertices,
                   plot::canvas& canvas,
                   analysis::vertex::tree& tree,
                   const reader::tracking_options& options);
//----------------------------------------------------------------------------------------------

//__Print Bar___________________________________________________________________________________
inline void print_bar(const size_t count=99,
                      const char b='=') {
  std::cout << "\n" << std::string(count, b) << "\n\n";
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
  std::cout << "  Track Count: "   << tracks.size() << "\n"
            << "  Track Density: " << tracks.size() / static_cast<type::real>(event.size())
                                        * 100.0L << " %\n";
}
//----------------------------------------------------------------------------------------------

} /* namespace box */ //////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* TRACKER__BOX__IO_HH */

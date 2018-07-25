/*
 * demo/box/logging.cc
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

#include "logging.hh"

#include <tracker/core/units.hh>

#include "geometry.hh"

namespace MATHUSLA {

namespace box { ////////////////////////////////////////////////////////////////////////////////

//__Draw Main Detector To Canvas________________________________________________________________
void draw_detector(plot::canvas& canvas) {
  for (const auto& volume : tracker_geometry::full_structure_except({"World",
                                                                     "Box",
                                                                     "SteelPlate",
                                                                     "ModifiedSandstone",
                                                                     "Marl",
                                                                     "Mix",
                                                                     "Earth"})) {
    const auto limits = tracker_geometry::limits_of(volume);
    canvas.add_box(limits.center,
                   limits.max.x - limits.min.x,
                   limits.max.y - limits.min.y,
                   limits.max.z - limits.min.z,
                   2.5,
                   plot::color::MAGENTA);
  }
}
//----------------------------------------------------------------------------------------------

//__Add Track and Intersecting Geometry to Canvas_______________________________________________
void draw_track(plot::canvas& canvas,
                const analysis::track& track) {
  const auto& full_event = track.full_event();
  const auto size = full_event.size();
  if (size == 0UL)
    return;
  uint_fast8_t brightness = 0, step = 230 / size;
  for (const auto& point : full_event) {
    const auto center = type::reduce_to_r3(point);
    canvas.add_box(center,
                   point.width.x, point.width.y, point.width.z,
                   2.5,
                   {brightness, brightness, brightness});
    brightness += step;
  }
  track.draw(canvas, 2, plot::color::RED);
}
//----------------------------------------------------------------------------------------------

//__Draw MC Tracks to Canvas____________________________________________________________________
void draw_mc_tracks(plot::canvas& canvas,
                    const analysis::mc::track_vector& tracks) {
  for (const auto& track : tracks) {
    const auto hits = track.hits;
    const auto size = hits.size();
    for (const auto& point : hits) {
      const auto center = type::reduce_to_r3(point);
      canvas.add_point(center, 0.5, plot::color::BLUE);
      const auto box = geometry::limits_of_volume(center);
      canvas.add_box(box.center,
                     box.max.x - box.min.x,
                     box.max.y - box.min.y,
                     box.max.z - box.min.z,
                     5,
                     plot::color::BLUE);
    }
    for (std::size_t i = 0; i < size - 1; ++i)
      canvas.add_line(type::reduce_to_r3(hits[i]), type::reduce_to_r3(hits[i+1]), 1, plot::color::BLUE);
  }
}
//----------------------------------------------------------------------------------------------

//__Add Track and Intersecting Geometry to Canvas_______________________________________________
void draw_vertex_and_guess(plot::canvas& canvas,
                           const analysis::vertex& vertex) {
  vertex.draw(canvas, 1.2, plot::color::GREEN);

  if (vertex.fit_converged()) {
    vertex.draw_guess(canvas, 1.2, plot::color::RED);
    const auto point = vertex.point();
    for (const auto& track : vertex.tracks())
      canvas.add_line(track.at_t(point.t), point, 2, plot::color::GREEN);
  }
}
//----------------------------------------------------------------------------------------------

//__Show and Add Tracks to Statistics___________________________________________________________
void save_tracks(const analysis::track_vector& tracks,
                 plot::canvas& canvas,
                 analysis::track::tree& tree,
                 const reader::tracking_options& options) {
  for (const auto& track : tracks) {
    if (options.verbose_output)
      std::cout << track << "\n";
    if (options.draw_events)
      draw_track(canvas, track);
  }
  tree.fill(tracks);
}
//----------------------------------------------------------------------------------------------

//__Show and Add Vertices to Statistics_________________________________________________________
void save_vertices(const analysis::vertex_vector& vertices,
                   plot::canvas& canvas,
                   analysis::vertex::tree& tree,
                   const reader::tracking_options& options) {
  for (const auto& vertex : vertices) {
    if (options.verbose_output)
      std::cout << vertex << "\n";
    if (options.draw_events)
      draw_vertex_and_guess(canvas, vertex);
  }
  tree.fill(vertices);
}
//----------------------------------------------------------------------------------------------

} /* namespace box */ //////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

/*
 * demo/prototype/logging.cc
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

//__Add Detector Centers to Canvas______________________________________________________________
void draw_detector_centers(plot::canvas& canvas) {
  for (const auto& name : prototype_geometry())
    canvas.add_point(geometry::limits_of(name).center, 0.25, plot::color::MAGENTA);
}
//----------------------------------------------------------------------------------------------

//__Add Track and Intersecting Geometry to Canvas_______________________________________________
void draw_track(plot::canvas& canvas,
                const analysis::track& track) {
  const auto& full_event = track.full_event();
  uint_fast8_t brightness = 0, step = 230 / full_event.size();
  for (const auto& point : full_event) {
    const auto center = type::reduce_to_r3(point);
    canvas.add_box(center,
                   point.width.x, point.width.y, point.width.z,
                   2.5,
                   {brightness, brightness, brightness});
    brightness += step;
  }
  canvas.add_line(track.front(), track.back(), 1, plot::color::RED);
  for (std::size_t i = 0; i < full_event.size() - 1; ++i) {
    canvas.add_line(type::reduce_to_r3(full_event[i]), type::reduce_to_r3(full_event[i+1]), 1, plot::color::BLACK);
  }
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
                     1,
                     plot::color::BLUE);
    }
    for (std::size_t i = 0; i < size - 1; ++i) {
      canvas.add_line(type::reduce_to_r3(hits[i]), type::reduce_to_r3(hits[i+1]), 1, plot::color::BLUE);
    }
  }
}
//----------------------------------------------------------------------------------------------

//__Add Track and Intersecting Geometry to Canvas_______________________________________________
void draw_vertex_and_guess(plot::canvas& canvas,
                           const analysis::vertex& vertex) {
  const auto point = vertex.point();
  const auto point_error = vertex.point_error();
  canvas.add_point(point, 1.5, plot::color::GREEN);
  canvas.add_box(point, point_error.x, point_error.y, point_error.z, 2.5, plot::color::GREEN);

  const auto guess = vertex.guess_fit();
  canvas.add_point(type::r3_point{guess.x.value, guess.y.value, guess.z.value}, 1.5, plot::color::RED);

  for (const auto& track : vertex.tracks())
    canvas.add_line(track(point.z), track.back(), 1.5, plot::color::GREEN);
}
//----------------------------------------------------------------------------------------------

//__Generate Histograms for Prototype___________________________________________________________
plot::histogram_collection generate_histograms() {
  return plot::histogram_collection({
    {"track_chi_squared",   "Track Chi-Squared Distribution",    "chi^2/dof",    "Track Count",  200, 0, 10},
    {"vertex_chi_squared",  "Vertex Chi-Squared Distribution",   "chi^2/dof",    "Vertex Count", 200, 0, 10},
    {"beta",                "Track Beta Distribution",           "beta",         "Track Count",  200, 0,  2},
    {"beta_error",          "Track Beta Error Distribution",     "beta error",   "Track Count",  200, 0,  2},
    {"beta_with_cut",       "Track Beta Distribution With Cut",  "beta",         "Track Count",  200, 0,  2},
    {"track_count",         "Track Count Distribution",          "Track Count",  "Event Count",   50, 0, 50},
    {"vertex_count",        "Vertex Count Distribution",         "Vertex Count", "Event Count",  100, 0, 50},
    {"track_size",          "Track Size Distribution",           "Hit Count",    "Track Count",   50, 0, 50},
    {"non_track_hit_count", "Non-Track Hit Count Distribution",  "Hit Count",    "Event Count",  100, 0, 50},
    {"vertex_t_error",      "Vertex T Error Distribution",       "t error (" + units::time_string   + ")", "Vertex Count", 100, 0, 20},
    {"vertex_x_error",      "Vertex X Error Distribution",       "x error (" + units::length_string + ")", "Vertex Count", 100, 0, 100},
    {"vertex_y_error",      "Vertex Y Error Distribution",       "y error (" + units::length_string + ")", "Vertex Count", 100, 0, 100},
    {"vertex_z_error",      "Vertex Z Error Distribution",       "z error (" + units::length_string + ")", "Vertex Count", 100, 0, 100},
  });
}
//----------------------------------------------------------------------------------------------

//__Show and Add Tracks to Statistics___________________________________________________________
void save_tracks(const analysis::track_vector& tracks,
                 plot::canvas& canvas,
                 plot::histogram_collection& histograms,
                 bool verbose) {
  histograms["track_count"].insert(tracks.size());
  for (const auto& track : tracks) {
    histograms["track_chi_squared"].insert(track.chi_squared_per_dof());
    const auto beta = track.beta();
    const auto beta_error = track.beta_error();
    histograms["beta"].insert(beta);
    if (beta + 3.0L * beta_error <= 1)
      histograms["beta_with_cut"].insert(beta);
    histograms["beta_error"].insert(beta_error);
    histograms["track_size"].insert(track.size());
    if (verbose) {
      draw_track(canvas, track);
      std::cout << track << "\n";
    }
  }
}
//----------------------------------------------------------------------------------------------

//__Show and Add Vertex to Statistics___________________________________________________________
void save_vertex(const analysis::vertex& vertex,
                 plot::canvas& canvas,
                 plot::histogram_collection& histograms,
                 bool verbose) {
  const auto size = vertex.size();
  if (size <= 1)
    return;

  histograms["vertex_chi_squared"].insert(vertex.chi_squared_per_dof());
  histograms["vertex_t_error"].insert(vertex.t_error() / units::time);
  histograms["vertex_x_error"].insert(vertex.x_error() / units::length);
  histograms["vertex_y_error"].insert(vertex.y_error() / units::length);
  histograms["vertex_z_error"].insert(vertex.z_error() / units::length);
  if (verbose) {
    draw_vertex_and_guess(canvas, vertex);
    std::cout << vertex << "\n";
  }
}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

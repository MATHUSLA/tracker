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
  const auto size = full_event.size();
  if (size == 0)
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
  vertex.draw(canvas, 1.2, plot::color::GREEN);
  vertex.draw_guess(canvas, 1.2, plot::color::RED);

  if (vertex.fit_converged()) {
    const auto point = vertex.point();
    for (const auto& track : vertex.tracks()) {
      canvas.add_line(track.at_t(point.t), point, 2, plot::color::GREEN);
    }
  }

}
//----------------------------------------------------------------------------------------------

//__Track Plotting Keys for Prototype___________________________________________________________
const analysis::track::plotting_keys& track_plotting_keys() {
  static analysis::track::plotting_keys _keys{
    "track_t0",
    "track_x0",
    "track_y0",
    "track_z0",
    "track_vx",
    "track_vy",
    "track_vz",
    "track_t0_error",
    "track_x0_error",
    "track_y0_error",
    "track_z0_error",
    "track_vx_error",
    "track_vy_error",
    "track_vz_error",
    "track_chi_squared",
    "track_chi_squared_per_dof",
    "track_chi_squared_p_value",
    "track_size",
    "track_beta",
    "track_beta_error",
    "track_angle",
    "track_angle_error"
  };
  return _keys;
}
//----------------------------------------------------------------------------------------------

//__Vertex Plotting Keys for Prototype__________________________________________________________
const analysis::vertex::plotting_keys& vertex_plotting_keys() {
  static analysis::vertex::plotting_keys _keys{
    "vertex_t",
    "vertex_x",
    "vertex_y",
    "vertex_z",
    "vertex_t_error",
    "vertex_x_error",
    "vertex_y_error",
    "vertex_z_error",
    "vertex_distance",
    "vertex_distance_error",
    "vertex_chi_squared",
    "vertex_chi_squared_per_dof",
    "vertex_chi_squared_p_value",
    "vertex_size",
  };
  return _keys;
}
//----------------------------------------------------------------------------------------------

//__Generate Histograms for Prototype___________________________________________________________
plot::histogram_collection generate_histograms() {
  static const auto& track = track_plotting_keys();
  static const auto& vertex = vertex_plotting_keys();
  static const auto time_unit = "(" + units::time_string + ")";
  static const auto length_unit = "(" + units::length_string + ")";
  static const auto velocity_unit = "(" + units::velocity_string + ")";
  return plot::histogram_collection({
    {track.t0, "Track T0 Distribution", "t0 " + time_unit,     "Track Count", 200,  300, 400},
    {track.x0, "Track X0 Distribution", "x0 " + length_unit,   "Track Count", 400, -600, 600},
    {track.y0, "Track Y0 Distribution", "y0 " + length_unit,   "Track Count", 400, -600, 600},
    {track.z0, "Track Z0 Distribution", "z0 " + length_unit,   "Track Count", 600, -600,   0},
    {track.vx, "Track VX Distribution", "vx " + velocity_unit, "Track Count", 200,  -35,  35},
    {track.vy, "Track VY Distribution", "vy " + velocity_unit, "Track Count", 200,  -35,  35},
    {track.vz, "Track VZ Distribution", "vz " + velocity_unit, "Track Count", 200,  -35,  35},

    {track.t0_error, "Track T0 Error Distribution", "t0 error " + time_unit,     "Track Count", 200, 0, 50},
    {track.x0_error, "Track X0 Error Distribution", "x0 error " + length_unit,   "Track Count", 200, 0, 50},
    {track.y0_error, "Track Y0 Error Distribution", "y0 error " + length_unit,   "Track Count", 200, 0, 50},
    {track.z0_error, "Track Z0 Error Distribution", "z0 error " + length_unit,   "Track Count", 200, 0, 50},
    {track.vx_error, "Track VX Error Distribution", "vx error " + velocity_unit, "Track Count", 200, 0, 10},
    {track.vy_error, "Track VY Error Distribution", "vy error " + velocity_unit, "Track Count", 200, 0, 10},
    {track.vz_error, "Track VZ Error Distribution", "vz error " + velocity_unit, "Track Count", 200, 0, 10},

    {track.chi_squared_per_dof, "Track Chi-Squared Distribution",           "chi^2/dof",    "Track Count",  200,  0,    12   },
    {track.beta,                "Track Beta Distribution",                  "#beta",        "Track Count",  200,  0,     2   },
    {"track_beta_with_cut",     "Track Beta Distribution With 1#sigma Cut", "#beta",        "Track Count",  200,  0,     2   },
    {track.beta_error,          "Track Beta Error Distribution",            "#beta error",  "Track Count",  200,  0,     2   },
    {track.angle,               "Track Angular Distribution",               "#theta",       "Track Count",  200,  0,     6.28},
    {track.angle_error,         "Track Angular Error Distribution",         "#theta error", "Track Count",  200,  0,     1   },
    {track.size,                "Track Size Distribution",                  "Hit Count",    "Track Count",   40,  0,    40   },

    {vertex.t, "Vertex T Distribution", "t " + time_unit,   "Vertex Count", 100,  300, 400},
    {vertex.x, "Vertex X Distribution", "x " + length_unit, "Vertex Count", 300, -100, 100},
    {vertex.y, "Vertex Y Distribution", "y " + length_unit, "Vertex Count", 300, -100, 100},
    {vertex.z, "Vertex Z Distribution", "z " + length_unit, "Vertex Count", 300, -600, 100},

    {vertex.t_error, "Vertex T Error Distribution", "t error " + time_unit,   "Vertex Count", 100, 0, 20},
    {vertex.x_error, "Vertex X Error Distribution", "x error " + length_unit, "Vertex Count", 100, 0, 20},
    {vertex.y_error, "Vertex Y Error Distribution", "y error " + length_unit, "Vertex Count", 100, 0, 20},
    {vertex.z_error, "Vertex Z Error Distribution", "z error " + length_unit, "Vertex Count", 100, 0, 20},

    {vertex.distance,       "Vertex Distance Distribution",       "distance "       + length_unit, "Vertex Count", 300, 0, 600},
    {vertex.distance_error, "Vertex Distance Error Distribution", "distance error " + length_unit, "Vertex Count", 500, 0, 1000},

    {vertex.chi_squared_per_dof, "Vertex Chi-Squared Distribution", "chi^2/dof",  "Vertex Count",  200, 0, 12},
    {vertex.size,                "Tracks per Vertex Distribution",  "Track Count", "Vertex Count",  10, 0, 10},

    {"non_track_hit_count", "Non-Track Hit Count Distribution", "Hit Count",    "Event Count", 100, 0, 50},
    {"track_count",         "Track Count Distribution",         "Track Count",  "Event Count",  50, 0, 50} //,
    // {"vertex_count",        "Vertex Count Distribution",        "Vertex Count", "Event Count", 100, 0, 50}
  });
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

} /* namespace MATHUSLA */

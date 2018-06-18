/*
 * demo/prototype/prototype_tracking.cc
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

#include <tracker/analysis.hh>
#include <tracker/geometry.hh>
#include <tracker/plot.hh>
#include <tracker/reader.hh>
#include <tracker/track.hh>
#include <tracker/units.hh>

#include "geometry.hh"

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace geometry = MATHUSLA::TRACKER::geometry;
namespace plot     = MATHUSLA::TRACKER::plot;
namespace reader   = MATHUSLA::TRACKER::reader;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

//__Find Tracks for Prototype___________________________________________________________________
const analysis::track_vector find_tracks(const analysis::event& event,
                                         const reader::tracking_options& options) {
  analysis::full_event combined_rpc_hits, original_rpc_hits;
  const auto optimized_event = combine_rpc_hits(event, combined_rpc_hits, original_rpc_hits);
  const auto layers          = analysis::partition(optimized_event, options.layer_axis, options.layer_depth);
  const auto seeds           = analysis::seed(options.seed_size, layers, options.line_width);
  const auto tracking_vector = reset_seeds(analysis::join_all(seeds), combined_rpc_hits, original_rpc_hits);
  return analysis::fit_seeds(tracking_vector);
}
//----------------------------------------------------------------------------------------------

//__Add Track and Intersecting Geometry to Canvas_______________________________________________
void draw_track(plot::canvas& canvas,
                const analysis::track& track) {
  const auto& full_event = track.full_event();
  uint_fast8_t brightness = 0, step = 230 / full_event.size();
  for (const auto& point : full_event) {
    const auto center = type::reduce_to_r3(point);
    const plot::color color{brightness, brightness, brightness};
    canvas.add_box(center, point.width.x, point.width.y, point.width.z, 2.5, color);
    canvas.add_point(center, 0.3, color);
    brightness += step;
  }
  canvas.add_line(type::reduce_to_r3(full_event.front()), type::reduce_to_r3(full_event.back()));
  canvas.add_line(track.front(), track.back(), 1, plot::color::RED);
}
//----------------------------------------------------------------------------------------------

//__Add Detector Centers to Canvas______________________________________________________________
void draw_detector_centers(plot::canvas& canvas) {
  for (const auto& name : prototype_geometry())
    canvas.add_point(geometry::limits_of(name).center, 0.25, plot::color::MAGENTA);
}
//----------------------------------------------------------------------------------------------

//__Prototype Tracking Algorithm________________________________________________________________
int prototype_tracking(int argc,
                       char* argv[]) {
  const auto options = reader::parse_input(argc, argv);
  const auto detector_map = reader::import_detector_map(options.geometry_map_file);
  const auto time_resolution_map = reader::import_time_resolution_map(options.geometry_time_file);

  plot::init();
  geometry::open(options.geometry_file, options.default_time_error, time_resolution_map);

  std::cout << "Begin Tracking in " << options.data_directory << ":\n\n";
  plot::histogram histogram("chi_squared", "Chi-Squared Distribution", "chi^2/dof", "Track Count", 100, 0, 10);
  for (const auto& path : reader::root::search_directory(options.data_directory)) {

    std::cout << "FILE: " << path << "\n\n";
    for (const auto& event : reader::root::import_events(path, options, detector_map)) {
      plot::canvas canvas(path);
      std::cout << "EVENT: " << canvas.name() << "\n";

      const auto compressed_event = analysis::compress(event, options.compression_size);
      const auto density = analysis::geometric_event_density(compressed_event);
      std::cout << "Geometric Event Density: " << density << "\n";

      if (density > 2) continue;

      canvas.add_points(compressed_event, 0.8, plot::color::BLUE);
      draw_detector_centers(canvas);

      const auto tracks = find_tracks(compressed_event, options);
      for (const auto& track : tracks) {
        draw_track(canvas, track);
        histogram.insert(track.chi_squared_per_dof());
        std::cout << track << "\n";
      }
      canvas.draw();

      std::cout << "\n" << std::string(99, '=') << "\n\n";
    }
  }
  histogram.draw();
  histogram.save(".temp/save_hist.root");
  geometry::close();
  plot::end();
  return 0;
}
//----------------------------------------------------------------------------------------------

//__Quiet Prototype Tracking Algorithm__________________________________________________________
int quiet_prototype_tracking(int argc,
                             char* argv[]) {
  const auto options = reader::parse_input(argc, argv);
  const auto detector_map = reader::import_detector_map(options.geometry_map_file);
  const auto time_resolution_map = reader::import_time_resolution_map(options.geometry_time_file);

  plot::init(false);
  geometry::open(options.geometry_file, options.default_time_error, time_resolution_map);

  std::cout << "Begin Tracking in " << options.data_directory << ":\n\n";
  const auto statistics_path_prefix = options.data_directory + "/statistics";

  for (const auto& path : reader::root::search_directory(options.data_directory)) {
    std::cout << "Path: " << path << "\n";

    const auto imported_events = reader::root::import_events(path, options, detector_map);
    if (imported_events.size() == 0)
      continue;

    const auto statistics_path = statistics_path_prefix /*std::to_string(path_counter)*/ + ".root";
    plot::histogram chi_squared_histogram("chi_squared",
      "Chi-Squared Distribution", "chi^2/dof", "Track Count",
      200, 0, 10);
    plot::histogram beta_histogram("beta",
      "Beta Distribution", "beta", "Track Count",
      200, 0, 2);
    plot::histogram event_density_histogram("event_density",
      "Event Density Distribution", "Track Count", "Event Count",
      100, 0, 100);

    uint_fast64_t event_counter{};
    for (const auto& event : imported_events) {
      std::cout << "Event " << event_counter << " with " << event.size() << " hits.\n";

      if (event.empty()) {
        ++event_counter;
        continue;
      }

      const auto compressed_event = analysis::compress(event, options.compression_size);
      const auto compression_gain = event.size() / static_cast<type::real>(compressed_event.size());
      std::cout << "  Compression Gain: " << compression_gain << "\n";

      const auto density = modified_geometry_event_density(compressed_event);
      std::cout << "  Event Density: " << density * 100.0L << "%\n";

      if (density >= options.event_density_limit) {
        ++event_counter;
        continue;
      }

      const auto tracks = find_tracks(compressed_event, options);
      for (const auto& track : tracks) {
        chi_squared_histogram.insert(track.chi_squared_per_dof());
        beta_histogram.insert(track.beta());
      }

      if (tracks.size() > 1000) {
        const auto name = "event" + std::to_string(event_counter);
        plot::histogram individual_chi_square_histogram(name + "_chi_squared",
          name + " Chi-Squared Distribution", "chi^2/dof", "Track Count",
          200, 0, 10);
        plot::canvas individual_canvas(name);
        for (const auto& track : tracks) {
          individual_chi_square_histogram.insert(track.chi_squared_per_dof());
          draw_track(individual_canvas, track);
        }
        individual_chi_square_histogram.save(statistics_path);
        individual_canvas.draw();
        individual_canvas.save(statistics_path);
      }

      event_density_histogram.insert(tracks.size());

      std::cout << "  Track Count: "
                << tracks.size() << "\n";
      std::cout << "  Track Density: "
                << tracks.size() / static_cast<type::real>(event.size()) * 100.0L << "%\n";

      ++event_counter;
    }
    chi_squared_histogram.save(statistics_path);
    beta_histogram.save(statistics_path);
    event_density_histogram.save(statistics_path);
  }

  geometry::close();
  plot::end();
  return 0;
}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

//__Main Function: Prototype Tracker____________________________________________________________
int main(int argc,
         char* argv[]) {
  return MATHUSLA::prototype_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------

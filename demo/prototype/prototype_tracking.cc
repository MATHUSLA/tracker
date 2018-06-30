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

#include <tracker/analysis/analysis.hh>
#include <tracker/analysis/monte_carlo.hh>
#include <tracker/analysis/track.hh>
#include <tracker/analysis/vertex.hh>
#include <tracker/geometry.hh>
#include <tracker/plot.hh>
#include <tracker/reader.hh>

#include <tracker/util/io.hh>

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
  analysis::event combined_rpc_hits;
  analysis::full_event original_rpc_hits;
  const auto optimized_event = combine_rpc_hits(event, combined_rpc_hits, original_rpc_hits);
  const auto layers          = analysis::partition(optimized_event, options.layer_axis, options.layer_depth);
  const auto seeds           = analysis::seed(options.seed_size, layers, options.line_width);
  const auto tracking_vector = reset_seeds(analysis::join_all(seeds), combined_rpc_hits, original_rpc_hits);
  return analysis::overlap_fit_seeds(tracking_vector);
}
//----------------------------------------------------------------------------------------------

//__Print Bar___________________________________________________________________________________
void print_bar(const size_t count=99,
               const char b='=') {
  std::cout << "\n" << std::string(count, b) << "\n\n";
}
//----------------------------------------------------------------------------------------------

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
}
//----------------------------------------------------------------------------------------------

//__Draw MC Tracks to Canvas____________________________________________________________________
void draw_mc_tracks(plot::canvas& canvas,
                    const analysis::mc::track_vector& tracks) {
  for (const auto& track : tracks) {
    const auto hits = track.hits;
    const auto size = hits.size();
    for (std::size_t i = 0; i < size - 1; ++i) {
      canvas.add_point(type::reduce_to_r3(hits[i]), 0.8, plot::color::BLUE);
      canvas.add_line(type::reduce_to_r3(hits[i]),
                      type::reduce_to_r3(hits[i+1]),
                      1,
                      plot::color::BLUE);
    }
    canvas.add_point(type::reduce_to_r3(hits.back()), 0.8, plot::color::BLUE);
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
  const type::r3_point guess_center{guess.x.value, guess.y.value, guess.z.value};
  canvas.add_point(guess_center, 1.5, plot::color::RED);

  for (const auto& track : vertex.tracks())
    canvas.add_line(track(point.z), track.back(), 1.5, plot::color::GREEN);
}
//----------------------------------------------------------------------------------------------

//__Show and Add Tracks to Statistics___________________________________________________________
void save_tracks(const analysis::track_vector& tracks,
                 plot::canvas& canvas,
                 plot::histogram_collection& histograms,
                 bool verbose) {
  for (const auto& track : tracks) {
    histograms["track_chi_squared"].insert(track.chi_squared_per_dof());
    histograms["beta"].insert(track.beta());
    histograms["beta_error"].insert(track.beta_error());
    histograms["track_size"].insert(track.size());
    draw_track(canvas, track);
    //if (verbose)
      //std::cout << track << "\n";
  }
}
//----------------------------------------------------------------------------------------------

//__Prototype Tracking Algorithm________________________________________________________________
int prototype_tracking(int argc,
                       char* argv[]) {
  const auto options = reader::parse_input(argc, argv);
  const auto detector_map = reader::import_detector_map(options.geometry_map_file);
  const auto time_resolution_map = reader::import_time_resolution_map(options.geometry_time_file);

  plot::init(options.verbose_output);
  geometry::open(options.geometry_file, options.default_time_error, time_resolution_map);

  std::cout << "Begin Tracking in " << options.data_directory << ":\n\n";
  const auto statistics_path_prefix = options.statistics_directory + "/" + options.statistics_file_prefix;

  uint_fast64_t path_counter{};
  for (const auto& path : reader::root::search_directory(options.data_directory, options.data_file_extension)) {
    const auto path_counter_string = std::to_string(path_counter++);

    print_bar();
    std::cout << "Read Path: " << path << "\n";

    const auto event_bundle = reader::root::import_event_mc_bundle(path, options, detector_map);
    const auto imported_events = event_bundle.events;
    const auto import_size = imported_events.size();
    if (import_size == 0)
      continue;

    const auto mc_imported_events = event_bundle.true_events;

    const auto statistics_path = statistics_path_prefix + path_counter_string + "." + options.statistics_file_extension;
    plot::histogram_collection histograms({
      {"track_chi_squared",   "Track Chi-Squared Distribution",    "chi^2/dof",    "Track Count",  200, 0, 10},
      {"vertex_chi_squared",  "Track Chi-Squared Distribution",    "chi^2/dof",    "Vertex Count", 200, 0, 10},
      {"beta",                "Track Beta Distribution",           "beta",         "Track Count",  200, 0,  2},
      {"beta_error",          "Track Beta Error Distribution",     "beta error",   "Track Count",  200, 0,  2},
      {"beta_with_cut",       "Track Beta Distribution With Cut",  "beta",         "Track Count",  200, 0,  2},
      {"track_count",         "Track Count Distribution",          "Track Count",  "Event Count",   51, 0, 50},
      {"vertex_count",        "Vertex Count Distribution",         "Vertex Count", "Event Count",  100, 0, 50},
      {"track_size",          "Track Size Distribution",           "Hit Count",    "Track Count",   11, 0, 10},
      {"non_track_hit_count", "Non-Track Hit Count Distribution",  "Hit Count",    "Event Count",  100, 0, 50},
      {"vertex_t_error",      "Vertex T Error Distribution",       "t error (" + units::time_string   + ")", "Vertex Count", 100, 0, 20},
      {"vertex_x_error",      "Vertex X Error Distribution",       "x error (" + units::length_string + ")", "Vertex Count", 100, 0, 100},
      {"vertex_y_error",      "Vertex Y Error Distribution",       "y error (" + units::length_string + ")", "Vertex Count", 100, 0, 100},
      {"vertex_z_error",      "Vertex Z Error Distribution",       "z error (" + units::length_string + ")", "Vertex Count", 100, 0, 100},
    });

    for (uint_fast64_t event_counter{}; event_counter < import_size; ++event_counter) {
      const auto& event = imported_events[event_counter];
      const auto event_size = event.size();
      const auto event_counter_string = std::to_string(event_counter);

      const auto compressed_event = analysis::compress(event, options.compression_size);
      const auto compression_gain = event_size / static_cast<type::real>(compressed_event.size());

      if (event.empty() || compression_gain == event_size)
        continue;

      std::cout << "Event " << event_counter << " with " << event_size << " hits.\n";

      const auto event_density = modified_geometry_event_density(compressed_event);
      if (options.verbose_output)
        std::cout << "  Compression Gain: " << compression_gain << "\n"
                  << "  Event Density: "    << event_density * 100.0L << " %\n";

      if (event_density >= options.event_density_limit)
        continue;

      plot::canvas canvas(path + event_counter_string);
      draw_detector_centers(canvas);
      draw_mc_tracks(canvas, analysis::mc::convert(mc_imported_events[event_counter]));

      const auto tracks = find_tracks(compressed_event, options);
      const auto tracks_size = tracks.size();

      histograms["track_count"].insert(tracks_size);

      if (tracks_size > 100) {
        /*
        const auto name = "event" + event_counter_string;
        plot::canvas individual_canvas(name);
        plot::histogram_collection individual_histograms(name + "_", {
          {"track_chi_squared", name + " Chi-Squared Distribution", "chi^2/dof", "Track Count", 200, 0, 10},
          {"beta", name + "Beta Distribution", "beta", "Track Count", 200, 0, 2}
        });
        save_tracks(tracks, individual_canvas, individual_histograms, options.verbose_output);
        // FIXME: should not need a draw before a save
        individual_canvas.draw();
        plot::save_all(statistics_path, individual_canvas, individual_histograms);
        */
      } else {
        canvas.add_points(compressed_event, 0.8, plot::color::BLACK);
        save_tracks(tracks, canvas, histograms, options.verbose_output);
        // TODO: to remove
        for (const auto& track : tracks) {
          util::io::print_range(track.event()) << "\n";
        }
      }

      histograms["track_count"].insert(tracks_size);

      if (options.verbose_output)
        std::cout << "  Track Count: "   << tracks.size() << "\n"
                  << "  Track Density: " << tracks.size() / static_cast<type::real>(event.size()) * 100.0L << " %\n";

      //const analysis::vertex vertex(tracks);
      //if (options.verbose_output) {
        //draw_vertex_and_guess(canvas, vertex);
        // std::cout << vertex << "\n";
      //}

      canvas.draw();

      std::clog << "event: " << event_counter << "\n";
      for (const auto& track : tracks) {
        util::io::print_range(track.event(), " ", "", std::clog) << "\n";
      }
    }
    histograms.draw_all();
    histograms.save_all(statistics_path);
  }

  print_bar();
  geometry::close();
  plot::end();
  return 0;
}
//----------------------------------------------------------------------------------------------

//__Silent Prototype Tracking Algorithm_________________________________________________________
int silent_prototype_tracking(int argc,
                              char* argv[]) {
  util::io::remove_buffer(std::cout, std::cerr /*, std::clog*/);
  return prototype_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

//__Main Function: Prototype Tracker____________________________________________________________
int main(int argc,
         char* argv[]) {
  return MATHUSLA::silent_prototype_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------

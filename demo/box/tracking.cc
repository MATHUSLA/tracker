/*
 * demo/box/tracking.cc
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

#include <tracker/util/io.hh>

#include "geometry.hh"
#include "io.hh"

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace mc       = analysis::mc;
namespace geometry = MATHUSLA::TRACKER::geometry;
namespace plot     = MATHUSLA::TRACKER::plot;
namespace reader   = MATHUSLA::TRACKER::reader;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

//__Find Primary Tracks for Box_________________________________________________________________
const analysis::track_vector find_primary_tracks(const analysis::full_event& event,
                                                 const reader::tracking_options& options,
                                                 plot::canvas& canvas,
                                                 analysis::full_event& non_track_points) {
  const auto layers = analysis::partition(event, options.layer_axis, options.layer_depth);
  const auto seeds = analysis::seed(options.seed_size, layers, analysis::seed_heuristic::double_cone{options.line_width});

  /*
  for (const auto& seed : seeds) {
    for (std::size_t i{}; i < seed.size() - 1; ++i) {
      canvas.add_line(type::reduce_to_r4(seed[i]), type::reduce_to_r4(seed[i+1]), 1, plot::color::BLACK);
    }
  }
  */

  auto first_tracks = analysis::independent_fit_seeds(analysis::join_all(seeds), options.layer_axis);

  for (auto& track : first_tracks)
    track.prune_on_chi_squared(20.0L);

  const auto out = analysis::overlap_fit_tracks(first_tracks, 1UL);
  non_track_points = analysis::non_tracked_points(event, out, true);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Find Tracks for Box_________________________________________________________________________
const analysis::track_vector find_tracks(const analysis::full_event& event,
                                         const reader::tracking_options& options,
                                         const type::real limit_chi_squared,
                                         const std::size_t overlap,
                                         plot::canvas& canvas,
                                         analysis::full_event& non_track_points) {
  const auto layers = analysis::partition(event, options.layer_axis, options.layer_depth);
  const auto seeds = analysis::seed(options.seed_size, layers, analysis::seed_heuristic::double_cone{options.line_width});

  /*
  for (const auto& seed : seeds) {
    for (std::size_t i{}; i < seed.size() - 1; ++i) {
      canvas.add_line(type::reduce_to_r4(seed[i]), type::reduce_to_r4(seed[i+1]), 1, plot::color::BLACK);
    }
  }
  */

  auto first_tracks = analysis::independent_fit_seeds(analysis::join_all(seeds), options.layer_axis);

  for (auto& track : first_tracks)
    track.prune_on_chi_squared(limit_chi_squared);

  const auto out = analysis::overlap_fit_tracks(first_tracks, overlap);
  non_track_points = analysis::non_tracked_points(event, out, true);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Box Tracking Algorithm______________________________________________________________________
int box_tracking(int argc,
                 char* argv[]) {
  box::extension_parser extension{};
  const auto options = reader::parse_input(argc, argv, extension);

  plot::init(options.draw_events);
  geometry::open(options.geometry_file, options.default_time_error);
  box::geometry::import(extension);

  std::cout << "Begin Tracking in " << options.data_directory << ":\n\n";
  const auto statistics_path_prefix = options.statistics_directory + "/" + options.statistics_file_prefix;
  const plot::value_tag filetype_tag("FILETYPE", "MATHUSLA TRACKING STATFILE");
  const plot::value_tag project_tag("PROJECT", "Box");

  std::size_t path_counter{};
  for (const auto& path : reader::root::search_directory(options.data_directory, options.data_file_extension)) {
    const auto path_counter_string = std::to_string(path_counter++);
    const auto statistics_save_path = statistics_path_prefix
                                    + path_counter_string
                                    + "." + options.statistics_file_extension;

    box::print_bar();
    std::cout << "Read Path: " << path << "\n";

    const auto event_bundle = reader::root::import_event_mc_bundle(path,
      options.data_track_id_key, options.data_t_key, options.data_x_key, options.data_y_key, options.data_z_key);

    const auto imported_events = event_bundle.events;
    const auto import_size = imported_events.size();
    if (import_size == 0UL)
      continue;

    const auto mc_imported_events = event_bundle.true_events;

    analysis::track::tree track_tree{"track_tree", "MATHUSLA Track Tree"};
    analysis::vertex::tree vertex_tree{"vertex_tree", "MATHUSLA Vertex Tree"};

    for (std::size_t event_counter{}; event_counter < import_size; ++event_counter) {
      const auto event = analysis::add_width<box::geometry>(imported_events[event_counter]);
      const auto event_size = event.size();
      const auto event_counter_string = std::to_string(event_counter);

      const auto compressed_event = options.time_smearing ? mc::time_smear<box::geometry>(mc::compress<box::geometry>(event))
                                                          : mc::compress<box::geometry>(event);
      const auto compression_size = event_size / static_cast<type::real>(compressed_event.size());

      if (event_size == 0UL || compression_size == event_size)
        continue;

      const auto altered_event =
        box::geometry::restrict_layer_count(
          mc::add_noise<box::geometry>(
            mc::use_efficiency(compressed_event, options.simulated_efficiency),
            options.simulated_noise_rate,
            options.event_time_window.begin,
            options.event_time_window.end),
          extension.layer_count);

      const auto event_density = box::geometry::event_density(altered_event);
      box::print_event_summary(event_counter, altered_event.size(), compression_size, event_density);

      if (event_density >= options.event_density_limit)
        continue;

      plot::canvas canvas("event" + event_counter_string, path + event_counter_string);
      if (options.draw_events) {
        box::draw_detector(canvas, extension.layer_count);
        box::draw_mc_tracks(canvas, mc::convert_events(mc_imported_events[event_counter]));
        for (const auto hit : altered_event)
          canvas.add_point(type::reduce_to_r3(hit), 0.8, plot::color::BLACK);
      }

      analysis::full_event non_primary_track_points, non_secondary_track_points;
      auto tracks = find_primary_tracks(altered_event, options, canvas, non_primary_track_points);

      auto secondary_tracks = find_tracks(non_primary_track_points,
                                          options,
                                          100.0L,
                                          0UL,
                                          canvas,
                                          non_secondary_track_points);
      tracks.reserve(tracks.size() + secondary_tracks.size());
      tracks.insert(tracks.cend(),
                    std::make_move_iterator(secondary_tracks.cbegin()),
                    std::make_move_iterator(secondary_tracks.cend()));

      box::save_tracks(tracks, canvas, track_tree, options);
      box::print_tracking_summary(event, tracks);

      analysis::track_vector converged_tracks;
      converged_tracks.reserve(tracks.size());
      std::copy_if(std::make_move_iterator(tracks.begin()),
                   std::make_move_iterator(tracks.end()),
                   std::back_inserter(converged_tracks),
                   [](const auto& track) { return track.fit_converged(); });

      box::save_vertices(analysis::pairwise_fit_tracks(converged_tracks), canvas, vertex_tree, options);

      canvas.draw();
    }

    plot::save_all(statistics_save_path,
      filetype_tag,
      project_tag,
      plot::value_tag{"DATAPATH", path},
      plot::value_tag{"EVENTS", std::to_string(import_size)},
      plot::value_tag{"EFFICIENCY", std::to_string(options.simulated_efficiency)},
      plot::value_tag{"NOISE", std::to_string(options.simulated_noise_rate * units::time) + " / " + units::time_string},
      box::geometry::value_tags());
    track_tree.save(statistics_save_path);
    vertex_tree.save(statistics_save_path);
  }

  box::print_bar();
  geometry::close();
  plot::end();
  return 0;
}
//----------------------------------------------------------------------------------------------

//__Silent Box Tracking Algorithm_______________________________________________________________
int silent_box_tracking(int argc,
                        char* argv[]) {
  util::io::remove_buffer(std::cout, std::cerr, std::clog);
  return box_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

//__Main Function: Box Tracker__________________________________________________________________
int main(int argc,
         char* argv[]) {
  return MATHUSLA::box_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------

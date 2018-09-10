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
#include <tracker/script.hh>

#include "geometry.hh"
#include "io.hh"

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace mc       = analysis::mc;
namespace geometry = MATHUSLA::TRACKER::geometry;
namespace plot     = MATHUSLA::TRACKER::plot;
namespace reader   = MATHUSLA::TRACKER::reader;
namespace script   = MATHUSLA::TRACKER::script;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

//__Alter Event_________________________________________________________________________________
const analysis::full_event alter_event(const analysis::full_event& event,
                                       const script::tracking_options& options,
                                       const box::io::extension_parser& extension) {
  return box::geometry::restrict_layer_count(
           mc::use_efficiency(
             mc::add_noise<box::geometry>(
               event,
               options.simulated_noise_rate,
               options.event_time_window.begin,
               options.event_time_window.end,
               box::geometry::full(extension.layer_count)),
             options.simulated_efficiency),
           extension.layer_count);
}
//----------------------------------------------------------------------------------------------

//__Find Tracks for Box_________________________________________________________________________
const analysis::track_vector find_tracks(const analysis::full_event& event,
                                         const script::tracking_options& options,
                                         const type::real limit_chi_squared,
                                         const std::size_t overlap,
                                         plot::canvas& canvas,
                                         analysis::full_event& non_track_points) {
  namespace ash = analysis::seed_heuristic;
  const auto layers = analysis::partition(event, options.layer_axis, options.layer_depth);
  const auto seeds = analysis::seed(
    options.seed_size,
    layers,
    ash::all<ash::double_cone, ash::piecewise_speed>{
      {options.line_width},
      {0.9L * units::speed_of_light}});

  /*
  for (const auto& seed : seeds) {
    for (std::size_t i{}; i < seed.size() - 1UL; ++i)
      canvas.add_line(type::reduce_to_r4(seed[i]), type::reduce_to_r4(seed[i + 1UL]), 1, plot::color::BLACK);
    for (const auto& point : seed)
      std::cout << type::reduce_to_r4(point) << " ";
    std::cout << "\n";
  }
  */

  const auto joined = analysis::join_all(seeds);
  auto first_tracks = analysis::independent_fit_seeds(joined, options.layer_axis);

  for (auto& track : first_tracks)
    track.prune_on_chi_squared(limit_chi_squared);

  const auto out = analysis::overlap_fit_tracks(first_tracks, overlap);
  non_track_points = analysis::non_tracked_points(event, out, true);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Track an Event Bundle_______________________________________________________________________
void track_event_bundle(const script::path_vector& paths,
                        const mc::event_vector_bundle& bundle,
                        const script::tracking_options& options,
                        const box::io::extension_parser& extension,
                        const script::path_type& save_path) {
  static const plot::value_tag filetype_tag("FILETYPE", "MATHUSLA TRACKING STATFILE");
  static const plot::value_tag project_tag("PROJECT", "Box");

  const auto imported_events = bundle.events;
  const auto import_size = imported_events.size();
  if (import_size == 0UL)
    return;

  const auto mc_imported_events = bundle.true_events;

  analysis::track::tree track_tree{"track_tree", "MATHUSLA Track Tree"};
  analysis::vertex::tree vertex_tree{"vertex_tree", "MATHUSLA Vertex Tree"};
  track_tree.add_friend(vertex_tree, "vertex");
  vertex_tree.add_friend(track_tree, "track");

  std::cout << "Event Count: " << import_size << "\n";

  for (std::size_t event_counter{}; event_counter < import_size; ++event_counter) {
    const auto event = analysis::add_width<box::geometry>(imported_events[event_counter]);
    const auto event_size = event.size();
    const auto event_counter_string = std::to_string(event_counter);

    const auto compressed_event = options.time_smearing ? mc::time_smear<box::geometry>(mc::compress<box::geometry>(event))
                                                        : mc::compress<box::geometry>(event);
    const auto compression_size = event_size / static_cast<type::real>(compressed_event.size());

    if (event_size == 0UL || compression_size == event_size) {
      track_tree.fill();
      vertex_tree.fill();
      continue;
    }

    const auto altered_event = alter_event(compressed_event, options, extension);
    const auto event_density = box::geometry::event_density(altered_event);
    box::io::print_event_summary(event_counter, altered_event.size(), compression_size, event_density);

    if (event_density >= options.event_density_limit) {
      track_tree.fill();
      vertex_tree.fill();
      continue;
    }

    // FIXME: better title for canvas
    std::string canvas_title{"event"};
    for (const auto& path : paths) {
      canvas_title.push_back(' ');
      canvas_title.insert(canvas_title.cend(), path.begin(), path.end());
    }

    plot::canvas canvas("event" + event_counter_string, canvas_title + " | " + event_counter_string);
    if (options.draw_events) {
      box::io::draw_detector(canvas, extension.layer_count);
      box::io::draw_mc_tracks(canvas, mc::convert_events(mc_imported_events[event_counter]));
      for (const auto& hit : altered_event)
        canvas.add_point(type::reduce_to_r3(hit), 0.8, plot::color::BLACK);
    }

    analysis::full_event non_primary_track_points, non_secondary_track_points;
    auto tracks = find_tracks(altered_event,
                              options,
                              20.0L,
                              1UL,
                              canvas,
                              non_primary_track_points);
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

    box::io::save_tracks(tracks, canvas, track_tree, options);
    box::io::print_tracking_summary(event, tracks);

    analysis::track_vector converged_tracks;
    converged_tracks.reserve(tracks.size());
    std::copy_if(std::make_move_iterator(tracks.begin()),
                 std::make_move_iterator(tracks.end()),
                 std::back_inserter(converged_tracks),
                 [](const auto& track) { return track.fit_converged(); });

    box::io::save_vertices(analysis::pairwise_fit_tracks(converged_tracks), canvas, vertex_tree, options);

    canvas.draw();
  }

  plot::save_all(save_path,
    filetype_tag,
    plot::value_tag{"TIMESTAMP", util::time::GetString("%c %Z")},
    project_tag,
    plot::value_tag{"EVENTS", std::to_string(import_size)},
    plot::value_tag{"EFFICIENCY", std::to_string(options.simulated_efficiency)},
    plot::value_tag{"NOISE", std::to_string(options.simulated_noise_rate * units::time) + " / " + units::time_string},
    box::geometry::value_tags(),
    box::io::data_paths_value_tags(paths, options.data_timing_offsets));

  box::io::save_files(save_path, track_tree, vertex_tree, paths, options.merge_input);
}
//----------------------------------------------------------------------------------------------

//__Box Tracking Algorithm______________________________________________________________________
int box_tracking(int argc,
                 char* argv[]) {
  box::io::extension_parser extension{};
  const auto options = script::parse_command_line(argc, argv, extension);

  plot::init(options.draw_events);
  geometry::open(options.geometry_file, options.default_time_error);
  box::geometry::import(extension);

  box::io::print_tracking_directories(options.data_directories);

  std::size_t path_counter{};
  const auto statistics_path_prefix = box::io::add_statistics_path(options);
  for (const auto& paths : reader::root::transpose_search_directories(options.data_directories, options.data_file_extension)) {
    const auto path_counter_string = std::to_string(path_counter++);
    const auto statistics_save_path = statistics_path_prefix
                                    + path_counter_string
                                    + "."
                                    + options.statistics_file_extension;

    box::io::print_bar();
    box::io::print_tracking_paths(paths);

    track_event_bundle(paths,
      reader::root::import_event_mc_bundle(paths,
        options.data_timing_offsets,
        options.data_track_id_key,
        options.data_t_key,
        options.data_x_key,
        options.data_y_key,
        options.data_z_key),
      options,
      extension,
      statistics_save_path);
  }

  box::io::print_bar();
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

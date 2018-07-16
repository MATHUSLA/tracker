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

#include <tracker/util/index_vector.hh>
#include <tracker/util/bit_vector.hh>
#include <tracker/util/io.hh>

#include "geometry.hh"
// #include "logging.hh"

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace geometry = MATHUSLA::TRACKER::geometry;
namespace plot     = MATHUSLA::TRACKER::plot;
namespace reader   = MATHUSLA::TRACKER::reader;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

//__Find Primary Tracks for Box_________________________________________________________________
const analysis::track_vector find_primary_tracks(const analysis::event& event,
                                                 const reader::tracking_options& options,
                                                 plot::canvas& canvas,
                                                 analysis::event& non_track_points) {
  const auto layers = analysis::partition(event, options.layer_axis, options.layer_depth);
  const auto seeds = analysis::seed(options.seed_size, layers, analysis::topology::double_cone{options.line_width});

  /*
  for (const auto seed : seeds) {
    for (std::size_t i{}; i < seed.size() - 1; ++i) {
      canvas.add_line(type::reduce_to_r4(seed[i]), type::reduce_to_r4(seed[i+1]), 1, plot::color::BLACK);
    }
  }
  */

  auto out = analysis::overlap_fit_seeds(analysis::join_all(seeds), options.layer_axis, 1UL);

  // TODO: improve efficiency
  const auto size = event.size();
  util::bit_vector save_list(size);
  for (auto& track : out) {
    for (const auto& point : track.event()) {
      const auto search = util::algorithm::range_binary_find_first(event,
                                                                   point,
                                                                   type::t_ordered<analysis::hit>{});
      if (search != event.cend())
        save_list.set(static_cast<std::size_t>(search - event.cbegin()));
    }
  }

  non_track_points.clear();
  non_track_points.reserve(size - save_list.count());
  save_list.unset_conditional_push_back(event, non_track_points);

  return out;
}
//----------------------------------------------------------------------------------------------

//__Box Tracking Algorithm______________________________________________________________________
int box_tracking(int argc,
                 char* argv[]) {
  const auto options = reader::parse_input(argc, argv);
  const auto detector_map = reader::import_detector_map(options.geometry_map_file);

  plot::init(options.draw_events);
  geometry::open(options.geometry_file,
                 options.default_time_error,
                 reader::import_time_resolution_map(options.geometry_time_file));

  std::size_t path_counter{};
  for (const auto& path : reader::root::search_directory(options.data_directory, options.data_file_extension)) {

    for (;;) {

    }

  }

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

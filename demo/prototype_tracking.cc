/*
 * demo/prototype_tracking.cc
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
#include <tracker/reader.hh>
#include <tracker/plot.hh>
#include <tracker/units.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/bit_vector.hh>
#include <tracker/util/io.hh>

namespace MATHUSLA {

//__Combine Pair of Hits if they Occur in Overlapping RPCs______________________________________
const TRACKER::geometry::box_volume combine_rpc_volume_pair(const TRACKER::geometry::box_volume& first,
                                                            const TRACKER::geometry::box_volume& second) {
  using namespace TRACKER;
  auto out = geometry::coordinatewise_intersection(first, second);
  const auto union_volume = geometry::coordinatewise_union(first, second);
  out.min.z = union_volume.min.z;
  out.center.z = union_volume.center.z;
  out.max.z = union_volume.max.z;
  return out;
}
//----------------------------------------------------------------------------------------------

//__Check If RPC Combine Created a Valid Strip Overlap__________________________________________
constexpr bool was_combine_successful(const TRACKER::geometry::box_volume& combined) {
  return (combined.min.x || combined.max.x) && (combined.min.y || combined.max.y);
}
//----------------------------------------------------------------------------------------------

//__Check If RPC Combine Created a Valid Strip Overlap__________________________________________
const TRACKER::analysis::full_hit construct_hit(const type::real top_time,
                                                const type::real bottom_time,
                                                const std::string& top_volume,
                                                const std::string& bottom_volume,
                                                const TRACKER::geometry::box_volume& combined) {
  using namespace TRACKER;
  const type::r4_point errors{
    std::hypot(geometry::time_resolution_of(top_volume),
               geometry::time_resolution_of(bottom_volume)) / 2.0L,
    combined.max.x - combined.min.x,
    combined.max.y - combined.min.y,
    combined.max.z - combined.min.z};
  return {(top_time + bottom_time) / 2.0L, combined.center.x, combined.center.y, combined.center.z,
          errors};
}
//----------------------------------------------------------------------------------------------

//__Combine Hits if they Occur in Overlapping RPCs______________________________________________
const TRACKER::analysis::full_event combine_rpc_hits(const TRACKER::analysis::event& points,
                                                     const type::real time_threshold) {
  using namespace TRACKER;
  using namespace util::math;

  static const analysis::real z_lower = 24.0L * units::length;
  static const analysis::real z_upper = 45.0L * units::length;

  const auto size = points.size();
  analysis::full_event event;
  event.reserve(size);

  const auto parts = analysis::partition(points, type::Coordinate::Z, z_lower).parts;
  const auto partition_size = parts.size();

  size_t layer_index = 0;
  for (; layer_index < partition_size - 1; ++layer_index) {
    const auto top = parts[layer_index];
    const auto bottom = parts[layer_index + 1];

    if (within(top.front().z, bottom.back().z, z_lower, z_upper)) {
      const auto top_size = top.size();
      const auto bottom_size = bottom.size();
      util::bit_vector discard_list(bottom_size);

      size_t top_index = 0, bottom_index = 0;
      while (top_index < top_size) {
        bottom_index = discard_list.first_unset(bottom_index);
        const auto top_point = top[top_index];
        const auto bottom_point = bottom[bottom_index];

        if (bottom_index == bottom_size) {
          event.push_back(analysis::find_errors(top_point));
          ++top_index;
          bottom_index = 0;
          continue;
        }

        const auto top_volume = geometry::volume(top_point);
        const auto bottom_volume = geometry::volume(bottom_point);
        const auto combined = combine_rpc_volume_pair(
          geometry::limits_of(top_volume),
          geometry::limits_of(bottom_volume));

        if (was_combine_successful(combined) && within(top_point.t, bottom_point.t, time_threshold)) {
          event.push_back(construct_hit(top_point.t, bottom_point.t, top_volume, bottom_volume, combined));
          discard_list.set(bottom_index);
          ++top_index;
          bottom_index = 0;
        } else {
          ++bottom_index;
        }
      }

      for (; top_index < top_size; ++top_index) {
        event.push_back(analysis::find_errors(top[top_index]));
      }
      for (bottom_index = 0; bottom_index < bottom_size; ++bottom_index) {
        if (!discard_list[bottom_index])
          event.push_back(analysis::find_errors(bottom[bottom_index]));
      }
      ++layer_index;
    } else {
      util::algorithm::back_insert_transform(top, event,
        [](const auto& part){ return analysis::find_errors(part); });
    }

  }

  if (layer_index == partition_size - 1) {
    util::algorithm::back_insert_transform(parts.back(), event,
      [](const auto& part){ return analysis::find_errors(part); });
  }

  util::algorithm::reverse(event);
  return event;
}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

//__Main Function: Prototype Tracker____________________________________________________________
int main(int argc, char* argv[]) {
  using namespace MATHUSLA;
  using namespace MATHUSLA::TRACKER;

  const auto options = reader::parse_input(argc, argv);
  const auto detector_map = reader::import_detector_map(options.geometry_map_file);
  const auto time_resolution_map = reader::import_time_resolution_map(options.geometry_time_file);

  plot::init();
  geometry::open(options.geometry_file, options.default_time_error, time_resolution_map);
  for (const auto& path : reader::root::search_directory(options.root_directory)) {
    plot::canvas canvas(path);
    std::cout << path << "\n";
    for (const auto& event : reader::root::import_events(path, options, detector_map)) {

      const auto collapsed_event = analysis::collapse(event, options.collapse_size);

      canvas.add_points(collapsed_event, 1.5, {90, 90, 90});
      for (const auto& name : geometry::full_structure_except({"world", "Sandstone", "Marl", "Mix", "Earth"})) {
        const auto limits = geometry::limits_of(name);
        canvas.add_point(limits.center, 0.25, plot::color::MAGENTA);
      }

      util::io::print_range(collapsed_event, "\n", "NEW ") << "\n\n";

      const auto full_event_points = combine_rpc_hits(collapsed_event, 2 * units::time);

      uint_fast8_t index = 0, step = 230 / full_event_points.size();
      for (const auto& point : full_event_points) {
        canvas.add_box(type::reduce_to_r3(point),
                       point.error.x, point.error.y, point.error.z,
                       3, {index, index, index});
        index += step;
      }

      /*
      const auto layers = analysis::partition(combine_rpc_hits(collapsed_event, 2*units::time),
                                              options.layer_axis,
                                              options.layer_depth);

      const auto seeds = analysis::seed(options.seed_size,
                                        layers,
                                        options.line_width);

      const auto tracks = analysis::fit_seeds(analysis::join_all(seeds));
      for (const auto& track : tracks) {
        const auto& event = track.event();
        canvas.add_points(event, 0.5);
        for (const auto& point : event) {
          const auto limits = geometry::limits_of_volume(point);
          canvas.add_point(limits.center, 0.25, plot::color::BLUE);
          canvas.add_box(limits.min, limits.max, 2, plot::color::BLUE);
        }
        canvas.add_line(event.front(), event.back());
        std::cout << track << "\n";
        canvas.add_line(track.front(), track.back(), 1, plot::color::RED);
      }

      OLD CODE: vvvv
      ////////////////////////////////////////////////////////////////////////////////

      const auto layered_event = analysis::partition(collapsed_event, options.layer_axis, options.layer_depth);

      canvas.add_points(collapsed_event, 1.5, {90, 90, 90});

      util::io::print_range(event, "\n", "OLD ") << "\n\n";
      util::io::print_range(collapsed_event, "\n", "NEW ") << "\n\n";
      for (const auto& layer : layered_event.parts)
        util::io::print_range(layer, "\n", "LAYER ") << "\n\n\n";

      const auto seeds = analysis::seed(options.seed_size,
                                        layered_event,
                                        options.line_width);

      std::cout << "seeds (" << seeds.size() << "):\n\n";
      for (const auto& seed : seeds)
        util::io::print_range(seed) << "\n";

      util::io::newline();

      const auto joined = analysis::join_all(seeds);

      std::cout << "joined (" << joined.size() << "):\n\n";
      for (const auto& seed : joined)
        util::io::print_range(seed) << "\n";

      const auto tracks = analysis::fit_seeds(joined, {"MIGRAD", {}, false, -1, 0.5, 150});
      std::cout << "\nTRACK FITTING (" << tracks.size() << "):\n\n";
      for (const auto& track : tracks) {
        const auto& event = track.event();
        canvas.add_points(event, 0.5);
        for (const auto& point : event) {
          const auto limits = geometry::limits_of_volume(point);
          canvas.add_point(limits.center, 0.25, plot::color::BLUE);
          canvas.add_box(limits.min, limits.max, 2, plot::color::BLUE);
        }
        canvas.add_line(event.front(), event.back());
        std::cout << track << "\n";
        canvas.add_line(track.front(), track.back(), 1, plot::color::RED);
      }
      */

      canvas.draw();
    }
  }

  geometry::close();
  plot::end();

  return 0;
}
//----------------------------------------------------------------------------------------------

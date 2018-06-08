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
#include <tracker/plot.hh>
#include <tracker/reader.hh>
#include <tracker/units.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/bit_vector.hh>
#include <tracker/util/io.hh>

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace geometry = MATHUSLA::TRACKER::geometry;
namespace plot     = MATHUSLA::TRACKER::plot;
namespace reader   = MATHUSLA::TRACKER::reader;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

//__Combine Pair of Hits if they Occur in Overlapping RPCs______________________________________
const geometry::box_volume combine_rpc_volume_pair(const geometry::box_volume& first,
                                                   const geometry::box_volume& second) {
  auto out = geometry::coordinatewise_intersection(first, second);
  const auto union_volume = geometry::coordinatewise_union(first, second);
  out.min.z = union_volume.min.z;
  out.center.z = union_volume.center.z;
  out.max.z = union_volume.max.z;
  return out;
}
//----------------------------------------------------------------------------------------------

//__Check If RPC Combine Created a Valid Strip Overlap__________________________________________
constexpr bool was_combine_successful(const geometry::box_volume& combined) {
  return (combined.min.x || combined.max.x) && (combined.min.y || combined.max.y);
}
//----------------------------------------------------------------------------------------------

//__Check If RPC Combine Created a Valid Strip Overlap__________________________________________
const analysis::full_hit construct_hit(const type::real top_time,
                                       const type::real bottom_time,
                                       const std::string& top_volume,
                                       const std::string& bottom_volume,
                                       const geometry::box_volume& combined) {
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
const analysis::full_event combine_rpc_hits(const analysis::event& points,
                                            const type::real time_threshold,
                                            analysis::full_event& seeding_rpc_hits,
                                            analysis::full_event& tracking_rpc_hits) {
  static const analysis::real z_lower = 24.0L * units::length;
  static const analysis::real z_upper = 45.0L * units::length;
  using namespace util::math;

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
          event.push_back(analysis::add_errors(top_point));
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
          const auto new_hit = construct_hit(top_point.t, bottom_point.t, top_volume, bottom_volume, combined);
          event.push_back(new_hit);
          seeding_rpc_hits.push_back(new_hit);
          tracking_rpc_hits.push_back(analysis::add_errors(top_point));
          tracking_rpc_hits.push_back(analysis::add_errors(bottom_point));
          discard_list.set(bottom_index);
          ++top_index;
          bottom_index = 0;
        } else {
          ++bottom_index;
        }
      }

      for (; top_index < top_size; ++top_index) {
        event.push_back(analysis::add_errors(top[top_index]));
      }
      for (bottom_index = 0; bottom_index < bottom_size; ++bottom_index) {
        if (!discard_list[bottom_index])
          event.push_back(analysis::add_errors(bottom[bottom_index]));
      }
      ++layer_index;
    } else {
      util::algorithm::back_insert_transform(top, event,
        [](const auto& part){ return analysis::add_errors(part); });
    }

  }

  if (layer_index == partition_size - 1) {
    util::algorithm::back_insert_transform(parts.back(), event,
      [](const auto& part){ return analysis::add_errors(part); });
  }

  util::algorithm::reverse(event);
  return event;
}
//----------------------------------------------------------------------------------------------

//__Reset Seed Vector Using RPC Combination Hits________________________________________________
const analysis::full_event_vector reset_seeds(const analysis::full_event_vector& joined_seeds,
                                              const analysis::full_event& seeding_rpc_hits,
                                              const analysis::full_event& tracking_rpc_hits) {
  const auto size = seeding_rpc_hits.size();
  if (size == 0) return {};

  analysis::full_event_vector out;
  for (const auto& seed : joined_seeds) {
    out.push_back({});
    auto& back = out.back();
    for (const auto& hit : seed) {
      size_t rpc_index = 0;
      for (; rpc_index < size; ++rpc_index) {
        if (hit == seeding_rpc_hits[rpc_index]) {
          back.push_back(tracking_rpc_hits[2 * rpc_index]);
          back.push_back(tracking_rpc_hits[2 * rpc_index + 1]);
        }
      }
      if (rpc_index == size)
        back.push_back(hit);
    }
    type::t_sort(back);
  }

  return out;
}
//----------------------------------------------------------------------------------------------

//__Add Track to Canvas_________________________________________________________________________
void draw_track(plot::canvas& canvas,
                const analysis::track& track) {
  const auto& full_event = track.full_event();
  uint_fast8_t brightness = 0, step = 230 / full_event.size();
  for (const auto& point : full_event) {
    const auto center = type::reduce_to_r3(point);
    const plot::color color{brightness, brightness, brightness};
    canvas.add_box(center, point.error.x, point.error.y, point.error.z, 2.5, color);
    canvas.add_point(center, 0.3, color);
    brightness += step;
  }
  canvas.add_line(type::reduce_to_r3(full_event.front()), type::reduce_to_r3(full_event.back()));
  canvas.add_line(track.front(), track.back(), 1, plot::color::RED);
}
//----------------------------------------------------------------------------------------------

//__Add Detector Centers to Canvas______________________________________________________________
void draw_detector_centers(plot::canvas& canvas) {
  for (const auto& name
      : geometry::full_structure_except({"world", "Sandstone", "Marl", "Mix", "Earth"})) {
    canvas.add_point(geometry::limits_of(name).center, 0.25, plot::color::MAGENTA);
  }
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
  for (const auto& path : reader::root::search_directory(options.root_directory)) {
    plot::canvas canvas(path);
    std::cout << path << "\n";
    for (const auto& event : reader::root::import_events(path, options, detector_map)) {

      const auto collapsed_event = analysis::collapse(event, options.collapse_size);
      util::io::print_range(collapsed_event, "\n") << "\n\n";
      canvas.add_points(collapsed_event, 0.8, plot::color::BLUE);
      draw_detector_centers(canvas);

      analysis::full_event seeding_rpc_hits;
      analysis::full_event tracking_rpc_hits;

      const auto full_event_points = combine_rpc_hits(collapsed_event, 2 * units::time, seeding_rpc_hits, tracking_rpc_hits);
      const auto layers = analysis::partition(full_event_points, options.layer_axis, options.layer_depth);
      const auto seeds = analysis::seed(options.seed_size, layers, options.line_width);
      const auto tracking_vector = reset_seeds(analysis::join_all(seeds), seeding_rpc_hits, tracking_rpc_hits);

      std::cout << tracking_vector.size() << "\n";

      for (const auto& seed : tracking_vector) {
        for (const auto& point : seed) {
          std::cout << type::reduce_to_r4(point) << " ";
        }
        std::cout << "\n";
      }
      std::cout << "\n\n";

      for (const auto& track : analysis::fit_seeds(tracking_vector)) {
        draw_track(canvas, track);
        std::cout << track << "\n";
      }
      canvas.draw();
    }
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

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
#include <tracker/util/io.hh>

namespace MATHUSLA {

/*
//__Parse Volume Name to Check if RPC___________________________________________________________
constexpr bool is_rpc(const std::string& name) {
  using namespace TRACKER;
  using namespace util::algorithm;
  const auto size = name.size();
  if (size != 4 && size != 5) {
    return false;
  } else {
    const auto strip = between(name[size - 1], '1', '8');
    const auto pad = (between(name[size - 2], '1', '9') && name[size - 3] == '0')
                   || (name[size - 2] == '0' && name[size - 3] == '1');
    const auto rpc = size == 4 ? between(name[size - 4], '1', '9')
                               : ((between(name[size - 4], '1', '9') && name[size - 5] == '0')
                               || (between(name[size - 4], '0', '2') && name[size - 5] == '1'));
    return strip && pad && rpc;
  }
}
//----------------------------------------------------------------------------------------------

//__Get RPC Volume if Point is Inside RPC_______________________________________________________
const TRACKER::geometry::box_volume get_rpc(const TRACKER::analysis::r4_point& point) {
  using namespace TRACKER;
  const auto name = geometry::volume(point);
  return is_rpc(name) ? geometry::limits_of(name) : {};
}
//----------------------------------------------------------------------------------------------
*/

//__Combine Pair of Hits if they Occur in Overlapping RPCs______________________________________
const TRACKER::geometry::box_volume combine_rpc_hit_pair(const TRACKER::analysis::hit& first,
                                                         const TRACKER::analysis::hit& second) {
  using namespace TRACKER;

  const auto first_limits = geometry::limits_of_volume(first);
  const auto second_limits = geometry::limits_of_volume(second);
  const auto intersection_volume = geometry::intersection_volume(first_limits, second_limits);
  const auto union_volume = geometry::union_volume(first_limits, second_limits);

  std::cout << first_limits << " & " << second_limits << "\n = " << intersection_volume << "\n";
  std::cout << first_limits << " | " << second_limits << "\n = " << union_volume << "\n\n";

  auto out = intersection_volume;
  out.center.z = union_volume.center.z;
  out.min.z = union_volume.min.z;
  out.max.z = union_volume.max.z;
  return out;
}
//----------------------------------------------------------------------------------------------

//__Combine Hits if they Occur in Overlapping RPCs______________________________________________
const TRACKER::geometry::box_volume_vector combine_rpc_hits(const TRACKER::analysis::event& points,
                                                            const TRACKER::analysis::real time_threshold) {
  using namespace TRACKER;
  using namespace util::math;

  static const analysis::real z_threshold = 10 * units::length;

  const auto size = points.size();
  geometry::box_volume_vector boxes;
  boxes.reserve(size);

  size_t i = 0;
  for (; i < size - 1; ++i) {
    const auto current = points[i];
    const auto next = points[i + 1];
    if (within(current.z, next.z, z_threshold) && within(current.t, next.t, time_threshold)) {
      boxes.push_back(combine_rpc_hit_pair(current, next));
      i += 1;
    } else {
      boxes.push_back(geometry::limits_of_volume(current));
    }
  }

  if (i == size)
    boxes.push_back(geometry::limits_of_volume(points[i-1]));

  return boxes;
}

//const TRACKER::analysis::full_event_vector combine_rpc_hits(const TRACKER::analysis::event& points,
//                                                            const TRACKER::analysis::real time_threshold) {
//}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

//__Main Function: Prototype Tracker____________________________________________________________
int main(int argc, char* argv[]) {
  using namespace MATHUSLA;
  using namespace MATHUSLA::TRACKER;

  const auto options = reader::parse_input(argc, argv);
  const auto detector_map = reader::import_detector_map(options.geometry_map_file);

  plot::init();
  geometry::open(options.geometry_file);
  for (const auto& path : reader::root::search_directory(options.root_directory)) {
    for (const auto& event : reader::root::import_events(path, options, detector_map)) {
      plot::canvas canvas(path);

      const auto collapsed_event = analysis::collapse(event, options.collapse_size);

      canvas.add_points(collapsed_event, 1.5, {90, 90, 90});
      for (const auto& name : geometry::full_structure_except({"world", "Sandstone", "Marl", "Mix", "Earth"})) {
        const auto limits = geometry::limits_of(name);
        canvas.add_point(limits.center, 0.25, plot::color::MAGENTA);
      }

      const auto box_volumes = combine_rpc_hits(collapsed_event, 2 * units::time);

      for (const auto& box : box_volumes)
        canvas.add_box(box.min, box.max, 3, plot::color::BLACK);

      /* TODO: finish vvvvv

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

/*
 * demo/prototype/geometry.cc
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

#include "geometry.hh"

#include <tracker/core/stat.hh>
#include <tracker/core/units.hh>
#include <tracker/analysis.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/bit_vector.hh>

//__Namespace Alias_____________________________________________________________________________
namespace stat = MATHUSLA::TRACKER::stat;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

//__Total Geometry of the Prototype Detector____________________________________________________
const geometry::structure_vector prototype_geometry() {
  return geometry::full_structure_except({"World", "Sandstone", "Marl", "Mix", "Earth"});
}
//----------------------------------------------------------------------------------------------

//__Hits per Total Geometry_____________________________________________________________________
type::real modified_geometry_event_density(const analysis::event& event) {
  return event.size() / static_cast<type::real>(prototype_geometry().size());
}
//----------------------------------------------------------------------------------------------

//__Check If Volume is NOT RPC__________________________________________________________________
bool is_not_rpc(const geometry::structure_value& name) {
  return name[0] == 'A' || name[1] == 'B';
}
//----------------------------------------------------------------------------------------------

//__Combine Pair of Hits if they Occur in Overlapping RPCs______________________________________
const geometry::box_volume combine_rpc_volume_pair(const geometry::structure_value& first,
                                                   const geometry::structure_value& second) {
  if (is_not_rpc(first) || is_not_rpc(second))
    return geometry::box_volume{};

  const auto first_box = geometry::limits_of(first);
  const auto second_box = geometry::limits_of(second);
  const auto union_volume = geometry::coordinatewise_union(first_box, second_box);
  auto out = geometry::coordinatewise_intersection(first_box, second_box);
  out.min.z = union_volume.min.z;
  out.center.z = union_volume.center.z;
  out.max.z = union_volume.max.z;
  return out;
}
//----------------------------------------------------------------------------------------------

//__Check If RPC Combine Created a Valid Strip Overlap__________________________________________
bool was_combine_successful(const geometry::box_volume& combined) {
  return (combined.min.x || combined.max.x) && (combined.min.y || combined.max.y);
}
//----------------------------------------------------------------------------------------------

//__Construct True Hit from RPC Hit Volumes_____________________________________________________
const analysis::hit construct_hit(const type::real top_time,
                                  const type::real bottom_time,
                                  const geometry::box_volume& combined) {
  return {0.5L * (top_time + bottom_time), combined.center.x, combined.center.y, combined.center.z};
}
//----------------------------------------------------------------------------------------------

//__Combine All Hits that Occur in Overlapping RPCs_____________________________________________
const analysis::event combine_rpc_hits(const analysis::event& points,
                                       analysis::event& combined_rpc_hits,
                                       analysis::full_event& original_rpc_hits) {

  static const type::real z_lower = 24.0L * units::length;
  static const type::real z_upper = 45.0L * units::length;
  static const type::real time_threshold = 3.0L * units::time;

  using namespace util::math;

  const auto size = points.size();
  if (size == 0UL)
    return analysis::event{};

  analysis::event event;
  event.reserve(size);

  const auto parts = analysis::partition(points, type::Coordinate::Z, z_lower).parts;
  const auto partition_size = parts.size();

  std::size_t layer_index{};
  for (; layer_index < partition_size - 1UL; ++layer_index) {
    const auto top = parts[layer_index];
    const auto bottom = parts[layer_index + 1UL];

    if (within(top.front().z, bottom.back().z, z_lower, z_upper)) {
      const auto top_size = top.size();
      const auto bottom_size = bottom.size();
      util::bit_vector discard_list(bottom_size);

      std::size_t top_index{}, bottom_index{};
      while (top_index < top_size) {
        bottom_index = discard_list.first_unset(bottom_index);

        const auto top_point = top[top_index];
        const auto top_point_name = geometry::volume(top_point);
        if (bottom_index == bottom_size || is_not_rpc(top_point_name)) {
          event.push_back(top_point);
          ++top_index;
          bottom_index = 0UL;
          continue;
        }
        const auto bottom_point = bottom[bottom_index];

        const auto combined = combine_rpc_volume_pair(top_point_name, geometry::volume(bottom_point));
        if (was_combine_successful(combined) && within(top_point.t, bottom_point.t, time_threshold)) {
          const auto new_hit = construct_hit(top_point.t, bottom_point.t, combined);
          event.push_back(new_hit);
          combined_rpc_hits.push_back(new_hit);
          original_rpc_hits.push_back(analysis::add_width(top_point));
          original_rpc_hits.push_back(analysis::add_width(bottom_point));
          discard_list.set(bottom_index);
          ++top_index;
          bottom_index = 0UL;
        } else {
          ++bottom_index;
        }
      }

      for (; top_index < top_size; ++top_index) {
        event.push_back(top[top_index]);
      }
      for (bottom_index = 0UL; bottom_index < bottom_size; ++bottom_index) {
        if (!discard_list[bottom_index])
          event.push_back(bottom[bottom_index]);
      }
      ++layer_index;
    } else {
      util::algorithm::back_insert_transform(top, event,
        [](const auto& part) { return part; });
    }
  }

  if (layer_index == partition_size - 1UL) {
    util::algorithm::back_insert_transform(parts.back(), event,
      [](const auto& part) { return part; });
  }

  // FIXME: because upward going tracks are aligned -Z <-> T
  return util::algorithm::reverse(event);;
}
//----------------------------------------------------------------------------------------------

//__Reset Seed Vector Using RPC Combination Hits________________________________________________
const analysis::full_event_vector reset_seeds(const analysis::event_vector& joined_seeds,
                                              const analysis::event& combined_rpc_hits,
                                              const analysis::full_event& original_rpc_hits) {
  // FIXME: should output be resorted in Time?
  const auto combined_size = combined_rpc_hits.size();
  if (combined_size == 0) {
    analysis::full_event_vector out;
    for (const auto& seed : joined_seeds) {
      analysis::full_event next;
      next.reserve(seed.size());
      for (const auto& hit : seed) {
        next.push_back(analysis::add_width(hit));
      }
      next.shrink_to_fit();
      out.push_back(type::t_sort(next));
    }
    out.shrink_to_fit();
    return out;
  }

  analysis::full_event_vector out;
  out.reserve(joined_seeds.size());
  for (const auto& seed : joined_seeds) {
    analysis::full_event next;
    for (const auto& hit : seed) {
      size_t rpc_index = 0;
      for (; rpc_index < combined_size; ++rpc_index) {
        if (hit == combined_rpc_hits[rpc_index]) {
          next.push_back(original_rpc_hits[2 * rpc_index]);
          next.push_back(original_rpc_hits[2 * rpc_index + 1]);
          break;
        }
      }
      if (rpc_index == combined_size)
        next.push_back(analysis::add_width(hit));
    }
    next.shrink_to_fit();
    out.push_back(type::t_sort(next));
  }
  out.shrink_to_fit();
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

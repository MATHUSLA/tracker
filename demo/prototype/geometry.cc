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

#include <tracker/stat.hh>
#include <tracker/units.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/bit_vector.hh>

//__Namespace Alias_____________________________________________________________________________
namespace stat = MATHUSLA::TRACKER::stat;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

//__Total Geometry of the Prototype Detector____________________________________________________
const geometry::structure_vector prototype_geometry() {
  return geometry::full_structure_except({"world", "Sandstone", "Marl", "Mix", "Earth"});
}
//----------------------------------------------------------------------------------------------

//__Hits per Total Geometry_____________________________________________________________________
type::real modified_geometry_event_density(const analysis::event& event) {
  return event.size() / static_cast<type::real>(prototype_geometry().size());
}
//----------------------------------------------------------------------------------------------

//__Combine Pair of Hits if they Occur in Overlapping RPCs______________________________________
const geometry::box_volume combine_rpc_volume_pair(const geometry::box_volume& first,
                                                   const geometry::box_volume& second) {
  const auto union_volume = geometry::coordinatewise_union(first, second);
  auto out = geometry::coordinatewise_intersection(first, second);
  out.min.z = union_volume.min.z;
  out.center.z = union_volume.center.z;
  out.max.z = union_volume.max.z;
  return out;
}
//----------------------------------------------------------------------------------------------

//__Construct True Hit from RPC Hit Volumes_____________________________________________________
const analysis::full_hit construct_hit(const type::real top_time,
                                       const type::real bottom_time,
                                       const geometry::structure_value& top_volume,
                                       const geometry::structure_value& bottom_volume,
                                       const geometry::box_volume& combined) {
  const type::r4_point errors{
    stat::error::propagate_average(
      geometry::time_resolution_of(top_volume),
      geometry::time_resolution_of(bottom_volume)),
    combined.max.x - combined.min.x,
    combined.max.y - combined.min.y,
    combined.max.z - combined.min.z};
  return {0.5L * (top_time + bottom_time), combined.center.x, combined.center.y, combined.center.z, errors};
}
//----------------------------------------------------------------------------------------------

//__Combine All Hits that Occur in Overlapping RPCs_____________________________________________
const analysis::full_event combine_rpc_hits(const analysis::event& points,
                                            analysis::full_event& combined_rpc_hits,
                                            analysis::full_event& original_rpc_hits) {

  static const type::real z_lower = 24.0L * units::length;
  static const type::real z_upper = 45.0L * units::length;
  static const type::real time_threshold = 2.0L * units::time;

  using namespace util::math;

  const auto size = points.size();
  if (size == 0)
    return {};

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
        if (bottom_index == bottom_size) {
          event.push_back(analysis::add_width(top_point));
          ++top_index;
          bottom_index = 0;
          continue;
        }
        const auto bottom_point = bottom[bottom_index];

        const auto top_volume = geometry::volume(top_point);
        const auto bottom_volume = geometry::volume(bottom_point);
        const auto combined = combine_rpc_volume_pair(
          geometry::limits_of(top_volume),
          geometry::limits_of(bottom_volume));

        if (was_combine_successful(combined) && within(top_point.t, bottom_point.t, time_threshold)) {
          const auto new_hit = construct_hit(top_point.t, bottom_point.t, top_volume, bottom_volume, combined);
          event.push_back(new_hit);
          combined_rpc_hits.push_back(new_hit);
          original_rpc_hits.push_back(analysis::add_width(top_point));
          original_rpc_hits.push_back(analysis::add_width(bottom_point));
          discard_list.set(bottom_index);
          ++top_index;
          bottom_index = 0;
        } else {
          ++bottom_index;
        }
      }

      for (; top_index < top_size; ++top_index) {
        event.push_back(analysis::add_width(top[top_index]));
      }
      for (bottom_index = 0; bottom_index < bottom_size; ++bottom_index) {
        if (!discard_list[bottom_index])
          event.push_back(analysis::add_width(bottom[bottom_index]));
      }
      ++layer_index;
    } else {
      util::algorithm::back_insert_transform(top, event,
        [](const auto& part){ return analysis::add_width(part); });
    }
  }

  if (layer_index == partition_size - 1) {
    util::algorithm::back_insert_transform(parts.back(), event,
      [](const auto& part){ return analysis::add_width(part); });
  }

  util::algorithm::reverse(event);
  return event;
}
//----------------------------------------------------------------------------------------------

//__Reset Seed Vector Using RPC Combination Hits________________________________________________
const analysis::full_event_vector reset_seeds(const analysis::full_event_vector& joined_seeds,
                                              const analysis::full_event& combined_rpc_hits,
                                              const analysis::full_event& original_rpc_hits) {
  const auto combined_size = combined_rpc_hits.size();
  if (combined_size == 0)
    return joined_seeds;

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
        next.push_back(hit);
    }
    next.shrink_to_fit();
    out.push_back(type::t_sort(next));
  }
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

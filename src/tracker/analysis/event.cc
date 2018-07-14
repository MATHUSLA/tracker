/*
 * src/tracker/analysis/event.cc
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

#include <tracker/analysis/event.hh>

#include <queue>

#include <tracker/core/stat.hh>
#include <tracker/util/algorithm.hh>
#include <tracker/geometry.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Calculate Number of Hits per unit Length____________________________________________________
template<class Event,
  typename Point = typename Event::value_type,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
const r4_point event_density(const Event& points) {
  const auto size = points.size();
  if (size == 0UL)
    return r4_point{};
  const auto begin = points.cbegin();
  const auto end = points.cend();
  const auto t_range = std::minmax_element(begin, end, t_ordered<Point>{});
  const auto x_range = std::minmax_element(begin, end, x_ordered<Point>{});
  const auto y_range = std::minmax_element(begin, end, y_ordered<Point>{});
  const auto z_range = std::minmax_element(begin, end, z_ordered<Point>{});
  return r4_point{util::math::abs(t_range.second->t - t_range.first->t),
                  util::math::abs(x_range.second->x - x_range.first->x),
                  util::math::abs(y_range.second->y - y_range.first->y),
                  util::math::abs(z_range.second->z - z_range.first->z)} / size;
}
const r4_point event_density(const event& points) {
  return event_density<>(points);
}
const r4_point event_density(const full_event& points) {
  return event_density<>(points);
}
//----------------------------------------------------------------------------------------------

//__Calculate Number of Hits per Geometric Element______________________________________________
template<class Event,
  typename Point = typename Event::value_type,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
real geometric_event_density(const Event& points) {
  const auto size = points.size();
  if (size == 0)
    return 0;

  geometry::structure_vector structures;
  util::algorithm::back_insert_transform(points, structures,
    [](const auto& point) { return geometry::volume(reduce_to_r4(point)); });

  std::sort(structures.begin(), structures.end());
  const auto begin = structures.begin();

  return size / static_cast<real>(std::distance(begin, std::unique(begin, structures.end())));
}
real geometric_event_density(const event& points) {
  return geometric_event_density<>(points);
}
real geometric_event_density(const full_event& points) {
  return geometric_event_density<>(points);
}
//----------------------------------------------------------------------------------------------

//__Find The Errors Associated with a Hit from Geometry_________________________________________
const full_hit add_width(const hit& point) {
  const auto volume = geometry::volume(point);
  const auto limits = geometry::limits_of(volume);
  const auto center = limits.center;
  const auto min = limits.min;
  const auto max = limits.max;
  return { point.t, center.x, center.y, center.z,
           { geometry::time_resolution_of(volume), max.x - min.x, max.y - min.y, max.z - min.z } };
}
const full_event add_width(const event& points) {
  full_event out;
  out.reserve(points.size());
  util::algorithm::back_insert_transform(points, out, [](const auto& point) { return add_width(point); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Compress Points by R4 Interval______________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event compress(const Event& points,
                     bool time_smearing) {
  const auto size = points.size();
  if (size <= 1)
    return points;

  Event out;
  out.reserve(size);

  const auto sorted_event = t_copy_sort(points);

  using size_type = typename Event::size_type;

  size_type index = 0;
  std::queue<size_type> marked_indices;
  while (index < size) {
    size_type collected = 1, missed_index = 0;
    const auto& point = sorted_event[index];
    const auto volume = geometry::volume(reduce_to_r4(point));
    const auto time_error = geometry::time_resolution_of(volume);
    const auto time_interval = point.t + time_error;
    auto sum = point;

    auto skipped = false;
    while (++index < size) {

      while (!marked_indices.empty() && index++ == marked_indices.front())
        marked_indices.pop();

      const auto& next = sorted_event[index];
      if (next.t > time_interval)
        break;

      if (volume == geometry::volume(reduce_to_r4(next))) {
        ++collected;
        sum += next;
        if (skipped)
          marked_indices.push(index);
      } else if (!skipped) {
        skipped = true;
        missed_index = index;
      }
    }

    if (skipped)
      index = missed_index;

    auto average = sum / collected;

    if (time_smearing) {
      using namespace stat;
      static auto current_error = geometry::default_time_resolution();
      static random::generator gen(random::normal(0, current_error));
      if (current_error != time_error) {
        gen.distribution(random::normal(0, time_error));
        current_error = time_error;
      }
      average.t += gen;
    }

    out.push_back(average);
  }

  return t_sort(out);
}
const event compress(const event& points,
                     bool time_smearing) {
  return compress<>(points, time_smearing);
}
const full_event compress(const full_event& points,
                          bool time_smearing) {
  return compress<>(points, time_smearing);
}
//----------------------------------------------------------------------------------------------

//__Partition Points by Coordinate______________________________________________________________
template<class EventPartition,
  typename Event = typename EventPartition::parts::value_type,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const EventPartition partition(const Event& points,
                               const Coordinate coordinate,
                               const real interval) {
  EventPartition out{{}, coordinate, interval};
  if (points.empty())
    return out;

  auto& parts = out.parts;
  const auto sorted_points = coordinate_stable_copy_sort(coordinate, points);
  const auto size = sorted_points.size();

  typename Event::size_type count = 0;
  auto point_iter = sorted_points.cbegin();
  while (count < size) {
    const auto& point = *point_iter;
    Event current_layer{point};
    ++count;

    while (count < size) {
      const auto& next = *(++point_iter);
      if ((coordinate == Coordinate::T && (next.t > point.t + interval)) ||
          (coordinate == Coordinate::X && (next.x > point.x + interval)) ||
          (coordinate == Coordinate::Y && (next.y > point.y + interval)) ||
          (coordinate == Coordinate::Z && (next.z > point.z + interval))) {
        break;
      }
      current_layer.push_back(next);
      ++count;
    }

    parts.push_back(t_sort(current_layer));
  }

  return out;
}
const event_partition partition(const event& points,
                                const Coordinate coordinate,
                                const real interval) {
  return partition<event_partition>(points, coordinate, interval);
}
const full_event_partition partition(const full_event& points,
                                     const Coordinate coordinate,
                                     const real interval) {
  return partition<full_event_partition>(points, coordinate, interval);
}
//----------------------------------------------------------------------------------------------

//__Reset Partition by new Interval_____________________________________________________________
// TODO: implement move repartition
template<class EventPartition,
  typename Event = typename EventPartition::parts::value_type,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const EventPartition repartition(const EventPartition& previous,
                                 const real interval) {
  return partition(reduce_partition(previous), previous.coordinate, interval);
}
const event_partition repartition(const event_partition& previous,
                                  const real interval) {
  return repartition<event_partition, event>(previous, interval);
}
const full_event_partition repartition(const full_event_partition& previous,
                                       const real interval) {
  return repartition<full_event_partition, full_event>(previous, interval);
}
//----------------------------------------------------------------------------------------------

//__Reduce Event Partition to Events____________________________________________________________
template<class EventPartition,
  typename Event = typename EventPartition::parts::value_type,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event reduce_partition(const EventPartition& previous) {
  const auto& parts = previous.parts;
  Event out;
  out.reserve(std::accumulate(parts.cbegin(), parts.cend(), 0ULL,
    [](const auto size, const auto& event) { return size + event.size(); }));

  for (const auto& event : parts)
    out.insert(out.cend(), event.cbegin(), event.cend());
  return out;
}
const event reduce_partition(const event_partition& previous) {
  return reduce_partition<event_partition, event>(previous);
}
const full_event reduce_partition(const full_event_partition& previous) {
  return reduce_partition<full_event_partition, full_event>(previous);
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

/*
 * src/tracker/analysis.cc
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

#include <numeric>
#include <queue>

#include <tracker/geometry.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/bit_vector.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Full Hit Stream Operator Overload___________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const full_hit& point) {
  return os << "[(" << point.t << ", " << point.x << ", " << point.y << ", " << point.z
            << ") +/- " << point.width << "]";
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

//__Centralize Events by Coordinate_____________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event centralize(const Event& points,
                       const Coordinate coordinate) {
  const auto size = points.size();
  if (size == 0) return {};
  auto out = coordinate_copy_sort(coordinate, points);
  //const hit min{points.front().t, 0, 0, 0};
  //std::transform(out.cbegin(), out.cend(), out.begin(),
  //  [&](const auto& point) { return point - min; });
  return out;
}
const event centralize(const event& points,
                       const Coordinate coordinate) {
  return centralize<>(points, coordinate);
}
const full_event centralize(const full_event& points,
                            const Coordinate coordinate) {
  return centralize<>(points, coordinate);
}
//----------------------------------------------------------------------------------------------

//__Collapse Points by R4 Interval______________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event collapse(const Event& points,
                     const r4_point& ds) {
  const auto size = points.size();
  if (size <= 1) return points;

  Event out;
  out.reserve(size);

  const auto& sorted_event = centralize(points, Coordinate::T);

  using size_type = typename Event::size_type;

  size_type index = 0;
  std::queue<size_type> marked_indicies;
  while (index < size) {
    size_type collected = 1, missed_index = 0;
    const auto& point = sorted_event[index];
    const auto time_interval = point.t + ds.t;
    auto sum = point;

    auto skipped = false;
    while (++index < size) {

      while (!marked_indicies.empty() && index++ == marked_indicies.front())
        marked_indicies.pop();

      const auto& next = sorted_event[index];
      if (next.t > time_interval)
        break;

      if (within_dr(point, next, ds)) {
        ++collected;
        sum += next;
        if (skipped)
          marked_indicies.push(index);
      } else if (!skipped) {
        skipped = true;
        missed_index = index;
      }
    }

    if (skipped)
      index = missed_index;

    out.push_back(sum / collected);
  }

  return out;
}
const event collapse(const event& points,
                     const r4_point& ds) {
  return collapse<>(points, ds);
}
const full_event collapse(const full_event& points,
                          const r4_point& ds) {
  return collapse<>(points, ds);
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
template<class EventPartition,
  typename Event = typename EventPartition::parts::value_type,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const EventPartition repartition(const EventPartition& previous,
                                 const real interval) {
  const auto& parts = previous.parts;
  Event reset_parts;
  reset_parts.reserve(std::accumulate(parts.cbegin(), parts.cend(), 0ULL,
    [](const auto size, const auto& event) { return size + event.size(); }));

  for (const auto& event : parts) {
    reset_parts.insert(reset_parts.cend(), event.cbegin(), event.cend());
  }

  return partition(reset_parts, previous.coordinate, interval);
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

//__Fast Check if Points Form a Line____________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
bool fast_line_check_2d(const Event& points,
                        const real threshold,
                        const Coordinate x1,
                        const Coordinate x2) {
  const auto& line_begin = points.front();
  const auto& line_end = points.back();
  return threshold >= std::accumulate(points.cbegin() + 1, points.cend() - 1, threshold,
    [&](const auto max, const auto& point) {
        return std::max(max, point_line_distance(point, line_begin, line_end, x1, x2)); });
}
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
bool fast_line_check_3d(const Event& points,
                        const real threshold,
                        const Coordinate x1,
                        const Coordinate x2,
                        const Coordinate x3) {
  const auto& line_begin = points.front();
  const auto& line_end = points.back();
  return threshold >= std::accumulate(points.cbegin() + 1, points.cend() - 1, threshold,
    [&](const auto max, const auto& point) {
        return std::max(max, point_line_distance(point, line_begin, line_end, x1, x2, x3)); });
}
bool fast_line_check(const event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2) {
  return fast_line_check_2d<>(points, threshold, x1, x2);
}
bool fast_line_check(const full_event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2) {
  return fast_line_check_2d<>(points, threshold, x1, x2);
}
bool fast_line_check(const event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2,
                     const Coordinate x3) {
  return fast_line_check_3d<>(points, threshold, x1, x2, x3);
}
bool fast_line_check(const full_event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2,
                     const Coordinate x3) {
  return fast_line_check_3d<>(points, threshold, x1, x2, x3);
}
//----------------------------------------------------------------------------------------------

//__Seeding Algorithm___________________________________________________________________________
template<class EventPartition,
  typename EventVector = typename EventPartition::parts,
  typename Event = typename EventVector::value_type,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const EventVector seed(const size_t n,
                       const EventPartition& partition,
                       const real line_threshold) {
  if (n <= 2)
    return {};

  const auto& layers = partition.parts;
  const auto layer_count = layers.size();

  // FIXME: find a close upper bound for out.reserve
  EventVector out{};

  // FIXME: unsure what to do here
  if (layer_count < n)
    return {};

  util::bit_vector_sequence layer_sequence;
  for (const auto& layer : layers)
    layer_sequence.emplace_back(1, layer.size());

  util::order2_permutations(n, layer_sequence, [&](const auto& chooser) {
    Event tuple;
    tuple.reserve(n);

    for (size_t index = 0; index < layer_count; ++index) {
      if (chooser[index]) {
        const auto& layer = layers[index];
        const auto& bits = layer_sequence[index];
        size_t j = 0;
        std::copy_if(layer.cbegin(), layer.cend(), std::back_inserter(tuple),
          [&](auto) { return bits[j++]; });
      }
    }

    if (fast_line_check(t_sort(tuple), line_threshold, Coordinate::X, Coordinate::Y, Coordinate::Z))
      out.push_back(tuple);
  });

  return out;
}
const event_vector seed(const size_t n,
                        const event_partition& partition,
                        const real line_threshold) {
  return seed<event_partition, event_vector>(n, partition, line_threshold);
}
const full_event_vector seed(const size_t n,
                             const full_event_partition& partition,
                             const real line_threshold) {
  return seed<full_event_partition, full_event_vector>(n, partition, line_threshold);
}
//----------------------------------------------------------------------------------------------

//__Check if Seeds can be Joined________________________________________________________________
bool seeds_compatible(const event& first,
                      const event& second,
                      const size_t difference) {
  return std::equal(first.cbegin() + difference, first.cend(), second.cbegin());
}
//----------------------------------------------------------------------------------------------

//__Join Two Seeds in Sequence__________________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event sequential_join(const Event& first,
                            const Event& second,
                            const size_t difference) {
  const auto second_size = second.size();
  const auto overlap = first.size() - difference;

  if (overlap <= 0 || second_size < overlap)
    return {};

  const auto size = difference + second_size;
  Event out;
  out.reserve(size);

  size_t index = 0;
  for (; index < difference; ++index) out.push_back(first[index]);
  for (; index < difference + overlap; ++index) {
    const auto& point = first[index];
    if (point != second[index - difference])
      return {};
    out.push_back(point);
  }
  for (; index < size; ++index) out.push_back(second[index - difference]);

  return out;
}
const event sequential_join(const event& first,
                 const event& second,
                 const size_t difference) {
  return sequential_join<>(first, second, difference);
}
const full_event sequential_join(const full_event& first,
                      const full_event& second,
                      const size_t difference) {
  return sequential_join<>(first, second, difference);
}
//----------------------------------------------------------------------------------------------

//__Join Two Seeds That Overlap_________________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event overlap_join(const Event& first,
                         const Event& second) {
  Event out;
  // TODO: implement
  return out;
}
const event overlap_join(const event& first,
                         const event& second);
const full_event overlap_join(const full_event& first,
                              const full_event& second);
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Index Representation of Event Vector________________________________________________________
using index_vector = std::vector<std::size_t>;
//----------------------------------------------------------------------------------------------

//__Join All Secondaries matching Seed__________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
void _join_secondaries(const size_t seed_index,
                       const size_t difference,
                       EventVector& seed_buffer,
                       const index_vector& indicies,
                       util::bit_vector& join_list,
                       index_vector& out) {
  const auto& seed = seed_buffer[indicies[seed_index]];
  const auto size = indicies.size();
  for (size_t i = 0; i < size; ++i) {
    const auto& next_seed = seed_buffer[indicies[i]];
    const auto& joined_seed = sequential_join(seed, next_seed, difference);
    if (!joined_seed.empty()) {
      seed_buffer.push_back(joined_seed);
      out.push_back(seed_buffer.size() - 1);
      join_list[i] = true;
      join_list[seed_index] = true;
    }
  }
}
//----------------------------------------------------------------------------------------------

//__Queue for Joinable Seeds____________________________________________________________________
using seed_queue = std::queue<index_vector>;
//----------------------------------------------------------------------------------------------

//__Partial Join Seeds from Seed Buffer_________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
bool _partial_join(EventVector& seed_buffer,
                   const index_vector& indicies,
                   const size_t difference,
                   seed_queue& joined,
                   seed_queue& singular,
                   EventVector& out) {
  const auto size = indicies.size();
  if (size <= 1) return false;

  util::bit_vector join_list(size);
  index_vector to_joined, to_singular;
  to_joined.reserve(size);
  to_singular.reserve(size);

  for (size_t index = 0; index < size; ++index)
    _join_secondaries(index, difference, seed_buffer, indicies, join_list, to_joined);

  if (!to_joined.empty()) {
    for (size_t i = 0; i < size; ++i)
      if (!join_list[i])
        to_singular.push_back(indicies[i]);

    joined.push(to_joined);
    singular.push(to_singular);
    return true;
  }

  for (size_t i = 0; i < size; ++i)
    out.push_back(seed_buffer[indicies[i]]);

  return false;
}
//----------------------------------------------------------------------------------------------

//__Join Seeds from Seed Queues_________________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
void _join_next_in_queue(seed_queue& queue,
                         EventVector& seed_buffer,
                         const size_t difference,
                         seed_queue& joined,
                         seed_queue& singular,
                         EventVector& out) {
  if (!queue.empty()) {
    const auto indicies = std::move(queue.front());
    queue.pop();
    _partial_join(seed_buffer, indicies, difference, joined, singular, out);
  }
}
//----------------------------------------------------------------------------------------------

//__Seed Join Loop______________________________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
void _full_join(EventVector& seed_buffer,
                const size_t difference,
                const size_t max_difference,
                seed_queue& joined,
                seed_queue& singular,
                EventVector& out) {
  do {
    _join_next_in_queue(joined, seed_buffer, difference, joined, singular, out);
    _join_next_in_queue(singular, seed_buffer, difference + 1, joined, singular, out);
  } while (!joined.empty() || !singular.empty());
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Seed Join___________________________________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
const EventVector join_all(const EventVector& seeds) {
  const auto size = seeds.size();
  EventVector out;
  out.reserve(size); // FIXME: bad estimate?

  EventVector seed_buffer;
  std::copy(seeds.cbegin(), seeds.cend(), std::back_inserter(seed_buffer));

  seed_queue joined, singular;

  index_vector initial_seeds;
  initial_seeds.reserve(size);
  for (size_t i = 0; i < size; ++i)
    initial_seeds.push_back(i);

  joined.push(initial_seeds);
  _full_join(seed_buffer, 1, 2, joined, singular, out);
  return out;
}
const event_vector join_all(const event_vector& seeds) {
  return join_all<>(seeds);
}
const full_event_vector join_all(const full_event_vector& seeds) {
  return join_all<>(seeds);
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

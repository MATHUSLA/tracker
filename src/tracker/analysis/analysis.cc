/*
 * src/tracker/analysis/analysis.cc
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

#include <tracker/analysis/analysis.hh>

#include <numeric>
#include <queue>

#include <tracker/core/stat.hh>
#include <tracker/util/algorithm.hh>
#include <tracker/util/bit_vector.hh>
#include <tracker/util/index_vector.hh>
#include <tracker/util/type.hh>
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
  const auto& parts = previous.parts;
  Event reset_parts;
  reset_parts.reserve(std::accumulate(parts.cbegin(), parts.cend(), 0ULL,
    [](const auto size, const auto& event) { return size + event.size(); }));

  for (const auto& event : parts)
    reset_parts.insert(reset_parts.cend(), event.cbegin(), event.cend());

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
// TODO: implement different topological types like cylinder (done), cone, parabola, ...etc
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
bool is_linear_2d(const Event& points,
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
bool is_linear_3d(const Event& points,
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
bool is_linear(const event& points,
               const real threshold,
               const Coordinate x1,
               const Coordinate x2) {
  return is_linear_2d<>(points, threshold, x1, x2);
}
bool is_linear(const full_event& points,
               const real threshold,
               const Coordinate x1,
               const Coordinate x2) {
  return is_linear_2d<>(points, threshold, x1, x2);
}
bool is_linear(const event& points,
               const real threshold,
               const Coordinate x1,
               const Coordinate x2,
               const Coordinate x3) {
  return is_linear_3d<>(points, threshold, x1, x2, x3);
}
bool is_linear(const full_event& points,
               const real threshold,
               const Coordinate x1,
               const Coordinate x2,
               const Coordinate x3) {
  return is_linear_3d<>(points, threshold, x1, x2, x3);
}
//----------------------------------------------------------------------------------------------

//__Seeding Algorithm___________________________________________________________________________
template<class EventPartition,
  typename EventVector = typename EventPartition::parts,
  typename Event = typename EventVector::value_type,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const EventVector seed(const std::size_t n,
                       const EventPartition& partition,
                       const real line_threshold) {
  if (n <= 1)
    return EventVector{};

  const auto& layers = partition.parts;
  const auto layer_count = layers.size();

  // FIXME: unsure what to do here
  if (layer_count < n)
    return EventVector{};

  // FIXME: find a close upper bound for out.reserve
  EventVector out{};

  util::bit_vector_sequence layer_sequence;
  for (const auto& layer : layers)
    layer_sequence.emplace_back(1, layer.size());

  util::order2_permutations(n, layer_sequence, [&](const auto& chooser) {
    Event tuple;
    tuple.reserve(n);
    for (size_t index = 0; index < layer_count; ++index)
      if (chooser[index])
        layer_sequence[index].set_conditional_push_back(layers[index], tuple);

    if (is_linear(t_sort(tuple), line_threshold, Coordinate::X, Coordinate::Y, Coordinate::Z))
      out.push_back(tuple);
  });

  return out;
}
const event_vector seed(const std::size_t n,
                        const event_partition& partition,
                        const real line_threshold) {
  return seed<event_partition, event_vector>(n, partition, line_threshold);
}
const full_event_vector seed(const std::size_t n,
                             const full_event_partition& partition,
                             const real line_threshold) {
  return seed<full_event_partition, full_event_vector>(n, partition, line_threshold);
}
//----------------------------------------------------------------------------------------------

//__Join Two Seeds in Sequence__________________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event sequential_join(const Event& first,
                            const Event& second,
                            const std::size_t difference) {
  const auto second_size = second.size();
  const auto overlap = first.size() - difference;

  if (overlap <= 0UL || second_size < overlap)
    return Event{};

  const auto size = difference + second_size;
  Event out;
  out.reserve(size);

  std::size_t index{};
  for (; index < difference; ++index) out.push_back(first[index]);
  for (; index < difference + overlap; ++index) {
    const auto& point = first[index];
    if (point != second[index - difference])
      return Event{};
    out.push_back(point);
  }

  // FIXME: index -= difference instead ?
  for (; index < size; ++index) out.push_back(second[index - difference]);

  return out;
}
const event sequential_join(const event& first,
                            const event& second,
                            const std::size_t difference) {
  return sequential_join<>(first, second, difference);
}
const full_event sequential_join(const full_event& first,
                                 const full_event& second,
                                 const std::size_t difference) {
  return sequential_join<>(first, second, difference);
}
//----------------------------------------------------------------------------------------------

//__Join Two Seeds Such That One is a Subset of the Other_______________________________________
template<class Event,
  typename Point = typename Event::value_type,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
const Event subset_join(const Event& first,
                        const Event& second) {
  if (first.size() >= second.size()) {
    return util::algorithm::range_includes(first, second, t_ordered<Point>{}) ? first : Event{};
  } else {
    return util::algorithm::range_includes(second, first, t_ordered<Point>{}) ? second : Event{};
  }
}
const event subset_join(const event& first,
                        const event& second) {
  return subset_join<>(first, second);
}
const full_event subset_join(const full_event& first,
                             const full_event& second) {
  return subset_join<>(first, second);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Traverse Loop_Join Seeds____________________________________________________________________
template<class Iter, class BackIter, class UnaryFunction>
std::size_t _loop_join_copy_traverse(Iter begin,
                                     Iter end,
                                     BackIter out,
                                     Iter& traversal,
                                     UnaryFunction f) {
  traversal = util::algorithm::copy_until(begin, end, out, f).first;
  return static_cast<std::size_t>(traversal - begin);
}
//----------------------------------------------------------------------------------------------

/* TODO: implement
//__Traverse Loop_Join Seeds____________________________________________________________________
template<class Iter, class BackIter, class UnaryFunction>
std::size_t _loop_join_move_traverse(Iter begin,
                                     Iter end,
                                     BackIter out,
                                     Iter& traversal,
                                     UnaryFunction f) {
  traversal = util::algorithm::move_until(begin, end, out, f).first;
  return static_cast<std::size_t>(traversal - begin);
}
//----------------------------------------------------------------------------------------------
*/

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Join Two Seeds Which form a Loop____________________________________________________________
template<class Event,
  typename Iter  = typename Event::const_iterator,
  typename RIter = typename Event::const_reverse_iterator,
  typename Point = typename Event::value_type,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
const Event loop_join(const Event& first,
                      const Event& second) {
  const auto first_size = first.size();
  const auto second_size = second.size();

  Event out;
  out.reserve(first_size + second_size);
  auto out_inserter = std::back_inserter(out);

  const auto first_begin = first.cbegin();
  const auto first_end = first.cend();
  const auto second_begin = second.cbegin();
  const auto second_end = second.cend();

  Iter first_iter;
  const auto second_front = second.front();
  _loop_join_copy_traverse(first_begin, first_end, out_inserter, first_iter,
    [&](const auto& point) { return point == second_front; });

  if (first_iter == first_end)
    return Event{};

  auto second_iter = second_begin;
  const auto front_overlap = _loop_join_copy_traverse(first_iter, first_end, out_inserter, first_iter,
    [&](const auto& point) { return point != *second_iter++; });

  if (first_iter == first_end)
    return out;

  auto first_stop = first_end - 1;
  auto second_stop = second_end - 1;
  if (static_cast<std::size_t>(first_end - first_iter) >= second_size) {
    const auto second_back = second.back();
    const auto first_rend = first.crend();
    RIter reverse_stop;
    first_stop -= _loop_join_copy_traverse(first.crbegin(), first_rend, out_inserter, reverse_stop,
                    [&](const auto& point) { return point == second_back; });
    first_stop -= _loop_join_copy_traverse(reverse_stop, first_rend, out_inserter, reverse_stop,
                    [&](const auto& point) { return point != *second_stop--; });
    ++second_stop;
  } else {
    const auto first_back = first.back();
    const auto second_rend = second.crend();
    RIter reverse_stop;
    second_stop -= _loop_join_copy_traverse(second.crbegin(), second_rend, out_inserter, reverse_stop,
                     [&](const auto& point) { return point == first_back; });
    second_stop -= _loop_join_copy_traverse(reverse_stop, second_rend, out_inserter, reverse_stop,
                     [&](const auto& point) { return point != *first_stop--; });
    ++first_stop;
  }

  if (++first_stop == first_begin || ++second_stop == second_begin)
    return Event{};

  std::copy(first_iter, first_stop, out_inserter);
  std::copy(second_begin + front_overlap, second_stop, out_inserter);

  util::algorithm::sort_range(out, t_ordered<Point>{});

  const auto out_end = out.end();
  out.erase(std::unique(out.begin(), out_end), out_end);
  out.shrink_to_fit();
  return out;
}
const event loop_join(const event& first,
                      const event& second) {
  return loop_join<>(first, second);
}
const full_event loop_join(const full_event& first,
                           const full_event& second) {
  return loop_join<>(first, second);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Join All Secondaries matching Seed__________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
void _sequential_join_secondaries(const std::size_t seed_index,
                                  const std::size_t difference,
                                  EventVector& seed_buffer,
                                  const util::index_vector<>& indices,
                                  util::bit_vector& join_list,
                                  util::index_vector<>& out) {
  const auto& seed = seed_buffer[indices[seed_index]];
  const auto size = indices.size();
  for (std::size_t index{}; index < size; ++index) {
    const auto& next_seed = seed_buffer[indices[index]];
    const auto joined_seed = sequential_join(seed, next_seed, difference);
    if (!joined_seed.empty()) {
      out.joint_push_back(seed_buffer, joined_seed);
      join_list.set(index);
      join_list.set(seed_index);
    }
  }
}
//----------------------------------------------------------------------------------------------

//__Queue for Joinable Seeds____________________________________________________________________
using seed_queue = std::queue<util::index_vector<>>;
//----------------------------------------------------------------------------------------------

//__Partial Join Seeds from Seed Buffer_________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
bool _sequential_partial_join(EventVector& seed_buffer,
                              const util::index_vector<>& indices,
                              const std::size_t difference,
                              seed_queue& joined,
                              seed_queue& singular,
                              EventVector& out) {
  const auto size = indices.size();

  if (size > 1) {
    util::bit_vector join_list(size);
    util::index_vector<> to_joined, to_singular;
    to_joined.reserve(size);
    to_singular.reserve(size);

    for (std::size_t index{}; index < size; ++index)
      _sequential_join_secondaries(index, difference, seed_buffer, indices, join_list, to_joined);

    if (!to_joined.empty()) {
      join_list.unset_conditional_push_back(indices, to_singular);
      joined.push(to_joined);
      singular.push(to_singular);
      return true;
    }
  }

  indices.conditional_push_back(seed_buffer, out);
  return false;
}
//----------------------------------------------------------------------------------------------

//__Join Seeds from Seed Queues_________________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
void _sequential_join_next_in_queue(seed_queue& queue,
                                    EventVector& seed_buffer,
                                    const std::size_t difference,
                                    seed_queue& joined,
                                    seed_queue& singular,
                                    EventVector& out) {
  if (!queue.empty()) {
    const auto indices = std::move(queue.front());
    queue.pop();
    _sequential_partial_join(seed_buffer, indices, difference, joined, singular, out);
  }
}
//----------------------------------------------------------------------------------------------

//__Seed Join Loop______________________________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
void _sequential_full_join(EventVector& seed_buffer,
                           const std::size_t difference,
                           // TODO: const std::size_t max_difference,
                           seed_queue& joined,
                           seed_queue& singular,
                           EventVector& out) {
  do {
    _sequential_join_next_in_queue(joined, seed_buffer, difference, joined, singular, out);
    _sequential_join_next_in_queue(singular, seed_buffer, difference + 1, joined, singular, out);
  } while (!joined.empty() || !singular.empty());
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Optimally Join All Seeds by Sequence________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
const EventVector sequential_join_all(const EventVector& seeds) {
  const auto size = seeds.size();
  if (size == 1)
    return seeds;

  EventVector out;
  out.reserve(size);

  seed_queue joined, singular;
  joined.emplace(size);

  auto seed_buffer = seeds;
  _sequential_full_join(seed_buffer, 1, joined, singular, out);
  out.shrink_to_fit();
  return out;
}
const event_vector sequential_join_all(const event_vector& seeds) {
  return sequential_join_all<>(seeds);
}
const full_event_vector sequential_join_all(const full_event_vector& seeds) {
  return sequential_join_all<>(seeds);
}
//----------------------------------------------------------------------------------------------

//__Optimally Join All Seeds by Subset__________________________________________________________
template<class EventVector,
  typename Event = typename EventVector::value_type,
  typename Point = typename Event::value_type,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
const EventVector subset_join_all(const EventVector& seeds) {
  const auto size = seeds.size();
  if (size == 1)
    return seeds;

  EventVector out;
  out.reserve(size);

  const auto sorted = util::algorithm::copy_sort_range(seeds, util::type::size_greater<Event>{});
  util::bit_vector joined_list(size);

  size_t top_index = 0, bottom_index = 1;
  while (top_index < size) {
    bottom_index = joined_list.first_unset(bottom_index);
    const auto& top_seed = sorted[top_index];
    if (bottom_index == size) {
      out.push_back(top_seed);
      joined_list.set(top_index);
      top_index = joined_list.first_unset(1 + top_index);
      bottom_index = 1 + top_index;
      continue;
    }
    const auto& bottom_seed = sorted[bottom_index];
    if (util::algorithm::range_includes(top_seed, bottom_seed, t_ordered<Point>{})) {
      joined_list.set(bottom_index);
      // FIXME: check this line vvvvvvvvvvvvvvvvvvvvvvvvv
      top_index = joined_list.first_unset(1 + top_index);
      bottom_index = 1 + top_index;
    } else {
      bottom_index = joined_list.first_unset(1 + bottom_index);
    }
  }

  joined_list.unset_conditional_push_back(sorted, out);
  out.shrink_to_fit();
  return out;
}
const event_vector subset_join_all(const event_vector& seeds) {
  return subset_join_all<>(seeds);
}
const full_event_vector subset_join_all(const full_event_vector& seeds) {
  return subset_join_all<>(seeds);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Check Time then Loop Join___________________________________________________________________
template<class Event,
  typename Point = typename Event::value_type,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
const Event _time_ordered_loop_join(const Event& first,
                                    const Event& second) {
  return first.front().t <= second.front().t
    ? loop_join(first, second)
    : loop_join(second, first);
}
//----------------------------------------------------------------------------------------------
} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Optimally Join All Seeds by Loop____________________________________________________________
template<class EventVector,
  typename Event = typename EventVector::value_type,
  typename Point = typename Event::value_type,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
const EventVector loop_join_all(const EventVector& seeds) {
  const auto size = seeds.size();
  if (size <= 1)
    return seeds;

  EventVector out;
  out.reserve(size);
  auto seed_buffer = util::algorithm::copy_sort_range(seeds, util::type::size_less<Event>{});

  util::bit_vector join_list(size);

  std::size_t top_rindex{}, bottom_rindex = 1UL;
  while (top_rindex < seed_buffer.size()) {
    const auto last_index = seed_buffer.size() - 1;
    const auto top_index = last_index - top_rindex;
    const auto bottom_index = join_list.last_unset(bottom_rindex);
    bottom_rindex = last_index - bottom_index;

    const auto& top_seed = seed_buffer[top_index];
    if (bottom_index == last_index + 1) {
      out.push_back(top_seed);
      join_list.set(top_index);
      top_rindex = last_index - join_list.last_unset(1 + top_rindex);
      bottom_rindex = 1 + top_rindex;
      continue;
    }
    const auto& bottom_seed = seed_buffer[bottom_index];
    const auto joined = _time_ordered_loop_join(top_seed, bottom_seed);
    if (!joined.empty()) {
      if (joined != seed_buffer.back()) {
        seed_buffer.push_back(joined);
        join_list.extend();
      }
      join_list.set(top_index);
      join_list.set(bottom_index);
      top_rindex = 0UL;
    }
    ++bottom_rindex;
  }

  join_list.unset_conditional_push_back(seed_buffer, out);
  out.shrink_to_fit();
  return out;
}
const event_vector loop_join_all(const event_vector& seeds) {
  return loop_join_all<>(seeds);
}
const full_event_vector loop_join_all(const full_event_vector& seeds) {
  return loop_join_all<>(seeds);
}
//----------------------------------------------------------------------------------------------

//__Seed Join___________________________________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
const EventVector join_all(const EventVector& seeds) {
  return loop_join_all(subset_join_all(sequential_join_all(seeds)));
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

/*
 * src/tracker/analysis/seed.cc
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

#include <tracker/analysis/seed.hh>

#include <tracker/util/bit_vector.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

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

//__Check if Points are Monotonic in One Coordinate_____________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
bool is_monotonic(const Event& points,
                  const Coordinate c) {
  if (points.size() <= 2)
    return true;

  if (select_r1(points.front(), c) <= select_r1(points.back(), c)) {
    return std::is_sorted(points.cbegin(), points.cend(),
      [&](const auto& left, const auto& right) { return select_r1(left, c) < select_r1(right, c); });
  } else {
    return std::is_sorted(points.crbegin(), points.crend(),
      [&](const auto& left, const auto& right) { return select_r1(left, c) < select_r1(right, c); });
  }
}
bool is_monotonic(const event& points,
                  const Coordinate c) {
  return is_monotonic<>(points, c);
}
bool is_monotonic(const full_event& points,
                  const Coordinate c) {
  return is_monotonic<>(points, c);
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
    for (std::size_t index{}; index < layer_count; ++index)
      if (chooser[index])
        layer_sequence[index].set_conditional_push_back(layers[index], tuple);

    /*
    if (is_linear(t_sort(tuple), line_threshold, Coordinate::X, Coordinate::Y, Coordinate::Z))
      out.push_back(tuple);
    */

    if (is_monotonic(tuple, Coordinate::T)
        && is_linear(tuple, line_threshold, Coordinate::X, Coordinate::Y, Coordinate::Z)) {
      out.push_back(t_sort(tuple));
    }

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

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

/*
 * include/tracker/analysis/seed.hh
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

#ifndef TRACKER__ANALYSIS__SEED_HH
#define TRACKER__ANALYSIS__SEED_HH
#pragma once

#include <tracker/analysis/event.hh>

#include <tracker/util/bit_vector.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

namespace seed_heuristic { /////////////////////////////////////////////////////////////////////

//__Check if Points are Monotonic in One Coordinate_____________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
bool is_monotonic(const Event& points,
                  const Coordinate c) {
  if (points.size() <= 2UL)
    return true;

  if (select_r1(points.front(), c) <= select_r1(points.back(), c)) {
    return std::is_sorted(points.cbegin(), points.cend(),
      [&](const auto& left, const auto& right) { return select_r1(left, c) < select_r1(right, c); });
  } else {
    return std::is_sorted(points.crbegin(), points.crend(),
      [&](const auto& left, const auto& right) { return select_r1(left, c) < select_r1(right, c); });
  }
}
//----------------------------------------------------------------------------------------------

//__Check if Points are Monotonic in One Coordinate_____________________________________________
template<Coordinate C, class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
bool is_monotonic(const Event& points) {
  if (points.size() <= 2UL)
    return true;

  if (select_r1(points.front(), C) <= select_r1(points.back(), C)) {
    return std::is_sorted(points.cbegin(), points.cend(), coordinate_ordered<C, Event>{});
  } else {
    return std::is_sorted(points.crbegin(), points.crend(), coordinate_ordered<C, Event>{});
  }
}
//----------------------------------------------------------------------------------------------

//__Fast Check if Points Form a Line____________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
bool is_linear(const Event& points,
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
bool is_linear(const Event& points,
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
//----------------------------------------------------------------------------------------------

//__Disjunction Heuristic Type__________________________________________________________________
template<class ...Ts>
struct any {
  template<class Event,
    typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
  bool operator()(const Event&) { return true; }
};
template<class T, class ...Ts>
struct any<T, Ts...> : T, any<Ts...> {
  any(T&& t,
      Ts&& ...ts) : T(std::forward<T>(t)), any<Ts...>(std::forward<Ts>(ts)...) {}
  template<class Event,
    typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
  bool operator()(const Event& points) {
    return T::operator()(points) || any<Ts...>::operator()(points);
  }
};
//----------------------------------------------------------------------------------------------

//__Conjuction Heuristic Type___________________________________________________________________
template<class ...Ts>
struct all {
  template<class Event,
    typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
  bool operator()(const Event&) { return true; }
};
template<class T, class ...Ts>
struct all<T, Ts...> : T, all<Ts...> {
  all(T&& t,
      Ts&& ...ts) : T(std::forward<T>(t)), all<Ts...>(std::forward<Ts>(ts)...) {}
  template<class Event,
    typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
  bool operator()(const Event& points) {
    return T::operator()(points) && all<Ts...>::operator()(points);
  }
};
//----------------------------------------------------------------------------------------------

//__Binary Disjunction Heuristic Type___________________________________________________________
template<class A, class B>
using either = any<A, B>;
//----------------------------------------------------------------------------------------------

//__Binary Conjuction Heuristic Type____________________________________________________________
template<class A, class B>
using both = all<A, B>;
//----------------------------------------------------------------------------------------------

//__Time-Ordered Heuristic Type_________________________________________________________________
struct time_ordered {
  template<class Event,
    typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
  bool operator()(const Event& points) { return is_monotonic<Coordinate::T>(points); }
};
//----------------------------------------------------------------------------------------------

//__Cylinder Heuristic Type_____________________________________________________________________
struct cylinder {
  real radius;
  Coordinate c1, c2, c3;
  cylinder(const real r,
           const Coordinate x1=Coordinate::X,
           const Coordinate x2=Coordinate::Y,
           const Coordinate x3=Coordinate::Z)
      : radius(r), c1(x1), c2(x2), c3(x3) {}

  template<class Event,
    typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
  bool operator()(const Event& points) {
    if (points.size() <= 2UL) return true;
    const auto end = points.cend() - 1;
    for (auto iter = points.cbegin() + 1; iter != end; ++iter)
      if (radius > point_line_distance(*iter, points.front(), points.back(), c1, c2, c3))
        return false;
    return true;
  }
};
//----------------------------------------------------------------------------------------------

//__Double Cone Heuristic Type__________________________________________________________________
struct double_cone {
  real radius;
  Coordinate c1, c2, c3;
  double_cone(const real r,
              const Coordinate x1=Coordinate::X,
              const Coordinate x2=Coordinate::Y,
              const Coordinate x3=Coordinate::Z)
      : radius(r), c1(x1), c2(x2), c3(x3) {}

  template<class Event,
    typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
  bool operator()(const Event& points) {
    if (points.size() <= 2UL) return true;

    const auto front = select_r3(points.front(), c1, c2, c3);
    const auto back = select_r3(points.back(), c1, c2, c3);
    const auto twice_radius = 2.0L * radius;
    const auto scale_factor = twice_radius / std::hypot(twice_radius, norm(back - front));

    const auto end = points.cend() - 1;
    for (auto iter = points.cbegin() + 1; iter != end; ++iter) {
      const auto point = select_r3(*iter, c1, c2, c3);
      if (std::min(norm(point - front),
                   norm(point - back)) * scale_factor
            < point_line_distance(point, front, back, c1, c2, c3))
        return false;
    }
    return true;
  }
};
//----------------------------------------------------------------------------------------------

} /* namespace seed_heuristic */ ///////////////////////////////////////////////////////////////

namespace detail { /////////////////////////////////////////////////////////////////////////////
//__Seeding Algorithm___________________________________________________________________________
template<class SeedHeuristic, class EventPartition,
  typename EventVector = typename EventPartition::parts,
  typename Event = typename EventVector::value_type,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const EventVector seed(const std::size_t n,
                       const EventPartition& partition,
                       SeedHeuristic heuristic) {
  if (n <= 1UL)
    return EventVector{};

  const auto& layers = partition.parts;
  const auto layer_count = layers.size();

  if (layer_count < n)
    return EventVector{};

  // TODO: find a close upper bound for EventVector::reserve
  EventVector out{};

  util::bit_vector_sequence layer_sequence;
  for (const auto& layer : layers)
    layer_sequence.emplace_back(1UL, layer.size());

  util::order2_permutations(n, layer_sequence, [&](const auto& chooser) {
    Event tuple;
    tuple.reserve(n);
    for (std::size_t index{}; index < layer_count; ++index)
      if (chooser[index])
        layer_sequence[index].set_conditional_push_back(layers[index], tuple);

    if (heuristic(tuple))
      out.push_back(t_sort(tuple));
  });

  return out;
}
//----------------------------------------------------------------------------------------------
} /* namespace detail */ ///////////////////////////////////////////////////////////////////////

//__Seeding Algorithm___________________________________________________________________________
template<class SeedHeuristic>
const event_vector seed(const std::size_t n,
                        const event_partition& partition,
                        SeedHeuristic heuristic) {
  return detail::seed<SeedHeuristic, event_partition, event_vector>(n, partition, heuristic);
}
template<class SeedHeuristic>
const full_event_vector seed(const std::size_t n,
                             const full_event_partition& partition,
                             SeedHeuristic heuristic) {
  return detail::seed<SeedHeuristic, full_event_partition, full_event_vector>(n, partition, heuristic);
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__SEED_HH */

/*
 * include/tracker/analysis/monte_carlo.hh
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

#ifndef TRACKER__ANALYSIS__MONTE_CARLO_HH
#define TRACKER__ANALYSIS__MONTE_CARLO_HH
#pragma once

#include <queue>

#include <tracker/analysis/type.hh>
#include <tracker/analysis/track.hh>
// TODO: add features for -> #include <tracker/analysis/vertex.hh>

#include <tracker/core/stat.hh>
#include <tracker/util/algorithm.hh>
#include <tracker/util/bit_vector.hh>
#include <tracker/util/index_vector.hh>
#include <tracker/geometry.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { namespace mc { ////////////////////////////////////////////////////////////

//__Monte-Carlo Truth Event Types_______________________________________________________________
struct hit { std::size_t track_id; real t, x, y, z; };
using event = std::vector<hit>;
using event_vector = std::vector<event>;
struct track { std::size_t track_id; std::vector<analysis::hit> hits; };
using track_vector = std::vector<track>;
//----------------------------------------------------------------------------------------------

// FIXME: implement. Current form insufficient to compare to analysis::track
// class track { ... };

//__Monte-Carlo and Analysis Event Bundle Type__________________________________________________
struct event_bundle { event true_hits; analysis::event hits; };
struct full_event_bundle { event true_hits; analysis::full_event hits; };
using event_bundle_vector = std::vector<event_bundle>;
using full_event_bundle_vector = std::vector<full_event_bundle>;
struct event_vector_bundle { event_vector true_events; analysis::event_vector events; };
struct full_event_vector_bundle { event_vector true_events; analysis::full_event_vector events; };
//----------------------------------------------------------------------------------------------

//__Reduce Event Vector to Event________________________________________________________________
const event reduce(const event_vector& events);
const event_vector reduce(const std::vector<event_vector>& event_vectors);
const event_bundle reduce(const event_bundle_vector& bundles);
const full_event_bundle reduce(const full_event_bundle_vector& bundles);
const event_bundle_vector reduce(const std::vector<event_bundle_vector>& bundle_vectors);
const full_event_bundle_vector reduce(const std::vector<full_event_bundle_vector>& bundle_vectors);
const event_bundle reduce(const event_vector_bundle& bundle);
const full_event_bundle reduce(const full_event_vector_bundle& bundle);
const event_vector_bundle reduce(const std::vector<event_vector_bundle>& bundle_vector);
const full_event_vector_bundle reduce(const std::vector<full_event_vector_bundle>& bundle_vector);
//----------------------------------------------------------------------------------------------

//__Align Bundles By Offsets____________________________________________________________________
void align(event_bundle_vector& bundles,
           const real_vector& offsets,
           const Coordinate coordinate=Coordinate::T);
void align(full_event_bundle_vector& bundles,
           const real_vector& offsets,
           const Coordinate coordinate=Coordinate::T);
void align(event_vector_bundle& bundle,
           const real_vector& offsets,
           const Coordinate coordinate=Coordinate::T);
void align(full_event_vector_bundle& bundle,
           const real_vector& offsets,
           const Coordinate coordinate=Coordinate::T);
void align(std::vector<event_vector_bundle>& bundles,
           const real_vector& offsets,
           const Coordinate coordinate=Coordinate::T);
void align(std::vector<full_event_vector_bundle>& bundles,
           const real_vector& offsets,
           const Coordinate coordinate=Coordinate::T);
//----------------------------------------------------------------------------------------------

//__Type Conversion Helper Functions____________________________________________________________
const track_vector convert_events(const event& points);
// TODO: const event convert(const track& points);
// TODO: const track_vector convert(const event_vector& points);
// TODO: const event_vector convert(const track_vector& points);
//----------------------------------------------------------------------------------------------

//__Monte-Carlo Truth Evaluation________________________________________________________________
class truth_evaluation {
public:
  truth_evaluation(const track_vector& true_tracks,
                   const analysis::track_vector& tracks);

  truth_evaluation(const event& points,
                   const analysis::track_vector& tracks);

  const track_vector truth_tracks() const { return _true_tracks; }
  const analysis::track_vector analysis_tracks() const { return _tracks; }

  integer size_difference() const { return _tracks.size() - _true_tracks.size(); }

private:
  track_vector _true_tracks;
  analysis::track_vector _tracks;
};
//----------------------------------------------------------------------------------------------

//__Monte-Carlo Truth Hit Stream Operator Overload______________________________________________
inline std::ostream& operator<<(std::ostream& os,
                                const hit& point) {
  return os << "[" << point.track_id << " : ("
            << point.t / units::time   << ", "
            << point.x / units::length << ", "
            << point.y / units::length << ", "
            << point.z / units::length << ")]";
}
//----------------------------------------------------------------------------------------------

//__Monte-Carlo Truth Hit Equality______________________________________________________________
constexpr bool operator==(const hit& left,
                          const hit& right) {
  return left.track_id == right.track_id
    && left.t == right.t && left.x == right.x && left.y == right.y && left.z == right.z;
}
//----------------------------------------------------------------------------------------------

//__Monte-Carlo Truth Track Stream Operator Overload____________________________________________
inline std::ostream& operator<<(std::ostream& os,
                                const track& points) {
  os << "[" << points.track_id << " :\n";
  for (const auto& h : points.hits)
    os << "  " << units::scale_r4_length(h) << "\n";
  return os << "]";
}
//----------------------------------------------------------------------------------------------

//__Monte-Carlo Truth Track Equality____________________________________________________________
constexpr bool operator==(const track& left,
                          const track& right) {
  return left.track_id == right.track_id && left.hits == right.hits;
}
//----------------------------------------------------------------------------------------------

namespace detail { /////////////////////////////////////////////////////////////////////////////

//__Compress Points by R4 Interval______________________________________________________________
template<class Geometry,
         class Event,
   typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event compress(const Event& points) {
  const auto size = points.size();
  if (size <= 1UL)
    return points;

  Event out;
  out.reserve(size);

  const auto sorted_event = t_copy_sort(points);

  std::size_t index{};
  std::queue<std::size_t> marked_indices;
  while (index < size) {
    std::size_t collected{1UL}, missed_index{};
    const auto& point = sorted_event[index];
    const auto volume = geometry::custom::volume<Geometry>(reduce_to_r4(point));
    const auto time_error = geometry::custom::time_resolution_of<Geometry>(volume);
    const auto time_interval = point.t + time_error;
    auto sum = point;

    auto skipped = false;
    while (++index < size) {

      while (!marked_indices.empty() && index++ == marked_indices.front())
        marked_indices.pop();

      const auto& next = sorted_event[index];
      if (next.t > time_interval)
        break;

      if (volume == geometry::custom::volume<Geometry>(reduce_to_r4(next))) {
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

    out.push_back(sum / collected);
  }

  return t_sort(out);
}
//----------------------------------------------------------------------------------------------

//__Time Smear Points by Detector Time Resolution_______________________________________________
template<class Geometry,
         class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event time_smear(const Event& points) {
  Event out;
  out.reserve(points.size());

  util::algorithm::back_insert_transform(points, out, [&](auto hit) {
    using namespace stat;
    static auto current_error = geometry::custom::default_time_resolution<Geometry>();
    static random::generator gen(random::normal(0.0L, current_error));

    const auto time_error = geometry::custom::time_resolution_of_volume<Geometry>(reduce_to_r4(hit));
    if (current_error != time_error) {
      gen.distribution(random::normal(0.0L, time_error));
      current_error = time_error;
    }
    hit.t += gen;
    return hit;
  });

  return t_sort(out);
}
//----------------------------------------------------------------------------------------------

//__Simulate the Detector Efficiency____________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event use_efficiency(const Event& points,
                           const real efficiency) {
  if (efficiency >= 1.0L)
    return points;

  Event out;
  out.reserve(points.size());
  util::algorithm::back_insert_copy_if(points, out, [&](const auto) {
    static stat::random::generator gen(stat::random::uniform_real(0.0L, 1.0L));
    return gen <= efficiency;
  });
  out.shrink_to_fit();
  return out;
}
//----------------------------------------------------------------------------------------------

//__Get Poisson Vector__________________________________________________________________________
inline integer_vector _get_poisson_vector(const real mean,
                                          const std::size_t count) {
  static stat::random::generator gen{stat::random::poisson{mean}};
  integer_vector out;
  out.reserve(count);
  std::size_t empty_count{};
  gen.fill_container(count, out, [&](auto& c, auto r) {
    if (r <= 0LL) {
      ++empty_count;
      r = 0LL;
    }
    c.push_back(r);
  });
  return empty_count == count ? integer_vector{} : out;
}
//----------------------------------------------------------------------------------------------

//__Get Ranges Which Share Volumes______________________________________________________________
template<std::size_t N=util::default_index_vector_size>
util::index_vector<N> _get_unique_volume_ranges(const geometry::structure_vector& volumes,
                                                const util::index_vector<N>& indices) {

  const auto size = indices.size();

  auto unique_indices = indices;
  unique_indices.erase(
    std::unique(unique_indices.begin(), unique_indices.end(),
      [&](const auto left, const auto right) { return volumes[left] == volumes[right]; }),
    unique_indices.cend());
  const auto unique_size = unique_indices.size();

  util::index_vector<> out;
  out.reserve(unique_size);
  std::size_t i{}, u{};
  for (; i < size && u < unique_size; ++i) {
    if (indices[i] == unique_indices[u]) {
      out.push_back(i);
      ++u;
    }
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Get Time Width of Point_____________________________________________________________________
template<class Geometry,
         class Point,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
real _3sigma_time_width(const Point& point) {
  return 3.0L * geometry::custom::time_resolution_of_volume<Geometry>(reduce_to_r3(point));
}
//----------------------------------------------------------------------------------------------

//__Check if Time Windows Overlap_______________________________________________________________
template<class Geometry,
         class Point,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
bool _time_window_overlap(const Point& first,
                          const Point& second) {
  return first.t + _3sigma_time_width<Geometry>(first) >= second.t - _3sigma_time_width<Geometry>(second);
}
//----------------------------------------------------------------------------------------------

//__Generate Weights for Distribution___________________________________________________________
inline const real_vector _generate_weights(const real_vector& intervals) {
  const auto size = intervals.size();
  real_vector out;
  out.reserve(size);
  for (std::size_t i{}; i < size - 1UL; ++i)
    out.push_back(!(i % 2));
  return out;
}
//----------------------------------------------------------------------------------------------

//__Get Center of Point From Index______________________________________________________________
template<class Geometry,
         class Event,
  typename Point = typename Event::value_type,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
Point _get_center_from_index(const std::size_t index,
                             const Event& points,
                             const geometry::structure_vector& volumes) {
  auto out = points[index];
  const auto center = geometry::custom::limits_of<Geometry>(volumes[index]).center;
  out.x = center.x;
  out.y = center.y;
  out.z = center.z;
  return out;
}
//----------------------------------------------------------------------------------------------

//__Generate Distribution_______________________________________________________________________
template<class Geometry,
         class Event,
  typename Iter = util::index_vector<>::const_iterator,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
stat::random::piecewise_constant _create_open_time_distribution(const Event& points,
                                                                Iter begin,
                                                                Iter end,
                                                                real begin_time,
                                                                real end_time) {
  const auto size = static_cast<std::size_t>(end - begin);
  real_vector intervals;
  intervals.reserve(size);
  intervals.push_back(begin_time);

  util::bit_vector join_list{size};

  std::size_t top_index{}, bottom_index{1UL};
  while (top_index < size) {
    const auto top = points[begin[top_index]];
    real active_min{top.t - _3sigma_time_width<Geometry>(top)},
         active_max{top.t + _3sigma_time_width<Geometry>(top)};

    bottom_index = join_list.first_unset(1UL + top_index);

    while (bottom_index < size) {
      const auto bottom = points[begin[bottom_index]];
      if (_time_window_overlap<Geometry>(top, bottom)) {
        const auto width = _3sigma_time_width<Geometry>(bottom);
        active_min = std::min(active_min, bottom.t - width);
        active_max = std::max(active_max, bottom.t + width);
        join_list.set(bottom_index);
      }
      bottom_index = join_list.first_unset(1UL + bottom_index);
    }

    join_list.set(top_index);
    top_index = join_list.first_unset(1UL + top_index);

    intervals.push_back(std::max(begin_time, active_min));
    intervals.push_back(std::min(end_time, active_max));
  }
  intervals.push_back(end_time);

  util::algorithm::sort_range(intervals);
  intervals.erase(std::unique(intervals.begin(), intervals.end()), intervals.cend());
  intervals.shrink_to_fit();
  return stat::random::piecewise_constant{intervals, _generate_weights(intervals)};
}
//----------------------------------------------------------------------------------------------

//__Merge Events by Move into Output____________________________________________________________
template<class Event,
  typename CIter = typename Event::const_iterator,
  typename Iter = typename Event::iterator,
  typename Point = typename Event::value_type,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
const Event _move_merge_events(CIter points_begin,
                               CIter points_end,
                               Iter move_begin,
                               Iter move_end) {
  Event out;
  out.reserve(points_end - points_begin + move_end - move_begin);
  std::merge(points_begin,
             points_end,
             std::make_move_iterator(move_begin),
             std::make_move_iterator(move_end),
             std::back_inserter(out),
             t_ordered<Point>{});
  return out;
}
//----------------------------------------------------------------------------------------------

//__Add Noise to Points_________________________________________________________________________
template<class Geometry,
         class Event,
  typename Point = typename Event::value_type,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
const Event add_noise(const Event& points,
                      const real rate,
                      const real begin_time,
                      const real end_time) {
  if (rate == 0.0L)
    return points;

  const auto time_width = end_time - begin_time;
  const auto size = points.size();

  const auto add_counts = _get_poisson_vector(size, rate * time_width);
  if (add_counts.empty())
    return points;

  geometry::structure_vector volumes;
  volumes.reserve(size);
  util::algorithm::back_insert_transform(points, volumes,
    [](const auto& point) { return geometry::custom::volume<Geometry>(reduce_to_r3(point)); });

  util::index_vector<> indices{size};
  util::algorithm::sort_range(indices,
    [&](const auto left, const auto right) { return volumes[left] < volumes[right]; });
  const auto indices_begin = indices.cbegin();
  const auto indices_end = indices.cend();

  auto unique_ranges = _get_unique_volume_ranges<>(volumes, indices);
  const auto unique_size = unique_ranges.size();

  Event added_points;
  added_points.reserve(size);

  static stat::random::generator gen;
  for (std::size_t i{}; i < unique_size; ++i) {
    const std::size_t add_count = add_counts[i];
    if (add_count == 0UL)
      continue;

    const auto range_begin = indices_begin + unique_ranges[i];
    const auto range_end = i + 1UL == unique_size ? indices_end : indices_begin + unique_ranges[i + 1UL];

    gen.distribution(_create_open_time_distribution<Geometry, Event>(points, range_begin, range_end, begin_time, end_time));

    auto new_point = _get_center_from_index<Geometry>(indices[unique_ranges[i]], points, volumes);
    for (std::size_t j{}; j < add_count; ++j) {
      new_point.t = gen;
      added_points.push_back(new_point);
    }
  }

  return _move_merge_events<Event>(points.cbegin(), points.cend(),
                                   added_points.begin(), added_points.end());
}
//----------------------------------------------------------------------------------------------

} /* namespace detail */ ///////////////////////////////////////////////////////////////////////

//__Compress Points by R4 Interval______________________________________________________________
template<class Geometry=void>
const analysis::event compress(const analysis::event& points) {
  return detail::compress<Geometry, analysis::event>(points);
}
template<class Geometry=void>
const analysis::full_event compress(const analysis::full_event& points) {
  return detail::compress<Geometry, analysis::full_event>(points);
}
//----------------------------------------------------------------------------------------------

//__Time Smear Points by Detector Time Resolution_______________________________________________
template<class Geometry=void>
const analysis::event time_smear(const analysis::event& points) {
  return detail::time_smear<Geometry, analysis::event>(points);
}
template<class Geometry=void>
const analysis::full_event time_smear(const analysis::full_event& points) {
  return detail::time_smear<Geometry, analysis::full_event>(points);
}
//----------------------------------------------------------------------------------------------

//__Simulate the Detector Efficiency____________________________________________________________
inline const analysis::event use_efficiency(const analysis::event& points,
                                            const real efficiency) {
  return detail::use_efficiency<analysis::event>(points, efficiency);
}
inline const analysis::full_event use_efficiency(const analysis::full_event& points,
                                                 const real efficiency) {
  return detail::use_efficiency<analysis::full_event>(points, efficiency);
}
//----------------------------------------------------------------------------------------------

//__Add Noise to Points_________________________________________________________________________
template<class Geometry=void>
const analysis::event add_noise(const analysis::event& points,
                                const real rate,
                                const real begin_time,
                                const real end_time) {
  return detail::add_noise<Geometry, analysis::event>(points, rate, begin_time, end_time);
}
template<class Geometry=void>
const analysis::full_event add_noise(const analysis::full_event& points,
                                     const real rate,
                                     const real begin_time,
                                     const real end_time) {
  return detail::add_noise<Geometry, analysis::full_event>(points, rate, begin_time, end_time);
}
//----------------------------------------------------------------------------------------------

} } /* namespace analysis::mc */ ///////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__MONTE_CARLO_HH */

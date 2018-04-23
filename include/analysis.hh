#ifndef TRACKER__ANALYSIS_HH
#define TRACKER__ANALYSIS_HH
#pragma once

#include <cmath>
#include <deque>

#include "types.hh"

namespace MATHUSLA { namespace TRACKER {

namespace Analysis { ///////////////////////////////////////////////////////////////////////////

//__ROOT Directory Search_______________________________________________________________________
std::vector<std::string> search_directory(const std::string& path);
//----------------------------------------------------------------------------------------------

//__Event Types_________________________________________________________________________________
using event_points = r4_point_vector;
using event_vector = std::vector<event_points>;

template<std::size_t N> using event_tuple = r4_point_array<N>;
using event_pair      = event_tuple<2>;
using event_triple    = event_tuple<3>;
using event_quadruple = event_tuple<4>;
template<std::size_t N> using event_tuple_vector = std::vector<event_tuple<N>>;

struct event_partition { event_vector parts; Coordinate coordinate; };
//----------------------------------------------------------------------------------------------

//__Event File Import___________________________________________________________________________
event_vector import_events(const std::string& path,
                           const std::array<std::string, 4>& keys={{"T", "X", "Y", "Z"}});
//----------------------------------------------------------------------------------------------

//__R3 Interval Check___________________________________________________________________________
inline bool within_dr(const r4_point& a,
                      const r4_point& b,
                      const r4_point& dr) {
  return (std::abs(a.x - b.x) <= dr.x) && (std::abs(a.y - b.y) <= dr.y) && (std::abs(a.z - b.z) <= dr.z);
}
//----------------------------------------------------------------------------------------------

//__R4 Interval Check___________________________________________________________________________
inline bool within_ds(const r4_point& a,
                      const r4_point& b,
                      const r4_point& ds) {
  return (std::abs(a.t - b.t) <= ds.t) && within_dr(a, b, ds);
}
//----------------------------------------------------------------------------------------------

//__Collapse Points by R4 Interval______________________________________________________________
inline event_points collapse(const event_points& event,
                             const r4_point& ds) {
  event_points out{};
  if (event.empty())
    return event;

  const auto& sorted_events = t_copy_sort(event);
  const auto&& event_size = sorted_events.size();

  event_points::size_type index = 0;
  std::deque<decltype(index)> marked_indicies{};
  while (index < event_size) {
    event_points::size_type collected = 1, missed_index = 0;
    const auto& point = sorted_events[index];
    const auto&& time_interval = point.t + ds.t;
    auto sum_t = point.t,
         sum_x = point.x,
         sum_y = point.y,
         sum_z = point.z;

    auto skipped = false;
    while (++index < event_size) {
      while (!marked_indicies.empty() && index == marked_indicies.front()) {
        ++index;
        marked_indicies.pop_front();
      }

      const auto& next = sorted_events[index];
      if (next.t > time_interval)
        break;

      if (within_dr(point, next, ds)) {
        ++collected;
        sum_t += next.t;
        sum_x += next.x;
        sum_y += next.y;
        sum_z += next.z;
        if (skipped)
          marked_indicies.push_back(index);
      } else if (!skipped) {
        skipped = true;
        missed_index = index;
      }
    }

    if (skipped)
      index = missed_index;

    out.push_back({sum_t / collected, sum_x / collected, sum_y / collected, sum_z / collected});
  }

  return out;
}
//----------------------------------------------------------------------------------------------

//__Partition Points by Coordinate______________________________________________________________
template<Coordinate coordinate=Coordinate::Z>
inline event_partition partition(const event_points& points,
                                 const real interval) {
  event_partition out{{}, coordinate};
  if (points.empty())
    return out;

  auto& parts = out.parts;

  const auto& sorted_points = coordinate_copy_sort(points, coordinate);
  const auto&& size = sorted_points.size();

  event_points::size_type count = 0;
  auto point_iter = sorted_points.cbegin();
  while (count < size) {
    parts.emplace_back();
    auto& current_layer = parts.back();
    const auto& point = *point_iter;
    current_layer.push_back(point);
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

    t_stable_sort(current_layer);
  }

  return out;
}
//----------------------------------------------------------------------------------------------

//__SeedN Algorithm_____________________________________________________________________________
template<std::size_t N>
inline event_tuple_vector<N> seed(const event_points& event,
                                  const r4_point& ds,
                                  const real layer_dz) {
  const auto& points = collapse(event, ds);
  const auto&& size = points.size();

  if (size <= N)
    return {to_array<N>(points)};

  const auto& layers = partition<>(points, layer_dz);

  event_tuple_vector<N> out{};
  out.reserve(std::pow(size, N) / std::pow(N/2.718, N));  // work on this limit

  //TODO: finish seeding algorithm

  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace Analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS_HH */

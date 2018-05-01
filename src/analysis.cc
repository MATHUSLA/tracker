#include "analysis.hh"

#include <deque>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Collapse Points by R4 Interval______________________________________________________________
event_points collapse(const event_points& event,
                      const r4_point& ds) {
  if (event.empty())
    return event;

  event_points out;

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
event_partition partition(const event_points& points,
                          const real interval,
                          const Coordinate coordinate) {
  event_partition out{{}, coordinate};
  if (points.empty())
    return out;

  auto& parts = out.parts;

  const auto& sorted_points = coordinate_copy_sort(points, coordinate);
  const auto&& size = sorted_points.size();

  event_points::size_type count = 0;
  auto point_iter = sorted_points.cbegin();
  while (count < size) {
    const auto& point = *point_iter;
    event_points current_layer{point};
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
    parts.push_back(current_layer);
  }

  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

/*
 * src/tracker/analysis/simulation.cc
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

#include <tracker/analysis/simulation.hh>

#include <queue>

#include <tracker/core/stat.hh>
#include <tracker/geometry.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { namespace simulation { ////////////////////////////////////////////////////

//__Compress Points by R4 Interval______________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event compress(const Event& points) {
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

    // TODO: move to time_smear
    if (true) {
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
const event compress(const event& points) {
  return compress<>(points);
}
const full_event compress(const full_event& points) {
  return compress<>(points);
}
//----------------------------------------------------------------------------------------------

//__Time Smear Points by Detector Time Resolution_______________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event time_smear(const Event& points) {
  // TODO: implement
  return points;
}
const event time_smear(const event& points) {
  return time_smear<>(points);
}
const full_event time_smear(const full_event& points) {
  return time_smear<>(points);
}
//----------------------------------------------------------------------------------------------

//__Simulate the Detector Efficiency____________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event use_efficiency(const Event& points,
                           const real efficiency) {
  // TODO: implement
  return points;
}
const event use_efficiency(const event& points,
                           const real efficiency) {
  return use_efficiency<>(points, efficiency);
}
const full_event use_efficiency(const full_event& points,
                                const real efficiency) {
  return use_efficiency<>(points, efficiency);
}
//----------------------------------------------------------------------------------------------

//__Add Noise to Points_________________________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event add_noise(const Event& points) {
  // TODO: implement
  return points;
}
const event add_noise(const event& points) {
  return add_noise<>(points);
}
const full_event add_noise(const full_event& points) {
  return add_noise<>(points);
}
//----------------------------------------------------------------------------------------------

} } /* namespace analysis::simulation */ ///////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

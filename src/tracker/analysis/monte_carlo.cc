/*
 * src/tracker/analysis/monte_carlo.cc
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

#include <tracker/analysis/monte_carlo.hh>

#include <tracker/util/algorithm.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { namespace mc { ////////////////////////////////////////////////////////////

//__Type Conversion Helper Functions____________________________________________________________
const track_vector convert_events(const event& points) {
  const auto size = points.size();
  if (size == 0)
    return track_vector{};

  track_vector out;
  out.reserve(size);
  const auto sorted = util::algorithm::copy_sort_range(points,
    [](const auto& left, const auto& right) { return left.track_id < right.track_id; });

  const auto& first = sorted.front();
  track next{first.track_id, {reduce_to_r4(first)}};
  next.hits.reserve(size);
  for (std::size_t i = 1; i < size; ++i) {
    const auto& point = sorted[i];
    if (next.track_id == point.track_id) {
      next.hits.push_back(reduce_to_r4(point));
      if (i == size - 1) {
        next.hits.shrink_to_fit();
        out.push_back(next);
      }
    } else {
      next.hits.shrink_to_fit();
      out.push_back(next);
      next = {point.track_id, {reduce_to_r4(point)}};
      next.hits.reserve(size);
    }
  }

  out.shrink_to_fit();
  return out;
}
//----------------------------------------------------------------------------------------------

//__Truth Evaluation Constructor________________________________________________________________
truth_evaluation::truth_evaluation(const track_vector& true_tracks,
                                   const analysis::track_vector& tracks)
    : _true_tracks(true_tracks), _tracks(tracks) {}
//----------------------------------------------------------------------------------------------

//__Truth Evaluation Constructor________________________________________________________________
truth_evaluation::truth_evaluation(const event& points,
                                   const analysis::track_vector& tracks)
    : truth_evaluation(convert_events(points), tracks) {}
//----------------------------------------------------------------------------------------------

} } /* namespace analysis::mc */ ///////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

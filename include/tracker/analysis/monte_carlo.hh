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

#include <tracker/analysis/type.hh>
#include <tracker/analysis/track.hh>
// #include <tracker/analysis/vertex.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { namespace mc { ////////////////////////////////////////////////////////////

//__Monte-Carlo Truth Event Types_______________________________________________________________
struct hit { std::size_t track_id; real t, x, y, z; };
using event = std::vector<hit>;
using event_vector = std::vector<event>;
struct track { std::size_t track_id; std::vector<analysis::hit> hits; };
using track_vector = std::vector<track>;
//----------------------------------------------------------------------------------------------

//__Monte-Carlo and Analysis Event Bundle Type__________________________________________________
struct event_bundle { event true_hits; analysis::event hits; };
struct full_event_bundle { event true_hits; analysis::full_event hits; };
using event_bundle_vector = std::vector<event_bundle>;
using full_event_bundle_vector = std::vector<full_event_bundle>;
struct event_vector_bundle { event_vector true_events; analysis::event_vector events; };
struct full_event_vector_bundle { event_vector true_events; analysis::full_event_vector events; };
//----------------------------------------------------------------------------------------------

//__Type Conversion Helper Functions____________________________________________________________
const track convert(const event& points);
const event convert(const track& points);
const track_vector convert(const event_vector& points);
const event_vector convert(const track_vector& points);
//----------------------------------------------------------------------------------------------

//__Monte-Carlo Truth Evaluation________________________________________________________________
class truth_evaluation {
public:
  truth_evaluation(const track_vector& true_tracks,
                   const analysis::track_vector& tracks);

  const track_vector truth_tracks() const { return _true_tracks; }
  const analysis::track_vector analysis_tracks() const { return _tracks; }

private:
  track_vector _true_tracks;
  analysis::track_vector _tracks;
};
//----------------------------------------------------------------------------------------------

//__Monte-Carlo Truth Hit Stream Operator Overload______________________________________________
inline std::ostream& operator<<(std::ostream& os,
                                const hit& point) {
  return os << "[" << point.track_id << " : ("
            << point.t << ", " << point.x << ", " << point.y << ", " << point.z << ")]";
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
    os << h << "\n";
  return os << "]";
}
//----------------------------------------------------------------------------------------------

//__Monte-Carlo Truth Track Equality____________________________________________________________
constexpr bool operator==(const track& left,
                          const track& right) {
  return left.track_id == right.track_id && left.hits == right.hits;
}
//----------------------------------------------------------------------------------------------

} } /* namespace analysis::mc */ ///////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__MONTE_CARLO_HH */

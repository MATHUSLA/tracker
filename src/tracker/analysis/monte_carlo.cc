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

//__Reduce Event Vector to Event________________________________________________________________
const event reduce(const event_vector& events) {
  event out;
  out.reserve(std::accumulate(events.cbegin(), events.cend(), 0ULL,
    [](const auto size, const auto& hits) { return size + hits.size(); }));
  for (const auto& hits : events)
    out.insert(out.cend(), hits.cbegin(), hits.cend());
  return out;
}
template<class EventBundleVector,
  typename EventBundle = typename EventBundleVector::value_type>
const EventBundle reduce_vector(const EventBundleVector& bundles) {
  EventBundle out;
  auto& true_hits = out.true_hits;
  auto& hits = out.hits;
  true_hits.reserve(std::accumulate(bundles.cbegin(), bundles.cend(), 0ULL,
    [](const auto size, const auto& bundle) { return size + bundle.true_hits.size(); }));
  hits.reserve(std::accumulate(bundles.cbegin(), bundles.cend(), 0ULL,
    [](const auto size, const auto& bundle) { return size + bundle.hits.size(); }));
  for (const auto& bundle : bundles) {
    true_hits.insert(true_hits.cend(), bundle.true_hits.cbegin(), bundle.true_hits.cend());
    hits.insert(hits.cend(), bundle.hits.cbegin(), bundle.hits.cend());
  }
  return out;
}
const event_bundle reduce(const event_bundle_vector& bundles) {
  return reduce_vector<>(bundles);
}
const full_event_bundle reduce(const full_event_bundle_vector& bundles) {
  return reduce_vector<>(bundles);
}
template<class EventVectorBundle, class EventBundle>
const EventBundle reduce_bundle(const EventVectorBundle& bundle) {
  EventBundle out;
  auto& true_hits = out.true_hits;
  auto& hits = out.hits;
  true_hits.reserve(std::accumulate(bundle.true_events.cbegin(), bundle.true_events.cend(), 0ULL,
    [](const auto size, const auto& event) { return size + event.size(); }));
  hits.reserve(std::accumulate(bundle.events.cbegin(), bundle.events.cend(), 0ULL,
    [](const auto size, const auto& event) { return size + event.size(); }));
  for (const auto& event : bundle.true_events)
    true_hits.insert(true_hits.cend(), event.cbegin(), event.cend());
  for (const auto& event : bundle.events)
    hits.insert(hits.cend(), event.cbegin(), event.cend());
  return out;
}
const event_bundle reduce(const event_vector_bundle& bundle) {
  return reduce_bundle<event_vector_bundle, event_bundle>(bundle);
}
const full_event_bundle reduce(const full_event_vector_bundle& bundle) {
  return reduce_bundle<full_event_vector_bundle, full_event_bundle>(bundle);
}
//----------------------------------------------------------------------------------------------

//__Align Bundles By Offsets____________________________________________________________________
template<class EventBundleVector>
void align(EventBundleVector& bundles,
           const real_vector& offsets,
           const Coordinate coordinate) {
  const auto size = bundles.size();
  for (std::size_t i{}; i < size; ++i) {
    auto& bundle = bundles[i];
    for (auto& true_point : bundle.true_hits)
      select(true_point, coordinate) += offsets[i];
    for (auto& point : bundle.hits)
      select(point, coordinate) += offsets[i];
  }
}
void align(event_bundle_vector& bundles,
           const real_vector& offsets,
           const Coordinate coordinate) {
  align<>(bundles, offsets, coordinate);
}
void align(full_event_bundle_vector& bundles,
           const real_vector& offsets,
           const Coordinate coordinate) {
  align<>(bundles, offsets, coordinate);
}
template<class EventVectorBundle, class EventBundle>
void align(EventVectorBundle& bundle,
           const real_vector& offsets,
           const Coordinate coordinate) {
  const auto size = bundle.true_events.size();
  for (std::size_t i{}; i < size; ++i) {
    for (auto& true_point : bundle.true_events[i])
      select(true_point, coordinate) += offsets[i];
    for (auto& point : bundle.events[i])
      select(point, coordinate) += offsets[i];
  }
}
void align(event_vector_bundle& bundle,
           const real_vector& offsets,
           const Coordinate coordinate) {
  align<event_vector_bundle, event_bundle>(bundle, offsets, coordinate);
}
void align(full_event_vector_bundle& bundle,
           const real_vector& offsets,
           const Coordinate coordinate) {
  align<full_event_vector_bundle, full_event_bundle>(bundle, offsets, coordinate);
}
template<class EventVectorBundleVector, class EventVectorBundle, class EventBundle>
void align(EventVectorBundleVector& bundles,
           const real_vector& offsets,
           const Coordinate coordinate) {
  const auto size = bundles.size();
  for (std::size_t i{}; i < size; ++i) {
    auto& bundle = bundles[i];
    const auto bundle_size = bundle.true_events.size();
    for (std::size_t j{}; j < bundle_size; ++j) {
      for (auto& true_point : bundle.true_events[j])
        select(true_point, coordinate) += offsets[i];
      for (auto& point : bundle.events[j])
        select(point, coordinate) += offsets[i];
    }
  }
}
void align(std::vector<event_vector_bundle>& bundles,
           const real_vector& offsets,
           const Coordinate coordinate) {
  align<decltype(bundles), event_vector_bundle, event_bundle>(bundles, offsets, coordinate);
}
void align(std::vector<full_event_vector_bundle>& bundles,
           const real_vector& offsets,
           const Coordinate coordinate) {
  align<decltype(bundles), full_event_vector_bundle, full_event_bundle>(bundles, offsets, coordinate);
}
//----------------------------------------------------------------------------------------------

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

/*
 * src/tracker/analysis/join.cc
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

#include <tracker/analysis/join.hh>

#include <queue>

#include <tracker/util/algorithm.hh>
#include <tracker/util/bit_vector.hh>
#include <tracker/util/index_vector.hh>
#include <tracker/util/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

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

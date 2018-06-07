/*
 * src/analysis.cc
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

#include <tracker/analysis.hh>

#include <numeric>
#include <queue>

#include <ROOT/TMinuit.h>

#include <tracker/geometry.hh>
#include <tracker/units.hh>

#include <tracker/util/bit_vector.hh>
#include <tracker/util/error.hh>
#include <tracker/util/io.hh>
#include <tracker/util/math.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Find The Errors Associated with a Hit from Geometry_________________________________________
const full_hit add_errors(const hit& point) {
  const auto volume = geometry::volume(point);
  const auto limits = geometry::limits_of(volume);
  const auto& center = limits.center;
  const auto& min = limits.min;
  const auto& max = limits.max;
  return { point.t, center.x, center.y, center.z,
           { geometry::time_resolution_of(volume),
             limits.max.x - limits.min.x,
             limits.max.y - limits.min.y,
             limits.max.z - limits.min.z } };
}
const full_event add_errors(const event& points) {
  full_event out;
  out.reserve(points.size());
  std::transform(points.cbegin(), points.cend(), std::back_inserter(out),
    [](const auto& point) { return add_errors(point); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Centralize Events by Coordinate_____________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event centralize(const Event& points,
                       const Coordinate coordinate) {
  const auto size = points.size();
  if (size == 0) return {};
  auto out = coordinate_copy_sort(coordinate, points);
  //const hit min{points.front().t, 0, 0, 0};
  //std::transform(out.cbegin(), out.cend(), out.begin(),
  //  [&](const auto& point) { return point - min; });
  return out;
}
const event centralize(const event& points,
                       const Coordinate coordinate) {
  return centralize<>(points, coordinate);
}
const full_event centralize(const full_event& points,
                            const Coordinate coordinate) {
  return centralize<>(points, coordinate);
}
//----------------------------------------------------------------------------------------------

//__Collapse Points by R4 Interval______________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event collapse(const Event& points,
                     const r4_point& ds) {
  const auto size = points.size();
  if (size <= 1) return points;

  Event out;
  out.reserve(size);

  const auto& sorted_event = centralize(points, Coordinate::T);

  using size_type = typename Event::size_type;

  size_type index = 0;
  std::queue<size_type> marked_indicies;
  while (index < size) {
    size_type collected = 1, missed_index = 0;
    const auto& point = sorted_event[index];
    const auto time_interval = point.t + ds.t;
    auto sum = point;

    auto skipped = false;
    while (++index < size) {

      while (!marked_indicies.empty() && index++ == marked_indicies.front())
        marked_indicies.pop();

      const auto& next = sorted_event[index];
      if (next.t > time_interval)
        break;

      if (within_dr(point, next, ds)) {
        ++collected;
        sum += next;
        if (skipped)
          marked_indicies.push(index);
      } else if (!skipped) {
        skipped = true;
        missed_index = index;
      }
    }

    if (skipped)
      index = missed_index;

    out.push_back(sum / collected);
  }

  return out;
}
const event collapse(const event& points,
                     const r4_point& ds) {
  return collapse<>(points, ds);
}
const full_event collapse(const full_event& points,
                          const r4_point& ds) {
  return collapse<>(points, ds);
}
//----------------------------------------------------------------------------------------------

//__Partition Points by Coordinate______________________________________________________________
template<class EventPartition,
  typename Event = typename EventPartition::parts::value_type,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const EventPartition partition(const Event& points,
                               const Coordinate coordinate,
                               const real interval) {
  EventPartition out{{}, coordinate, interval};
  if (points.empty())
    return out;

  auto& parts = out.parts;
  const auto sorted_points = coordinate_stable_copy_sort(coordinate, points);
  const auto size = sorted_points.size();

  typename Event::size_type count = 0;
  auto point_iter = sorted_points.cbegin();
  while (count < size) {
    const auto& point = *point_iter;
    Event current_layer{point};
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

    parts.push_back(t_sort(current_layer));
  }

  return out;
}
const event_partition partition(const event& points,
                                const Coordinate coordinate,
                                const real interval) {
  return partition<event_partition>(points, coordinate, interval);
}
const full_event_partition partition(const full_event& points,
                                     const Coordinate coordinate,
                                     const real interval) {
  return partition<full_event_partition>(points, coordinate, interval);
}
//----------------------------------------------------------------------------------------------

//__Reset Partition by new Interval_____________________________________________________________
template<class EventPartition,
  typename Event = typename EventPartition::parts::value_type,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const EventPartition repartition(const EventPartition& previous,
                                 const real interval) {
  const auto& parts = previous.parts;
  Event reset_parts;
  reset_parts.reserve(std::accumulate(parts.cbegin(), parts.cend(), 0ULL,
    [](const auto size, const auto& event) { return size + event.size(); }));

  for (const auto& event : parts) {
    reset_parts.insert(reset_parts.cend(), event.cbegin(), event.cend());
  }

  return partition(reset_parts, previous.coordinate, interval);
}
const event_partition repartition(const event_partition& previous,
                                  const real interval) {
  return repartition<event_partition, event>(previous, interval);
}
const full_event_partition repartition(const full_event_partition& previous,
                                       const real interval) {
  return repartition<full_event_partition, full_event>(previous, interval);
}
//----------------------------------------------------------------------------------------------

//__Fast Check if Points Form a Line____________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
bool fast_line_check_2d(const Event& points,
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
bool fast_line_check_3d(const Event& points,
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
bool fast_line_check(const event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2) {
  return fast_line_check_2d<>(points, threshold, x1, x2);
}
bool fast_line_check(const full_event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2) {
  return fast_line_check_2d<>(points, threshold, x1, x2);
}
bool fast_line_check(const event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2,
                     const Coordinate x3) {
  return fast_line_check_3d<>(points, threshold, x1, x2, x3);
}
bool fast_line_check(const full_event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2,
                     const Coordinate x3) {
  return fast_line_check_3d<>(points, threshold, x1, x2, x3);
}
//----------------------------------------------------------------------------------------------

//__Seeding Algorithm___________________________________________________________________________
template<class EventPartition,
  typename EventVector = typename EventPartition::parts,
  typename Event = typename EventVector::value_type,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const EventVector seed(const size_t n,
                       const EventPartition& partition,
                       const real line_threshold) {
  if (n <= 2)
    return {};

  const auto& layers = partition.parts;
  const auto layer_count = layers.size();

  // FIXME: find a close upper bound for out.reserve
  EventVector out{};

  // FIXME: unsure what to do here
  if (layer_count < n)
    return {};

  util::bit_vector_sequence layer_sequence;
  for (const auto& layer : layers)
    layer_sequence.emplace_back(1, layer.size());

  util::order2_permutations(n, layer_sequence, [&](const auto& chooser) {
    Event tuple;
    tuple.reserve(n);

    for (size_t index = 0; index < layer_count; ++index) {
      if (chooser[index]) {
        const auto& layer = layers[index];
        const auto& bits = layer_sequence[index];
        size_t j = 0;
        std::copy_if(layer.cbegin(), layer.cend(), std::back_inserter(tuple),
          [&](auto) { return bits[j++]; });
      }
    }

    if (fast_line_check(t_sort(tuple), line_threshold, Coordinate::X, Coordinate::Y, Coordinate::Z))
      out.push_back(tuple);
  });

  return out;
}
const event_vector seed(const size_t n,
                        const event_partition& partition,
                        const real line_threshold) {
  return seed<event_partition, event_vector>(n, partition, line_threshold);
}
const full_event_vector seed(const size_t n,
                             const full_event_partition& partition,
                             const real line_threshold) {
  return seed<full_event_partition, full_event_vector>(n, partition, line_threshold);
}
//----------------------------------------------------------------------------------------------

//__Check if Seeds can be Joined________________________________________________________________
bool seeds_compatible(const event& first,
                      const event& second,
                      const size_t difference) {
  return std::equal(first.cbegin() + difference, first.cend(), second.cbegin());
}
//----------------------------------------------------------------------------------------------

//__Join Two Seeds______________________________________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const Event join(const Event& first,
                 const Event& second,
                 const size_t difference) {
  const auto second_size = second.size();
  const auto overlap = first.size() - difference;

  if (overlap <= 0 || second_size < overlap)
    return {};

  const auto size = difference + second_size;
  Event out;
  out.reserve(size);

  size_t index = 0;
  for (; index < difference; ++index) out.push_back(first[index]);
  for (; index < difference + overlap; ++index) {
    const auto& point = first[index];
    if (point != second[index - difference])
      return {};
    out.push_back(point);
  }
  for (; index < size; ++index) out.push_back(second[index - difference]);

  return out;
}
const event join(const event& first,
                 const event& second,
                 const size_t difference) {
  return join<>(first, second, difference);
}
const full_event join(const full_event& first,
                      const full_event& second,
                      const size_t difference) {
  return join<>(first, second, difference);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Index Representation of Event Vector________________________________________________________
using index_vector = std::vector<std::size_t>;
//----------------------------------------------------------------------------------------------

//__Join All Secondaries matching Seed__________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
void _join_secondaries(const size_t seed_index,
                       const size_t difference,
                       EventVector& seed_buffer,
                       const index_vector& indicies,
                       util::bit_vector& join_list,
                       index_vector& out) {
  const auto& seed = seed_buffer[indicies[seed_index]];
  const auto size = indicies.size();
  for (size_t i = 0; i < size; ++i) {
    const auto& next_seed = seed_buffer[indicies[i]];
    const auto& joined_seed = join(seed, next_seed, difference);
    if (!joined_seed.empty()) {
      seed_buffer.push_back(joined_seed);
      out.push_back(seed_buffer.size() - 1);
      join_list[i] = true;
      join_list[seed_index] = true;
    }
  }
}
//----------------------------------------------------------------------------------------------

//__Queue for Joinable Seeds____________________________________________________________________
using seed_queue = std::queue<index_vector>;
//----------------------------------------------------------------------------------------------

//__Partial Join Seeds from Seed Buffer_________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
bool _partial_join(EventVector& seed_buffer,
                   const index_vector& indicies,
                   const size_t difference,
                   seed_queue& joined,
                   seed_queue& singular,
                   EventVector& out) {
  const auto size = indicies.size();
  if (size <= 1) return false;

  util::bit_vector join_list(size);
  index_vector to_joined, to_singular;
  to_joined.reserve(size);
  to_singular.reserve(size);

  for (size_t index = 0; index < size; ++index)
    _join_secondaries(index, difference, seed_buffer, indicies, join_list, to_joined);

  if (!to_joined.empty()) {
    for (size_t i = 0; i < size; ++i)
      if (!join_list[i])
        to_singular.push_back(indicies[i]);

    joined.push(to_joined);
    singular.push(to_singular);
    return true;
  }

  for (size_t i = 0; i < size; ++i)
    out.push_back(seed_buffer[indicies[i]]);

  return false;
}
//----------------------------------------------------------------------------------------------

//__Join Seeds from Seed Queues_________________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
void _join_next_in_queue(seed_queue& queue,
                         EventVector& seed_buffer,
                         const size_t difference,
                         seed_queue& joined,
                         seed_queue& singular,
                         EventVector& out) {
  if (!queue.empty()) {
    const auto indicies = std::move(queue.front());
    queue.pop();
    _partial_join(seed_buffer, indicies, difference, joined, singular, out);
  }
}
//----------------------------------------------------------------------------------------------

//__Seed Join Loop______________________________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
void _full_join(EventVector& seed_buffer,
                const size_t difference,
                const size_t max_difference,
                seed_queue& joined,
                seed_queue& singular,
                EventVector& out) {
  do {
    _join_next_in_queue(joined, seed_buffer, difference, joined, singular, out);
    _join_next_in_queue(singular, seed_buffer, difference + 1, joined, singular, out);
  } while (!joined.empty() || !singular.empty());
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Seed Join___________________________________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
const EventVector join_all(const EventVector& seeds) {
  const auto size = seeds.size();

  EventVector out;
  out.reserve(size); // FIXME: bad estimate?

  EventVector seed_buffer;
  std::copy(seeds.cbegin(), seeds.cend(), std::back_inserter(seed_buffer));

  seed_queue joined, singular;

  index_vector initial_seeds;
  initial_seeds.reserve(size);
  for (size_t i = 0; i < size; ++i)
    initial_seeds.push_back(i);

  joined.push(initial_seeds);
  _full_join(seed_buffer, 1, 2, joined, singular, out);
  return out;
}
const event_vector join_all(const event_vector& seeds) {
  return join_all<>(seeds);
}
const full_event_vector join_all(const full_event_vector& seeds) {
  return join_all<>(seeds);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Calculate Squared Residual of Track wrt Full Hit____________________________________________
real _track_squared_residual(const real t0,
                             const real x0,
                             const real y0,
                             const real z0,
                             const real vx,
                             const real vy,
                             const real vz,
                             const full_hit& point) {
  const auto dt = (point.z - z0) / vz;
  const auto t_res = (dt + t0 - point.t) / point.error.t;
  const auto x_res = (std::fma(dt, vx, x0) - point.x) / point.error.x;
  const auto y_res = (std::fma(dt, vy, y0) - point.y) / point.error.y;
  return t_res*t_res + 12*x_res*x_res + 12*y_res*y_res;
}
//----------------------------------------------------------------------------------------------

//__Track Parameter Type________________________________________________________________________
struct _track_parameters { fit_parameter t0, x0, y0, z0, vx, vy, vz; };
//----------------------------------------------------------------------------------------------

//__Fast Guess of Initial Track Parameters______________________________________________________
_track_parameters _guess_track(const full_event& points) {
  const auto& first = points.front();
  const auto& last = points.back();
  const auto dt = last.t - first.t;
  const auto time_error = first.error.t;
  return {{first.t,                 time_error,         0, 0},
          {first.x,                 100*units::length,  0, 0},
          {first.y,                 100*units::length,  0, 0},
          {first.z,                 100*units::length,  0, 0},
          {(last.x - first.x) / dt, 50*units::velocity, 0, 0},
          {(last.y - first.y) / dt, 50*units::velocity, 0, 0},
          {(last.z - first.z) / dt, 50*units::velocity, 0, 0}};
}
//----------------------------------------------------------------------------------------------

//__Gaussian Negative Log Likelihood Calculation________________________________________________
thread_local full_event&& _nll_fit_event = {};
void _gaussian_nll(Int_t&, Double_t*, Double_t& out, Double_t* x, Int_t) {
  out = 0.5L * std::accumulate(_nll_fit_event.cbegin(), _nll_fit_event.cend(), 0.0L,
    [&](const auto sum, const auto& point) {
      return sum + _track_squared_residual(x[0], x[1], x[2], x[3], x[4], x[5], x[6], point); });
}
//----------------------------------------------------------------------------------------------

//__MINUIT Gaussian Fitter______________________________________________________________________
_track_parameters& _fit_event_minuit(const full_event& points,
                                     _track_parameters& parameters,
                                     const fit_settings& settings,
                                     const Coordinate fixed=Coordinate::Z) {
  TMinuit minuit;
  minuit.SetGraphicsMode(settings.graphics_on);
  minuit.SetPrintLevel(settings.print_level);
  minuit.SetErrorDef(settings.error_def);
  minuit.SetMaxIterations(settings.max_iterations);

  minuit.Command("SET STR 2");

  auto& t0 = parameters.t0;
  minuit.DefineParameter(0, "T0", t0.value, t0.error, t0.min, t0.max);
  auto& x0 = parameters.x0;
  minuit.DefineParameter(1, "X0", x0.value, x0.error, x0.min, x0.max);
  auto& y0 = parameters.y0;
  minuit.DefineParameter(2, "Y0", y0.value, y0.error, y0.min, y0.max);
  auto& z0 = parameters.z0;
  minuit.DefineParameter(3, "Z0", z0.value, z0.error, z0.min, z0.max);
  auto& vx = parameters.vx;
  minuit.DefineParameter(4, "VX", vx.value, vx.error, vx.min, vx.max);
  auto& vy = parameters.vy;
  minuit.DefineParameter(5, "VY", vy.value, vy.error, vy.min, vy.max);
  auto& vz = parameters.vz;
  minuit.DefineParameter(6, "VZ", vz.value, vz.error, vz.min, vz.max);

  switch (fixed) {
    case Coordinate::T: minuit.FixParameter(0); break;
    case Coordinate::X: minuit.FixParameter(1); break;
    case Coordinate::Y: minuit.FixParameter(2); break;
    case Coordinate::Z: minuit.FixParameter(3); break;
  }

  _nll_fit_event = points;
  minuit.SetFCN(_gaussian_nll);

  Int_t error_flag;
  auto command_parameters = settings.command_parameters;
  minuit.mnexcm(
    settings.command_name.c_str(),
    command_parameters.data(),
    command_parameters.size(),
    error_flag);

  switch (error_flag) {
    case 1:
    case 2:
    case 3: util::error::exit("[FATAL ERROR] Unknown MINUIT Command \"", settings.command_name,
                              "\". Exited with Error Code ", error_flag, ".\n");
    //case 4: util::error::exit("[FATAL ERROR] MINUIT Exited Abnormally ",
    //                          "with Error Code ", error_flag, ".\n");
    default: break;
  }

  Double_t value, error;
  minuit.GetParameter(0, value, error);
  t0.value = value;
  t0.error = error;
  minuit.GetParameter(1, value, error);
  x0.value = value;
  x0.error = error;
  minuit.GetParameter(2, value, error);
  y0.value = value;
  y0.error = error;
  minuit.GetParameter(3, value, error);
  z0.value = value;
  z0.error = error;
  minuit.GetParameter(4, value, error);
  vx.value = value;
  vx.error = error;
  minuit.GetParameter(5, value, error);
  vy.value = value;
  vy.error = error;
  minuit.GetParameter(6, value, error);
  vz.value = value;
  vz.error = error;

  // save covariance matrix
  // mnemat

  return parameters;
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

track::track(const std::vector<hit>& points,
             const fit_settings& settings)
    : track(add_errors(points), settings) {}

//__Track Constructor___________________________________________________________________________
track::track(const std::vector<full_hit>& points,
             const fit_settings& settings)
    : _full_event(points), _settings(settings) {

  auto fit_track = _guess_track(_full_event);
  _fit_event_minuit(_full_event, fit_track, _settings);

  _t0 = std::move(fit_track.t0);
  _x0 = std::move(fit_track.x0);
  _y0 = std::move(fit_track.y0);
  _z0 = std::move(fit_track.z0);
  _vx = std::move(fit_track.vx);
  _vy = std::move(fit_track.vy);
  _vz = std::move(fit_track.vz);

  const auto& full_event_begin = _full_event.cbegin();
  const auto& full_event_end = _full_event.cend();

  std::transform(full_event_begin, full_event_end, std::back_inserter(_delta_chi_squared),
    [&](const auto& point) {
      return _track_squared_residual(
        _t0.value, _x0.value, _y0.value, _z0.value, _vx.value, _vy.value, _vz.value,
        point);
    });

  std::transform(full_event_begin, full_event_end, std::back_inserter(_detectors),
    [&](const auto& point) { return geometry::volume(reduce_to_r3(point)); });
}
//----------------------------------------------------------------------------------------------

//__Get Position of Track at Fixed Z____________________________________________________________
const hit track::operator()(const real z) const {
  const auto dt = (z - _z0.value) / _vz.value;
  return { dt + _t0.value, std::fma(dt, _vx.value, _x0.value), std::fma(dt, _vy.value, _y0.value), z };
}
//----------------------------------------------------------------------------------------------

//__Relativistic Beta for the Track_____________________________________________________________
real track::beta() const {
  return std::sqrt(util::math::fused_product(
    _vx.value, _vx.value, _vy.value, _vy.value, _vz.value, _vz.value)) / units::speed_of_light;
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared Test Statistic__________________________________________________________________
real track::chi_squared() const {
  return std::accumulate(_delta_chi_squared.cbegin(), _delta_chi_squared.cend(), 0.0L);
}
//----------------------------------------------------------------------------------------------

//__Track Degrees of Freedom____________________________________________________________________
size_t track::degrees_of_freedom() const {
  return 3 * _full_event.size() - 6;
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared per Degree of Freedom___________________________________________________________
real track::chi_squared_per_dof() const {
  return chi_squared() / degrees_of_freedom();
}
//----------------------------------------------------------------------------------------------

//__Get Event from Track________________________________________________________________________
const event track::event() const {
  std::vector<hit> out;
  out.reserve(_full_event.size());
  std::transform(_full_event.cbegin(), _full_event.cend(), std::back_inserter(out),
    [](const auto& full_point) { return hit{full_point.t, full_point.x, full_point.y, full_point.z}; });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Get Front of Event from Track_______________________________________________________________
const hit track::front() const {
  return (*this)(_full_event.front().z);
}
//----------------------------------------------------------------------------------------------

//__Get Back of Event from Track________________________________________________________________
const hit track::back() const {
  return (*this)(_full_event.back().z);
}
//----------------------------------------------------------------------------------------------

//__Track Output Stream Operator________________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const track& track) {
  os.precision(7);
  os << "Track Parameters: \n"
     << "  T0: " << track.t0_value() << "  (+/- " << track.t0_error() << ")\n"
     << "  X0: " << track.x0_value() << "  (+/- " << track.x0_error() << ")\n"
     << "  Y0: " << track.y0_value() << "  (+/- " << track.y0_error() << ")\n"
     << "  Z0: " << track.z0_value() << "  (+/- " << track.z0_error() << ")\n"
     << "  VX: " << track.vx_value() << "  (+/- " << track.vx_error() << ")\n"
     << "  VY: " << track.vy_value() << "  (+/- " << track.vy_error() << ")\n"
     << "  VZ: " << track.vz_value() << "  (+/- " << track.vz_error() << ")\n";

  os.precision(6);
  os << "Event: \n";
  const auto points = track.event();
  const auto detectors = track.detectors();
  const auto size = points.size();
  for (size_t i = 0; i < size; ++i) {
    os << "  " << detectors[i] << " " << points[i] << "\n";
  }

  os.precision(7);
  os << "Statistics: \n"
     << "  chi2:     " << track.chi_squared() << " = ";
  util::io::print_range(track.chi_squared_vector(), " + ", "", os) << "\n";
  os << "  dof:      " << track.degrees_of_freedom()               << "\n"
     << "  chi2/dof: " << track.chi_squared_per_dof()              << "\n";

  os.precision(6);
  os << "Dynamics: \n"
     << "  beta:  " << track.beta()  << "\n"
     << "  front: " << track.front() << "\n"
     << "  back:  " << track.back()  << "\n";

  return os;
}
//----------------------------------------------------------------------------------------------

//__Fit all Seeds to Tracks_____________________________________________________________________
template<class EventVector,
    typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
const track_vector fit_seeds(const EventVector& seeds,
                             const fit_settings& settings) {
  track_vector out;
  out.reserve(seeds.size());
  for (const auto& seed : seeds)
    out.emplace_back(seed, settings);
  return out;
}
const track_vector fit_seeds(const event_vector& seeds,
                             const fit_settings& settings) {
  return fit_seeds<>(seeds, settings);
}
const track_vector fit_seeds(const full_event_vector& seeds,
                             const fit_settings& settings) {
  return fit_seeds<>(seeds, settings);
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

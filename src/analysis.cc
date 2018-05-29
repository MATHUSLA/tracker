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

#include "analysis.hh"

#include <numeric>
#include <queue>

#include "ROOT/TMinuit.h"

#include "geometry.hh"
#include "units.hh"

#include "util/bit_vector.hh"
#include "util/io.hh"
#include "util/math.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Average Point_______________________________________________________________________________
const r4_point mean(const event_points& points) {
  const auto size = points.size();
  return (size == 0) ? r4_point{} : std::accumulate(points.cbegin(), points.cend(), r4_point{}) / size;
}
//----------------------------------------------------------------------------------------------

//__Row-Major Covariance Matrix_________________________________________________________________
const r4_point_vector covariance_matrix(const event_points& points) {
  const auto size = points.size();
  if (size == 0) return {};

  const auto inv_size = 1.0L / size;

  const auto& average = inv_size * std::accumulate(points.cbegin(), points.cend(), r4_point{});
  const auto& points_T = transpose(points);

  const real covTT = inv_size * (points_T.ts * points_T.ts) - average.t * average.t;
  const real covXX = inv_size * (points_T.xs * points_T.xs) - average.x * average.x;
  const real covYY = inv_size * (points_T.ys * points_T.ys) - average.y * average.y;
  const real covZZ = inv_size * (points_T.zs * points_T.zs) - average.z * average.z;

  const real covTX = inv_size * (points_T.ts * points_T.xs) - average.t * average.x;
  const real covTY = inv_size * (points_T.ts * points_T.ys) - average.t * average.y;
  const real covTZ = inv_size * (points_T.ts * points_T.zs) - average.t * average.z;

  const real covXY = inv_size * (points_T.xs * points_T.ys) - average.x * average.y;
  const real covXZ = inv_size * (points_T.xs * points_T.zs) - average.x * average.z;

  const real covYZ = inv_size * (points_T.ys * points_T.zs) - average.y * average.z;

  return {{covTT, covTX, covTY, covTZ},
          {covTX, covXX, covXY, covXZ},
          {covTY, covXY, covYY, covYZ},
          {covTZ, covXZ, covYZ, covZZ}};
}
//----------------------------------------------------------------------------------------------

//__Time Normalize Events_______________________________________________________________________
const event_points time_normalize(const event_points& event) {
  const auto size = event.size();
  if (size == 0) return {};
  auto out = t_copy_sort(event);
  const auto min = r4_point{event.front().t, 0, 0, 0};
  std::transform(out.cbegin(), out.cend(), out.begin(),
    [&](const auto& point) { return point - min; });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Collapse Points by R4 Interval______________________________________________________________
const event_points collapse(const event_points& event,
                            const r4_point& ds) {
  const auto size = event.size();
  if (size == 0) return {};

  event_points out;
  out.reserve(size);

  const auto& sorted_event = time_normalize(event);

  using size_type = event_points::size_type;

  size_type index = 0;
  std::queue<size_type> marked_indicies;
  while (index < size) {
    size_type collected = 1, missed_index = 0;
    const auto& point = sorted_event[index];
    const auto time_interval = point.t + ds.t;

    r4_point sum = point;

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
//----------------------------------------------------------------------------------------------

//__Partition Points by Coordinate______________________________________________________________
const event_partition partition(const event_points& points,
                                const real interval,
                                const Coordinate coordinate) {
  event_partition out{{}, coordinate};
  if (points.empty())
    return out;

  auto& parts = out.parts;
  const auto& sorted_points = coordinate_stable_copy_sort(points, coordinate);
  const auto size = sorted_points.size();

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

    parts.push_back(t_sort(current_layer));
  }

  return out;
}
//----------------------------------------------------------------------------------------------

//__Center of Geometric Object for each Point___________________________________________________
const event_points find_centers(const event_points& points) {
  event_points out;
  out.reserve(points.size());
  std::transform(points.cbegin(), points.cend(), std::back_inserter(out),
    [&](const auto& point) { return geometry::find_center(point); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Fast Check if Points Form a Line____________________________________________________________
bool fast_line_check(const event_points& points, const real threshold) {
  const auto& line_begin = points.front();
  const auto& line_end = points.back();
  return threshold >= std::accumulate(points.cbegin() + 1, points.cend() - 1, threshold,
    [&](const auto& max, const auto& point) {
        return std::max(max, point_line_distance(point, line_begin, line_end)); });
}
//----------------------------------------------------------------------------------------------

//__Seeding Algorithm___________________________________________________________________________
const event_vector seed(const size_t n,
                        const event_points& event,
                        const r4_point& collapse_ds,
                        const real layer_dz,
                        const real line_dr) {
  if (n <= 2) return {};

  const auto& points = collapse(event, collapse_ds);
  const auto size = points.size();

  if (size <= n) return { points };

  event_vector out;
  out.reserve(std::pow(size, n) / std::pow(n/2.718, n));  // FIXME: work on this limit (Stirling's approximation)

  const auto& layers = partition(points, layer_dz).parts;
  const auto layer_count = layers.size();
  if (layer_count < n) return {}; // FIXME: unsure what to do here

  util::bit_vector_sequence layer_sequence;
  for (const auto& layer : layers)
    layer_sequence.emplace_back(1, layer.size());

  util::order2_permutations(n, layer_sequence, [&](const auto& chooser) {
    event_points tuple;
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

    if (fast_line_check(t_sort(tuple), line_dr))
      out.push_back(tuple);
  });

  return out;
}
//----------------------------------------------------------------------------------------------

//__Check if Seeds can be Joined________________________________________________________________
bool seeds_compatible(const event_points& first, const event_points& second, const size_t difference) {
  return std::equal(first.cbegin() + difference, first.cend(), second.cbegin(), second.cend() - difference);
}
//----------------------------------------------------------------------------------------------

//__Join Two Seeds______________________________________________________________________________
const event_points join(const event_points& first,
                        const event_points& second,
                        const size_t difference) {
  const auto first_size = first.size();
  const auto second_size = second.size();
  const auto overlap = first_size - difference;

  if (overlap <= 0 || second_size < overlap)
    return {};

  const auto size = difference + second_size;
  event_points out;
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
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Index Representation of Event Vector________________________________________________________
using index_vector = std::vector<std::size_t>;
//----------------------------------------------------------------------------------------------

//__Join All Secondaries matching Seed__________________________________________________________
void _join_secondaries(const size_t seed_index,
                       const size_t difference,
                       event_vector& seed_buffer,
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
bool _partial_join(event_vector& seed_buffer,
                   const index_vector& indicies,
                   const size_t difference,
                   seed_queue& joined,
                   seed_queue& singular,
                   event_vector& out) {
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
void _join_next_in_queue(seed_queue& queue,
                         event_vector& seed_buffer,
                         const size_t difference,
                         seed_queue& joined,
                         seed_queue& singular,
                         event_vector& out) {
  if (!queue.empty()) {
    const auto indicies = std::move(queue.front());
    queue.pop();
    _partial_join(seed_buffer, indicies, difference, joined, singular, out);
  }
}
//----------------------------------------------------------------------------------------------

//__Seed Join Loop______________________________________________________________________________
void _full_join(event_vector& seed_buffer,
                const size_t difference,
                const size_t max_difference,
                seed_queue& joined,
                seed_queue& singular,
                event_vector& out) {
  do {
    _join_next_in_queue(joined, seed_buffer, difference, joined, singular, out);
    _join_next_in_queue(singular, seed_buffer, difference + 1, joined, singular, out);
  } while (!joined.empty() || !singular.empty());
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Seed Join___________________________________________________________________________________
const event_vector join_all(const event_vector& seeds) {
  const auto size = seeds.size();

  event_vector out;
  out.reserve(size);  // FIXME: bad estimate?

  event_vector seed_buffer;
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
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Calculate Squared Residual of Track Through a Volume________________________________________
real _track_squared_residual(const real t0,
                             const real x0,
                             const real y0,
                             const real z0,
                             const real vx,
                             const real vy,
                             const real vz,
                             const r4_point& point) {
  const auto volume = geometry::volume(point);
  const auto limits = geometry::limits_of(volume);
  const auto& center = limits.center;
  const auto& min = limits.min;
  const auto& max = limits.max;
  const auto dt = (center.z - z0) / vz;
  const auto t_res = (dt + t0 - point.t) / (2 * units::time);
  const auto x_res = (std::fma(dt, vx, x0) - center.x) / (max.x - min.x);
  const auto y_res = (std::fma(dt, vy, y0) - center.y) / (max.y - min.y);
  return t_res*t_res + 12*x_res*x_res + 12*y_res*y_res;
}
//----------------------------------------------------------------------------------------------

//__Track Parameter Type________________________________________________________________________
struct _track_parameters { fit_parameter t0, x0, y0, z0, vx, vy, vz; };
//----------------------------------------------------------------------------------------------

//__Fast Guess of Initial Track Parameters______________________________________________________
_track_parameters _guess_track(const event_points& event) {
  const auto& first = event.front();
  const auto& last = event.back();
  const auto dt = last.t - first.t;
  return {{first.t, 2*units::time, 0, 0},
          {first.x, 100*units::length, 0, 0},
          {first.y, 100*units::length, 0, 0},
          {first.z, 100*units::length, 0, 0},
          {(last.x - first.x) / dt, 0.01*units::speed_of_light, 0, 0},
          {(last.y - first.y) / dt, 0.01*units::speed_of_light, 0, 0},
          {(last.z - first.z) / dt, 0.01*units::speed_of_light, 0, 0}};
}
//----------------------------------------------------------------------------------------------

//__Gaussian Negative Log Likelihood Calculation________________________________________________
thread_local event_points&& _nll_fit_event = {};
void _gaussian_nll(Int_t&, Double_t*, Double_t& out, Double_t* parameters, Int_t) {
  out = 0.5L * std::accumulate(_nll_fit_event.cbegin(), _nll_fit_event.cend(), 0.0L,
    [&](const auto& sum, const auto& point) {
      return sum + _track_squared_residual(
        parameters[0],
        parameters[1],
        parameters[2],
        parameters[3],
        parameters[4],
        parameters[5],
        parameters[6],
        point); });
}
//----------------------------------------------------------------------------------------------

//__MINUIT Gaussian Fitter______________________________________________________________________
_track_parameters& _fit_event(const event_points& event,
                              _track_parameters& parameters,
                              const fit_settings& settings) {
  TMinuit minuit;
  minuit.SetPrintLevel(-1);
  minuit.SetErrorDef(settings.error_def);
  minuit.SetMaxIterations(settings.max_iterations);

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

  _nll_fit_event = event;
  minuit.SetFCN(_gaussian_nll);

  Int_t error_flag;
  auto command_parameters = settings.command_parameters;
  minuit.mnexcm(
    settings.command_name.c_str(),
    command_parameters.data(),
    command_parameters.size(),
    error_flag);

  // TODO: read out error_flag

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

  return parameters;
}
//----------------------------------------------------------------------------------------------

//__4x4 Cholesky Factorization (only for Positive Definite Symmetric Matricies)_________________
real_array<16> _4x4_cholesky_factor(const real_array<16> m) {
  using namespace util::math;
  real_array<16> out{};
  out[0]  = std::sqrt(m[0]);
  out[4]  = m[4] / out[0];
  out[5]  = std::sqrt(m[5] - out[4] * out[4]);
  out[8]  = m[8] / out[0];
  out[9]  = (m[9] - out[8] * out[4]) / out[5];
  out[10] = std::sqrt(m[10] - out[9] * out[9] - out[8] * out[8]);
  out[12] = m[12] / out[0];
  out[13] = (m[13] - out[12] * out[4] ) / out[5];
  out[14] = (m[14] - out[13] * out[9] - out[12] * out[8]) / out[10];
  out[15] = std::sqrt(m[15] - out[14] * out[14] - out[13] * out[13] - out[12] * out[12]);
  // std::cout.precision(15);
  // std::cout << m[10] << " " << out[9] << " " << out[8] << ": " << m[10] - fused_product(out[9], out[9], out[8], out[8]) << "\n";
  // util::io::print_range(m) << "\n";
  return out;
}
//----------------------------------------------------------------------------------------------

//__4x4 Inverse Using Cholesky Factorization (only for Positive Definite Symmetric Matricies)___
real_array<16> _4x4_cholesky_inverse(const real_array<16> m) {
  const auto L = _4x4_cholesky_factor(m);

  real_array<16> L_inv{};
  L_inv[0]  = 1.0L / L[0];
  L_inv[5]  = 1.0L / L[5];
  L_inv[10] = 1.0L / L[10];
  L_inv[15] = 1.0L / L[15];

  L_inv[4]  = -L_inv[5]  *  L[4]  * L_inv[0];
  L_inv[8]  = -L_inv[10] * (L[8]  * L_inv[0]
                         +  L[9]  * L_inv[4]);
  L_inv[9]  = -L_inv[10] *  L[9]  * L_inv[5];
  L_inv[12] = -L_inv[15] * (L[12] * L_inv[0]
                         +  L[13] * L_inv[4]
                         +  L[14] * L_inv[8]);
  L_inv[13] = -L_inv[15] * (L[13] * L_inv[5]
                         +  L[14] * L_inv[9]);
  L_inv[14] = -L_inv[15] *  L[14] * L_inv[10];

  real_array<16> out;
  out[0]  = L_inv[0]  * L_inv[0]
          + L_inv[4]  * L_inv[4]
          + L_inv[8]  * L_inv[8]
          + L_inv[12] * L_inv[12];
  out[5]  = L_inv[5]  * L_inv[5]
          + L_inv[9]  * L_inv[9]
          + L_inv[13] * L_inv[13];
  out[10] = L_inv[10] * L_inv[10]
          + L_inv[14] * L_inv[14];
  out[15] = L_inv[15] * L_inv[15];

  out[4]  = L_inv[4]  * L_inv[5]
          + L_inv[8]  * L_inv[9]
          + L_inv[12] * L_inv[13];
  out[8]  = L_inv[8]  * L_inv[10]
          + L_inv[12] * L_inv[14];
  out[9]  = L_inv[9]  * L_inv[10]
          + L_inv[13] * L_inv[14];
  out[12] = L_inv[12] * L_inv[15];
  out[13] = L_inv[13] * L_inv[15];
  out[14] = L_inv[14] * L_inv[15];

  out[1]  = out[4];
  out[2]  = out[8];
  out[3]  = out[12];
  out[6]  = out[9];
  out[7]  = out[13];
  out[11] = out[14];

  return out;
}
//----------------------------------------------------------------------------------------------

//__4x4 Naive Matrix Inverse____________________________________________________________________
real_array<16> _4x4_naive_inverse(const real_array<16> m) {
  using namespace util::math;

  real det;
  real_array<16> out;

  const auto m_04_09 = m[4]  * m[9];
  const auto m_04_13 = m[4]  * m[13];
  const auto m_05_08 = m[5]  * m[8];
  const auto m_05_12 = m[5]  * m[12];
  const auto m_06_11 = m[6]  * m[11];
  const auto m_06_15 = m[6]  * m[15];
  const auto m_07_10 = m[7]  * m[10];
  const auto m_07_14 = m[7]  * m[14];
  const auto m_08_13 = m[8]  * m[13];
  const auto m_09_12 = m[9]  * m[12];
  const auto m_10_15 = m[10] * m[15];
  const auto m_11_14 = m[11] * m[14];

  out[0]  = fused_product(m[5], m_10_15 - m_11_14, m[9],  m_07_14 - m_06_15, m[13], m_06_11 - m_07_10);
  out[4]  = fused_product(m[4], m_11_14 - m_10_15, m[8],  m_06_15 - m_07_14, m[12], m_07_10 - m_06_11);
  out[8]  = fused_product(m[7], m_08_13 - m_09_12, m[11], m_05_12 - m_04_13, m[15], m_04_09 - m_05_08);
  out[12] = fused_product(m[6], m_09_12 - m_08_13, m[10], m_04_13 - m_05_12, m[14], m_05_08 - m_04_09);

  det = m[0] * out[0] + m[1] * out[4] + m[2] * out[8] + m[3] * out[12];
  if (det == 0) return {};

  const auto m_00_05 = m[0] * m[5];
  const auto m_00_09 = m[0] * m[9];
  const auto m_00_13 = m[0] * m[13];
  const auto m_01_04 = m[1] * m[4];
  const auto m_01_08 = m[1] * m[8];
  const auto m_01_12 = m[1] * m[12];
  const auto m_02_07 = m[2] * m[7];
  const auto m_02_11 = m[2] * m[11];
  const auto m_02_15 = m[2] * m[15];
  const auto m_03_06 = m[3] * m[6];
  const auto m_03_10 = m[3] * m[10];
  const auto m_03_14 = m[3] * m[14];
  const auto m_05_10 = m[5] * m[10];

  out[1]  = fused_product(m[1],  m_11_14 - m_10_15, m[9],  m_02_15 - m_03_14, m[13], m_03_10 - m_02_11);
  out[5]  = fused_product(m[0],  m_10_15 - m_11_14, m[8],  m_03_14 - m_02_15, m[12], m_02_11 - m_03_10);
  out[2]  = fused_product(m[1],  m_06_15 - m_07_14, m[5],  m_03_14 - m_02_15, m[13], m_02_07 - m_03_06);
  out[6]  = fused_product(m[0],  m_07_14 - m_06_15, m[4],  m_02_15 - m_03_14, m[12], m_03_06 - m_02_07);
  out[3]  = fused_product(m[1],  m_07_10 - m_06_11, m[5],  m_02_11 - m_03_10, m[9],  m_03_06 - m_02_07);
  out[7]  = fused_product(m[0],  m_06_11 - m_07_10, m[4],  m_03_10 - m_02_11, m[8],  m_02_07 - m_03_06);
  out[9]  = fused_product(m[3],  m_09_12 - m_08_13, m[11], m_00_13 - m_01_12, m[15], m_01_08 - m_00_09);
  out[13] = fused_product(m[2],  m_08_13 - m_09_12, m[10], m_01_12 - m_00_13, m[14], m_00_09 - m_01_08);
  out[10] = fused_product(m[3],  m_04_13 - m_05_12, m[7],  m_01_12 - m_00_13, m[15], m_00_05 - m_01_04);
  out[14] = fused_product(m[2],  m_05_12 - m_04_13, m[6],  m_00_13 - m_01_12, m[14], m_01_04 - m_00_05);
  out[11] = fused_product(m[3],  m_05_08 - m_04_09, m[7],  m_00_09 - m_01_08, m[11], m_01_04 - m_00_05);
  out[15] = fused_product(m[2],  m_04_09 - m_05_08, m[6],  m_01_08 - m_00_09, m[10], m_00_05 - m_01_04);

  const auto inv_det = 1.0L / det;;
  std::transform(out.cbegin(), out.cend(), out.begin(), [&](const auto& value) { return value * inv_det; });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Compute Weighted Inner Product______________________________________________________________
real _weighted_inner_product(const real_array<16> weight,
                             const real_array<4> left,
                             const real_array<4> right) {
  real out;
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      out += weight[4*i+j] * left[i] * right[j];
    }
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Compute Weighted Length_____________________________________________________________________
real _weighted_length(const real_array<16> weight,
                      const real_array<4> vector) {
  return _weighted_inner_product(weight, vector, vector);
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Track Constructor___________________________________________________________________________
track::track(const event_points& event, const fit_settings& settings)
    : _event(event), _settings(settings) {

  auto fit_track = _guess_track(_event);
  _fit_event(_event, fit_track, _settings);

  _t0 = std::move(fit_track.t0);
  _x0 = std::move(fit_track.x0);
  _y0 = std::move(fit_track.y0);
  _z0 = std::move(fit_track.z0);
  _vx = std::move(fit_track.vx);
  _vy = std::move(fit_track.vy);
  _vz = std::move(fit_track.vz);

  const auto& covariance = covariance_matrix(_event);
  const auto& row0 = covariance[0];
  const auto& row1 = covariance[1];
  const auto& row2 = covariance[2];
  const auto& row3 = covariance[3];

  const auto& inverse_covariance = _4x4_cholesky_inverse({
    row0.t, row0.x, row0.y, row0.z,
    row1.t, row1.x, row1.y, row1.z,
    row2.t, row2.x, row2.y, row2.z,
    row3.t, row3.x, row3.y, row3.z});

  const auto& event_begin = _event.cbegin();
  const auto& event_end = _event.cend();

  std::transform(event_begin, event_end, std::back_inserter(_delta_chi_squared),
    [&](const auto& point) {
      const auto delta = point - (*this)(point.z);
      return _weighted_length(inverse_covariance, {delta.t, delta.x, delta.y, delta.z});
    });

  std::transform(event_begin, event_end, std::back_inserter(_detectors),
    [&](const auto& point) { return geometry::volume(point); });

  std::transform(event_begin, event_end, std::back_inserter(_squared_residuals),
    [&](const auto& point) {
      return _track_squared_residual(
        _t0.value,
        _x0.value,
        _y0.value,
        _z0.value,
        _vx.value,
        _vy.value,
        _vz.value,
        point);
    });
}
//----------------------------------------------------------------------------------------------

//__Get Position of Track at Fixed Z____________________________________________________________
const r4_point track::operator()(const real z) const {
  const auto dt = (z - _z0.value) / _vz.value;
  return { dt + _t0.value, std::fma(dt, _vx.value, _x0.value), std::fma(dt, _vy.value, _y0.value), z };
}
//----------------------------------------------------------------------------------------------

//__Total Residual______________________________________________________________________________
real track::residual() const {
  return std::sqrt(squared_residual());
}
//----------------------------------------------------------------------------------------------

//__Total Squared-Residual______________________________________________________________________
real track::squared_residual() const {
  return std::accumulate(_squared_residuals.cbegin(), _squared_residuals.cend(), 0.0L);
}
//----------------------------------------------------------------------------------------------

//__Residual at Each Volume_____________________________________________________________________
const real_vector track::residual_vector() const {
  real_vector out;
  out.reserve(_squared_residuals.size());
  std::transform(_squared_residuals.cbegin(), _squared_residuals.cend(), std::back_inserter(out),
    [&](const auto& sq_res) { return std::sqrt(sq_res); });
  return out;
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
  return 4 * _event.size() - 6;
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared per Degree of Freedom___________________________________________________________
real track::chi_squared_per_dof() const {
  return chi_squared() / degrees_of_freedom();
}
//----------------------------------------------------------------------------------------------

//__Output Stream Operator______________________________________________________________________
std::ostream& operator<<(std::ostream& os, const track& t) {
  return os << "Track Parameters: \n"
    << "  T0: " << t.t0_value() << "  (+/- " << t.t0_error() << ")\n"
    << "  X0: " << t.x0_value() << "  (+/- " << t.x0_error() << ")\n"
    << "  Y0: " << t.y0_value() << "  (+/- " << t.y0_error() << ")\n"
    << "  Z0: " << t.z0_value() << "  (+/- " << t.z0_error() << ")\n"
    << "  VX: " << t.vx_value() << "  (+/- " << t.vx_error() << ")\n"
    << "  VY: " << t.vy_value() << "  (+/- " << t.vy_error() << ")\n"
    << "  VZ: " << t.vz_value() << "  (+/- " << t.vz_error() << ")\n";
}
//----------------------------------------------------------------------------------------------

//__Add Track from Seed to Track Vector_________________________________________________________
track_vector& operator+=(track_vector& tracks, const event_points& seed) {
  if (tracks.empty()) tracks.emplace_back(seed);
  else tracks.emplace_back(seed, tracks.front().settings());
  return tracks;
}
//----------------------------------------------------------------------------------------------

//__Fit all Seeds to Tracks_____________________________________________________________________
const track_vector fit_seeds(const event_vector& seeds,
                             const fit_settings& settings) {
  track_vector out;
  out.reserve(seeds.size());
  std::transform(seeds.cbegin(), seeds.cend(), std::back_inserter(out),
    [&](const auto& seed) { return track(seed, settings); });
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

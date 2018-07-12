/*
 * src/tracker/analysis/track.cc
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

#include <tracker/analysis/track.hh>

#include <tracker/core/stat.hh>
#include <tracker/core/units.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/bit_vector.hh>
#include <tracker/util/error.hh>
#include <tracker/util/index_vector.hh>
#include <tracker/util/io.hh>
#include <tracker/util/math.hh>

#include "../helper/analysis.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

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
  const auto t_res = (dt + t0              - point.t) / point.width.t;
  const auto x_res = (std::fma(dt, vx, x0) - point.x) / point.width.x;
  const auto y_res = (std::fma(dt, vy, y0) - point.y) / point.width.y;
  return t_res*t_res + 12.0L*x_res*x_res + 12.0L*y_res*y_res;
}
//----------------------------------------------------------------------------------------------

//__Fast Guess of Initial Track Parameters______________________________________________________
track::fit_parameters _guess_track(const full_event& points) {
  using namespace stat::type;
  using namespace stat::error;

  const auto& first = points.front();
  const auto& last = points.back();

  const uncertain_real first_t(first.t, first.width.t);
  const auto first_x = uncertain_real::from_uniform(first.x, first.width.x);
  const auto first_y = uncertain_real::from_uniform(first.y, first.width.y);
  const auto first_z = uncertain_real::from_uniform(first.z, first.width.z);

  const auto dt = uncertain_real(last.t, last.width.t)               - first_t;
  const auto dx = uncertain_real::from_uniform(last.x, last.width.x) - first_x;
  const auto dy = uncertain_real::from_uniform(last.y, last.width.y) - first_y;
  const auto dz = uncertain_real::from_uniform(last.z, last.width.z) - first_z;

  const auto vx = dx / dt;
  const auto vy = dy / dt;
  const auto vz = dz / dt;

  // NOTE: should V be constrained?
  return {{first_t, first_t.error, 0, 0},
          {first_x, first_x.error, 0, 0},
          {first_y, first_y.error, 0, 0},
          {first_z, first_z.error, 0, 0},
          {vx, vx.error, 0, 0},
          {vy, vy.error, 0, 0},
          {vz, vz.error, 0, 0}};
}
//----------------------------------------------------------------------------------------------

//__Fix V to be C_______________________________________________________________________________
real _vz_from_c(const real vx,
                const real vy) {
  static constexpr const auto c2 = units::speed_of_light * units::speed_of_light;
  return std::sqrt(c2 - vx * vx - vy * vy);
}
//----------------------------------------------------------------------------------------------

//__Gaussian Negative Log Likelihood Calculation________________________________________________
thread_local full_event&& _nll_fit_event = {};
void _gaussian_nll(Int_t&, Double_t*, Double_t& out, Double_t* x, Int_t) {
  out = 0.5L * std::accumulate(_nll_fit_event.cbegin(), _nll_fit_event.cend(), 0.0L,
    [&](const auto sum, const auto& point) {
      return sum + _track_squared_residual(x[0], x[1], x[2], x[3], x[4], x[5], x[6], point); });
}
void _gaussian_nll_two_hit_track(Int_t&, Double_t*, Double_t& out, Double_t* x, Int_t) {
  out = 0.5L * std::accumulate(_nll_fit_event.cbegin(), _nll_fit_event.cend(), 0.0L,
    [&](const auto sum, const auto& point) {
      return sum + _track_squared_residual(x[0], x[1], x[2], x[3], x[4], x[5], _vz_from_c(x[4], x[5]), point); });
}
//----------------------------------------------------------------------------------------------

//__MINUIT Gaussian Fitter______________________________________________________________________
bool _fit_event_minuit(const full_event& points,
                       const Coordinate direction,
                       track::fit_parameters& parameters,
                       track::covariance_matrix_type& covariance_matrix) {
  using namespace helper::minuit;

  auto& t0 = parameters.t0;
  auto& x0 = parameters.x0;
  auto& y0 = parameters.y0;
  auto& z0 = parameters.z0;
  auto& vx = parameters.vx;
  auto& vy = parameters.vy;
  auto& vz = parameters.vz;
  _nll_fit_event = points;

  TMinuit minuit;
  initialize(minuit, "T0", t0, "X0", x0, "Y0", y0, "Z0", z0, "VX", vx, "VY", vy);
  switch (direction) {
    case Coordinate::T: minuit.FixParameter(0); break;
    case Coordinate::X: minuit.FixParameter(1); break;
    case Coordinate::Y: minuit.FixParameter(2); break;
    case Coordinate::Z: minuit.FixParameter(3); break;
  }

  if (points.size() == 2) {
    // FIXME: error handling for -> execute(minuit, _gaussian_nll_two_hit_track);
    if (execute(minuit, _gaussian_nll_two_hit_track) == error::diverged)
      return false;
    vz.value = _vz_from_c(vx.value, vy.value);
    vz.error = stat::error::propagate(
      real_array<2>{-vx.value / vz.value, -vy.value / vz.value},
      real_array<4>{vx.error * vx.error, vx.error * vy.error,
                    vx.error * vy.error, vy.error * vy.error});
  } else {
    set_parameters(minuit, 6UL, "VZ", vz);
    // FIXME: error handling for -> execute(minuit, _gaussian_nll);
    if (execute(minuit, _gaussian_nll) == error::diverged)
      return false;
    get_parameters(minuit, 6UL, vz);
  }

  get_parameters(minuit, t0, x0, y0, z0, vx, vy);
  get_covariance<track::free_parameter_count>(minuit, covariance_matrix);

  if (points.size() == 2) {
    for (std::size_t i{}; i < 5UL; ++i) {
      covariance_matrix[30UL + i] = 0.0L;
      covariance_matrix[6UL * i + 5UL] = 0.0L;
    }
    covariance_matrix[35UL] = vz.error * vz.error;
  }

  return true;
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Track Constructor___________________________________________________________________________
track::track(const analysis::event& points,
             const Coordinate direction)
    : track(add_width(points), direction) {}
//----------------------------------------------------------------------------------------------

//__Track Constructor___________________________________________________________________________
track::track(const analysis::full_event& points,
             const Coordinate direction) {
  reset(points, direction);
}
//----------------------------------------------------------------------------------------------

//__Get Position of Track at Fixed Z____________________________________________________________
const r4_point track::operator()(const real p) const {
  switch (_direction) {
    case Coordinate::T: return at_t(p);
    case Coordinate::X: return at_x(p);
    case Coordinate::Y: return at_y(p);
    case Coordinate::Z: return at_z(p);
  }
}
//----------------------------------------------------------------------------------------------

//__Get Position of Track at Fixed T____________________________________________________________
const r4_point track::point(const real t) const {
  return at_t(t);
}
//----------------------------------------------------------------------------------------------

//__Get Error in Position of Track at Fixed T___________________________________________________
const r4_point track::point_error(const real t) const {
  return error_at_t(t);
}
//----------------------------------------------------------------------------------------------

//__Get Position of Track at Fixed T____________________________________________________________
const r4_point track::at_t(const real t) const {
  const auto dt = t - _final.t0.value;
  return {t,
          std::fma(dt, _final.vx.value, _final.x0.value),
          std::fma(dt, _final.vy.value, _final.y0.value),
          std::fma(dt, _final.vz.value, _final.z0.value)};
}
//----------------------------------------------------------------------------------------------

//__Get Position of Track at Fixed X____________________________________________________________
const r4_point track::at_x(const real x) const {
  const auto dt = (x - _final.x0.value) / _final.vx.value;
  return {dt + _final.t0.value,
          x,
          std::fma(dt, _final.vy.value, _final.y0.value),
          std::fma(dt, _final.vz.value, _final.z0.value)};
}
//----------------------------------------------------------------------------------------------

//__Get Position of Track at Fixed Y____________________________________________________________
const r4_point track::at_y(const real y) const {
  const auto dt = (y - _final.y0.value) / _final.vy.value;
  return {dt + _final.t0.value,
          std::fma(dt, _final.vx.value, _final.x0.value),
          y,
          std::fma(dt, _final.vz.value, _final.z0.value)};
}
//----------------------------------------------------------------------------------------------

//__Get Position of Track at Fixed Z____________________________________________________________
const r4_point track::at_z(const real z) const {
  const auto dt = (z - _final.z0.value) / _final.vz.value;
  return {dt + _final.t0.value,
          std::fma(dt, _final.vx.value, _final.x0.value),
          std::fma(dt, _final.vy.value, _final.y0.value),
          z};
}
//----------------------------------------------------------------------------------------------

//__Get Position of Track at Fixed Parameter____________________________________________________
const r4_point track::at(const Coordinate c,
                         const real r) const {
  switch (c) {
    case Coordinate::T: return at_t(r);
    case Coordinate::X: return at_x(r);
    case Coordinate::Y: return at_y(r);
    case Coordinate::Z: return at_z(r);
  }
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Get 3x3 Submatrix of Covariance Matrix______________________________________________________
const real_array<9> _3x3_covariance(const track& t,
                                    const track::parameter p1,
                                    const track::parameter p2,
                                    const track::parameter p3) {
  return {t.variance(p1),       t.covariance(p1, p2), t.covariance(p1, p3),
          t.covariance(p1, p2), t.variance(p2),       t.covariance(p2, p3),
          t.covariance(p1, p3), t.covariance(p2, p3), t.variance(p3)};
}
//----------------------------------------------------------------------------------------------

//__Get 4x4 Submatrix of Covariance Matrix______________________________________________________
const real_array<16> _4x4_covariance(const track& t,
                                     const track::parameter p1,
                                     const track::parameter p2,
                                     const track::parameter p3,
                                     const track::parameter p4) {
  return {t.variance(p1),       t.covariance(p1, p2), t.covariance(p1, p3), t.covariance(p1, p4),
          t.covariance(p1, p2), t.variance(p2),       t.covariance(p2, p3), t.covariance(p2, p4),
          t.covariance(p1, p3), t.covariance(p2, p3), t.variance(p3),       t.covariance(p3, p4),
          t.covariance(p1, p4), t.covariance(p2, p4), t.covariance(p3, p4), t.variance(p4)};
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Get Error in Position of Track at Fixed T___________________________________________________
const r4_point track::error_at_t(const real t) const {
  using tp = track::parameter;
  const auto dt = t - _final.t0.value;

  const real_array<3> x_gradient{-_final.vx.value, 1.0L, dt};
  const real_array<3> y_gradient{-_final.vy.value, 1.0L, dt};
  const real_array<3> z_gradient{-_final.vz.value, 1.0L, dt};

  return {0.0L,
          stat::error::propagate(x_gradient, _3x3_covariance(*this, tp::T0, tp::X0, tp::VX)),
          stat::error::propagate(y_gradient, _3x3_covariance(*this, tp::T0, tp::Y0, tp::VY)),
          stat::error::propagate(z_gradient, _3x3_covariance(*this, tp::T0, tp::Z0, tp::VZ))};
}
//----------------------------------------------------------------------------------------------

//__Get Error in Position of Track at Fixed X___________________________________________________
const r4_point track::error_at_x(const real x) const {
  using tp = track::parameter;
  const auto vx_inv = 1.0L / _final.vx.value;
  const auto dt = (x - _final.x0.value) * vx_inv;

  const real_array<3> t_gradient{1.0L, -vx_inv, -vx_inv * dt};

  const auto vy_by_vx = _final.vy.value * vx_inv;
  const real_array<4> y_gradient{-vy_by_vx, 1.0L, -vy_by_vx * dt, dt};

  const auto vz_by_vx = _final.vz.value * vx_inv;
  const real_array<4> z_gradient{-vz_by_vx, 1.0L, -vz_by_vx * dt, dt};

  return {stat::error::propagate(t_gradient, _3x3_covariance(*this, tp::T0, tp::X0, tp::VX)),
          0.0L,
          stat::error::propagate(y_gradient, _4x4_covariance(*this, tp::X0, tp::Y0, tp::VX, tp::VY)),
          stat::error::propagate(z_gradient, _4x4_covariance(*this, tp::X0, tp::Z0, tp::VX, tp::VZ))};
}
//----------------------------------------------------------------------------------------------

//__Get Error in Position of Track at Fixed Y___________________________________________________
const r4_point track::error_at_y(const real y) const {
  using tp = track::parameter;
  const auto vy_inv = 1.0L / _final.vy.value;
  const auto dt = (y - _final.y0.value) * vy_inv;

  const real_array<3> t_gradient{1.0L, -vy_inv, -vy_inv * dt};

  const auto vx_by_vy = _final.vx.value * vy_inv;
  const real_array<4> x_gradient{-vx_by_vy, 1.0L, -vx_by_vy * dt, dt};

  const auto vz_by_vy = _final.vz.value * vy_inv;
  const real_array<4> z_gradient{-vz_by_vy, 1.0L, -vz_by_vy * dt, dt};

  return {stat::error::propagate(t_gradient, _3x3_covariance(*this, tp::T0, tp::Y0, tp::VY)),
          stat::error::propagate(x_gradient, _4x4_covariance(*this, tp::X0, tp::Y0, tp::VX, tp::VY)),
          0.0L,
          stat::error::propagate(z_gradient, _4x4_covariance(*this, tp::Y0, tp::Z0, tp::VY, tp::VZ))};
}
//----------------------------------------------------------------------------------------------

//__Get Error in Position of Track at Fixed Z___________________________________________________
const r4_point track::error_at_z(const real z) const {
  using tp = track::parameter;
  const auto vz_inv = 1.0L / _final.vz.value;
  const auto dt = (z - _final.z0.value) * vz_inv;

  const real_array<3> t_gradient{1.0L, -vz_inv, -vz_inv * dt};

  const auto vx_by_vz = _final.vx.value * vz_inv;
  const real_array<4> x_gradient{-vx_by_vz, 1.0L, -vx_by_vz * dt, dt};

  const auto vz_by_vy = _final.vz.value * vz_inv;
  const real_array<4> y_gradient{-vz_by_vy, 1.0L, -vz_by_vy * dt, dt};

  return {stat::error::propagate(t_gradient, _3x3_covariance(*this, tp::T0, tp::Z0, tp::VZ)),
          stat::error::propagate(x_gradient, _4x4_covariance(*this, tp::X0, tp::Z0, tp::VX, tp::VZ)),
          stat::error::propagate(y_gradient, _4x4_covariance(*this, tp::Y0, tp::Z0, tp::VY, tp::VZ)),
          0.0L};
}
//----------------------------------------------------------------------------------------------

//__Get Error in Position of Track at Fixed Parameter___________________________________________
const r4_point track::error_at(const Coordinate c,
                               const real r) const {
  switch (c) {
    case Coordinate::T: return error_at_t(r);
    case Coordinate::X: return error_at_x(r);
    case Coordinate::Y: return error_at_y(r);
    case Coordinate::Z: return error_at_z(r);
  }
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter from Track________________________________________________________________
const fit_parameter track::fit_of(const track::parameter p) const {
  switch (p) {
    case track::parameter::T0: return _final.t0;
    case track::parameter::X0: return _final.x0;
    case track::parameter::Y0: return _final.y0;
    case track::parameter::Z0: return _final.z0;
    case track::parameter::VX: return _final.vx;
    case track::parameter::VY: return _final.vy;
    case track::parameter::VZ: return _final.vz;
  }
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter Value from Track__________________________________________________________
real track::value(const track::parameter p) const {
  return fit_of(p).value;
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter Error from Track__________________________________________________________
real track::error(const track::parameter p) const {
  return fit_of(p).error;
}
//----------------------------------------------------------------------------------------------

//__Check if Track Fit Diverged_________________________________________________________________
bool track::fit_diverged() const noexcept {
  return _guess != _final && _final == track::fit_parameters{};
}
//----------------------------------------------------------------------------------------------

//__Relativistic Beta for the Track_____________________________________________________________
real track::beta() const {
  return util::math::hypot(_final.vx.value, _final.vy.value, _final.vz.value) / units::speed_of_light;
}
//----------------------------------------------------------------------------------------------

//__Error in Relativistic Beta for the Track____________________________________________________
real track::beta_error() const {
  static constexpr const auto c_inv = 1.0L / units::speed_of_light;
  static constexpr const auto c_inv2 = c_inv * c_inv;
  static constexpr const auto twice_c_inv = 2.0L * c_inv;
  const auto covariance = c_inv2 * _3x3_covariance(*this, parameter::VX, parameter::VY, parameter::VZ);
  const real_array<3> gradient{
    twice_c_inv * _final.vx.value, twice_c_inv * _final.vy.value, twice_c_inv * _final.vz.value};
  return stat::error::propagate(gradient, covariance);
}
//----------------------------------------------------------------------------------------------

//__Unit Vector along Track_____________________________________________________________________
const r3_point track::unit() const {
  return r3_point{_final.vx.value, _final.vy.value, _final.vz.value}
    / util::math::hypot(_final.vx.value, _final.vy.value, _final.vz.value);
}
//----------------------------------------------------------------------------------------------

//__Error in Unit Vector along Track____________________________________________________________
const r3_point track::unit_error() const {
  const auto base = 1.0L / std::pow(util::math::sum_squares(_final.vx.value, _final.vy.value, _final.vz.value), 1.5L);
  const auto covariance = _3x3_covariance(*this, parameter::VX, parameter::VY, parameter::VZ);
  const auto vxvy = _final.vx.value * _final.vy.value;
  const auto vyvz = _final.vy.value * _final.vz.value;
  const auto vzvx = _final.vz.value * _final.vx.value;
  const real_array<3> x_gradient{
    base * util::math::sum_squares(_final.vy.value, _final.vz.value),
    base * -vxvy,
    base * -vzvx};
  const real_array<3> y_gradient{
    base * -vxvy,
    base * util::math::sum_squares(_final.vz.value, _final.vx.value),
    base * -vyvz};
  const real_array<3> z_gradient{
    base * -vzvx,
    base * -vyvz,
    base * util::math::sum_squares(_final.vx.value, _final.vy.value)};
  return {stat::error::propagate(x_gradient, covariance),
          stat::error::propagate(y_gradient, covariance),
          stat::error::propagate(z_gradient, covariance)};
}
//----------------------------------------------------------------------------------------------

//__Angle of Track with Respect to Parameter____________________________________________________
real track::angle() const {
  switch (_direction) {
    case Coordinate::T: return 0.0L; // FIXME: how to handle T direction parameterization
    case Coordinate::X: return std::acos(unit().x);
    case Coordinate::Y: return std::acos(unit().y);
    case Coordinate::Z: return std::acos(unit().z);
  }
}
//----------------------------------------------------------------------------------------------

//__Error in Angle of Track with Respect to Parameter___________________________________________
real track::angle_error() const {
  const auto unit_vector = unit();
  switch (_direction) {
    case Coordinate::T: return 0.0L; // FIXME: how to handle T direction parameterization
    case Coordinate::X: return unit_error().x / std::sqrt(1.0L - unit_vector.x * unit_vector.x);
    case Coordinate::Y: return unit_error().y / std::sqrt(1.0L - unit_vector.y * unit_vector.y);
    case Coordinate::Z: return unit_error().z / std::sqrt(1.0L - unit_vector.z * unit_vector.z);
  }
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared Test Statistic__________________________________________________________________
real track::chi_squared() const {
  return std::accumulate(_delta_chi2.cbegin(), _delta_chi2.cend(), 0.0L);
}
//----------------------------------------------------------------------------------------------

//__Track Degrees of Freedom____________________________________________________________________
std::size_t track::degrees_of_freedom() const {
  const auto s = size();
  return s == 2UL ? 1UL : 3UL * s - 6UL;
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared per Degree of Freedom___________________________________________________________
real track::chi_squared_per_dof() const {
  return chi_squared() / degrees_of_freedom();
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared P-Value_________________________________________________________________________
real track::chi_squared_p_value() const {
  return stat::chi_squared_p_value(chi_squared(), degrees_of_freedom());
}
//----------------------------------------------------------------------------------------------

//__Get Variance of a Track Parameter___________________________________________________________
real track::variance(const track::parameter p) const {
  return covariance(p, p);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Get Shift Index of Track Parameters for Covariance Matrix___________________________________
constexpr std::size_t _shift_covariance_index(const track::parameter p) {
  switch (p) {
    case track::parameter::T0: return 0;
    case track::parameter::X0: return 1;
    case track::parameter::Y0: return 2;
    case track::parameter::Z0: return 2;
    case track::parameter::VX: return 3;
    case track::parameter::VY: return 4;
    case track::parameter::VZ: return 5;
  }
}
//----------------------------------------------------------------------------------------------
} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Get Covariance between Track Parameters_____________________________________________________
real track::covariance(const track::parameter p,
                       const track::parameter q) const {
  switch (_direction) {
    case Coordinate::T: {
      if (p == track::parameter::T0 || q == track::parameter::T0) {
        return 0;
      } else {
        const auto p_index = _shift_covariance_index(p)
                           - (p == track::parameter::X0 || p == track::parameter::Y0);
        const auto q_index = _shift_covariance_index(q)
                           - (q == track::parameter::X0 || q == track::parameter::Y0);
        return _covariance[6 * p_index + q_index];
      }
    } case Coordinate::X: {
      if (p == track::parameter::X0 || q == track::parameter::X0) {
        return 0;
      } else {
        const auto p_index = _shift_covariance_index(p) - (p == track::parameter::Y0);
        const auto q_index = _shift_covariance_index(q) - (q == track::parameter::Y0);
        return _covariance[6 * p_index + q_index];
      }
    } case Coordinate::Y: {
      if (p == track::parameter::Y0 || q == track::parameter::Y0) {
        return 0;
      } else {
        return _covariance[6 * _shift_covariance_index(p) + _shift_covariance_index(q)];
      }
   } case Coordinate::Z: {
      if (p == track::parameter::Z0 || q == track::parameter::Z0) {
        return 0;
      } else {
        return _covariance[6 * _shift_covariance_index(p) + _shift_covariance_index(q)];
      }
    }
  }
}
//----------------------------------------------------------------------------------------------

//__Get Front of Event from Track_______________________________________________________________
const hit track::front() const {
  switch (_direction) {
    case Coordinate::T: return at_t(full_front().t);
    case Coordinate::X: return at_x(full_front().x);
    case Coordinate::Y: return at_y(full_front().y);
    case Coordinate::Z: return at_z(full_front().z);
  }
}
//----------------------------------------------------------------------------------------------

//__Get Back of Event from Track________________________________________________________________
const hit track::back() const {
  switch (_direction) {
    case Coordinate::T: return at_t(full_back().t);
    case Coordinate::X: return at_x(full_back().x);
    case Coordinate::Y: return at_y(full_back().y);
    case Coordinate::Z: return at_z(full_back().z);
  }
}
//----------------------------------------------------------------------------------------------

//__Get Event from Track________________________________________________________________________
const std::vector<hit> track::event() const {
  std::vector<hit> out;
  out.reserve(size());
  util::algorithm::back_insert_transform(_full_event, out,
    [](const auto& point) { return reduce_to_r4(point); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Reset Track_________________________________________________________________________________
std::size_t track::reset(const analysis::event& points) {
  return reset(add_width(points));
}
//----------------------------------------------------------------------------------------------

//__Reset Track_________________________________________________________________________________
std::size_t track::reset(const analysis::full_event& points) {
  _full_event = points;
  _guess = _guess_track(_full_event);
  _final = _guess;
  if (_fit_event_minuit(_full_event, _direction, _final, _covariance)) {
    const auto& full_event_begin = _full_event.cbegin();
    const auto& full_event_end = _full_event.cend();

    std::transform(full_event_begin, full_event_end, std::back_inserter(_delta_chi2),
      [&](const auto& point) {
        return _track_squared_residual(
          t0_value(), x0_value(), y0_value(), z0_value(), vx_value(), vy_value(), vz_value(), point);
      });

    std::transform(full_event_begin, full_event_end, std::back_inserter(_detectors),
      [&](const auto& point) { return geometry::volume(reduce_to_r3(point)); });
  } else {
    _final = {};
  }
  return size();
}
//----------------------------------------------------------------------------------------------

//__Reset Track_________________________________________________________________________________
std::size_t track::reset(const analysis::event& points,
                         const Coordinate direction) {
  return reset(add_width(points), direction);
}
//----------------------------------------------------------------------------------------------

//__Reset Track_________________________________________________________________________________
std::size_t track::reset(const analysis::full_event& points,
                         const Coordinate direction) {
  _direction = direction;
  return reset(points);
}
//----------------------------------------------------------------------------------------------

//__Insert Hit into Track and Refit_____________________________________________________________
std::size_t track::insert(const hit& point) {
  return insert(add_width(point));
}
//----------------------------------------------------------------------------------------------

//__Insert Hits into Track and Refit____________________________________________________________
std::size_t track::insert(const analysis::event& points) {
  return insert(add_width(points));
}
//----------------------------------------------------------------------------------------------

//__Insert Full Hit into Track and Refit________________________________________________________
std::size_t track::insert(const full_hit& point) {
  const auto search = util::algorithm::range_binary_find_first(_full_event, point, t_ordered<analysis::full_hit>{});
  if (search != _full_event.cend()) {
    _full_event.insert(search, point);
    _full_event.shrink_to_fit();
    return reset(_full_event);
  }
  return size();
}
//----------------------------------------------------------------------------------------------

//__Insert Full Hits into Track and Refit_______________________________________________________
std::size_t track::insert(const analysis::full_event& points) {
  const auto s = size();
  analysis::full_event saved_hits;
  saved_hits.reserve(s + points.size());
  const auto sorted = util::algorithm::copy_sort_range(points, t_ordered<analysis::full_hit>{});
  std::set_union(_full_event.cbegin(), _full_event.cend(),
                 sorted.cbegin(), sorted.cend(),
                 std::back_inserter(saved_hits),
                 t_ordered<analysis::full_hit>{});
  if (saved_hits.size() != s) {
    saved_hits.shrink_to_fit();
    return reset(saved_hits);
  }
  return s;
}
//----------------------------------------------------------------------------------------------

//__Remove Hit from Track and Refit_____________________________________________________________
std::size_t track::remove(const std::size_t index) {
  const auto s = size();
  if (index >= s)
    return s;

  analysis::full_event saved_hits;
  saved_hits.reserve(s - 1);
  const auto begin = _full_event.cbegin();
  const auto end = _full_event.cend();
  saved_hits.insert(saved_hits.cend(), begin, begin + index);
  saved_hits.insert(saved_hits.cend(), begin + index + 1, end);
  return reset(saved_hits);
}
//----------------------------------------------------------------------------------------------

//__Remove Hits from Track and Refit____________________________________________________________
std::size_t track::remove(const std::vector<std::size_t>& indices) {
  const auto sorted = util::algorithm::copy_sort_range(indices);
  const auto s = size();
  analysis::full_event saved_hits;
  saved_hits.reserve(s);
  for (std::size_t hit_index{}, removal_index{}; hit_index < s; ++hit_index) {
    if (sorted[removal_index] == hit_index) {
      ++removal_index;
    } else {
      saved_hits.push_back(_full_event[hit_index]);
    }
  }
  saved_hits.shrink_to_fit();
  return reset(saved_hits);
}
//----------------------------------------------------------------------------------------------

//__Remove Hits from Track if Exceed Maximum Chi-Squared and Refit______________________________
std::size_t track::prune_on_chi_squared(const real max_chi_squared) {
  const auto s = size();
  std::vector<std::size_t> indices;
  indices.reserve(s);
  for (std::size_t i{}; i < s; ++i) {
    if (chi_squared_vector()[i] > max_chi_squared)
      indices.push_back(i);
  }
  return remove(indices);
}
//----------------------------------------------------------------------------------------------

//__Reparameterize Track in New Direction_______________________________________________________
void track::reparameterize(const Coordinate direction) {
  if (direction != _direction)
    reset(_full_event, direction);
}
//----------------------------------------------------------------------------------------------

//__Track Data Tree Constructor_________________________________________________________________
track::tree::tree(const std::string& name)
    : tree(name, name) {}
//----------------------------------------------------------------------------------------------

//__Track Data Tree Constructor_________________________________________________________________
track::tree::tree(const std::string& name,
                  const std::string& title)
    : analysis::tree(name, title),
      t0(emplace_branch<branch_value_type>("t0")),
      x0(emplace_branch<branch_value_type>("x0")),
      y0(emplace_branch<branch_value_type>("y0")),
      z0(emplace_branch<branch_value_type>("z0")),
      vx(emplace_branch<branch_value_type>("vx")),
      vy(emplace_branch<branch_value_type>("vy")),
      vz(emplace_branch<branch_value_type>("vz")),
      t0_error(emplace_branch<branch_value_type>("t0_error")),
      x0_error(emplace_branch<branch_value_type>("x0_error")),
      y0_error(emplace_branch<branch_value_type>("y0_error")),
      z0_error(emplace_branch<branch_value_type>("z0_error")),
      vx_error(emplace_branch<branch_value_type>("vx_error")),
      vy_error(emplace_branch<branch_value_type>("vy_error")),
      vz_error(emplace_branch<branch_value_type>("vz_error")),
      chi_squared(emplace_branch<branch_value_type>("chi_squared")),
      chi_squared_per_dof(emplace_branch<branch_value_type>("chi_squared_per_dof")),
      chi_squared_p_value(emplace_branch<branch_value_type>("chi_squared_p_value")),
      size(emplace_branch<branch_value_type>("size")),
      beta(emplace_branch<branch_value_type>("beta")),
      beta_error(emplace_branch<branch_value_type>("beta_error")),
      angle(emplace_branch<branch_value_type>("angle")),
      angle_error(emplace_branch<branch_value_type>("angle_error")),
      _count(emplace_branch<decltype(_count)::value_type>("N")),
      _vector_branches({t0, x0, y0, z0, vx, vy, vz,
                        t0_error, x0_error, y0_error, z0_error, vx_error, vy_error, vz_error,
                        chi_squared, chi_squared_per_dof, chi_squared_p_value,
                        size, beta, beta_error, angle, angle_error}) {}
//----------------------------------------------------------------------------------------------

//__Track Data Tree Insertion___________________________________________________________________
void track::tree::insert(const track& track) {
  t0.get().push_back(track.t0_value() / units::time);
  x0.get().push_back(track.x0_value() / units::length);
  y0.get().push_back(track.y0_value() / units::length);
  z0.get().push_back(track.z0_value() / units::length);
  vx.get().push_back(track.vx_value() / units::velocity);
  vy.get().push_back(track.vy_value() / units::velocity);
  vz.get().push_back(track.vz_value() / units::velocity);
  t0_error.get().push_back(track.t0_error() / units::time);
  x0_error.get().push_back(track.x0_error() / units::length);
  y0_error.get().push_back(track.y0_error() / units::length);
  z0_error.get().push_back(track.z0_error() / units::length);
  vx_error.get().push_back(track.vx_error() / units::velocity);
  vy_error.get().push_back(track.vy_error() / units::velocity);
  vz_error.get().push_back(track.vz_error() / units::velocity);
  chi_squared.get().push_back(track.chi_squared());
  chi_squared_per_dof.get().push_back(track.chi_squared_per_dof());
  chi_squared_p_value.get().push_back(track.chi_squared_p_value());
  size.get().push_back(track.size());
  beta.get().push_back(track.beta());
  beta_error.get().push_back(track.beta_error());
  angle.get().push_back(track.angle());
  angle_error.get().push_back(track.angle_error());
  ++_count;
}
//----------------------------------------------------------------------------------------------

//__Clear Track Data Tree_______________________________________________________________________
void track::tree::clear() {
  _count = 0UL;
  for (auto& entry : _vector_branches)
    entry.get().get().clear();
}
//----------------------------------------------------------------------------------------------

//__Reserve Space for Track Data Tree___________________________________________________________
void track::tree::reserve(std::size_t capacity) {
  for (auto& entry : _vector_branches)
    entry.get().get().reserve(capacity);
}
//----------------------------------------------------------------------------------------------

//__Fill Plots with Tracking Variables__________________________________________________________
void track::fill_plots(plot::histogram_collection& collection,
                       const track::plotting_keys& keys) const {
  if (collection.count(keys.t0)) collection[keys.t0].insert(t0_value() / units::time);
  if (collection.count(keys.x0)) collection[keys.x0].insert(x0_value() / units::length);
  if (collection.count(keys.y0)) collection[keys.y0].insert(y0_value() / units::length);
  if (collection.count(keys.z0)) collection[keys.z0].insert(z0_value() / units::length);
  if (collection.count(keys.vx)) collection[keys.vx].insert(vx_value() / units::velocity);
  if (collection.count(keys.vy)) collection[keys.vy].insert(vy_value() / units::velocity);
  if (collection.count(keys.vz)) collection[keys.vz].insert(vz_value() / units::velocity);
  if (collection.count(keys.t0_error)) collection[keys.t0_error].insert(t0_error() / units::time);
  if (collection.count(keys.x0_error)) collection[keys.x0_error].insert(x0_error() / units::length);
  if (collection.count(keys.y0_error)) collection[keys.y0_error].insert(y0_error() / units::length);
  if (collection.count(keys.z0_error)) collection[keys.z0_error].insert(z0_error() / units::length);
  if (collection.count(keys.vx_error)) collection[keys.vx_error].insert(vx_error() / units::velocity);
  if (collection.count(keys.vy_error)) collection[keys.vy_error].insert(vy_error() / units::velocity);
  if (collection.count(keys.vz_error)) collection[keys.vz_error].insert(vz_error() / units::velocity);
  if (collection.count(keys.chi_squared_per_dof)) collection[keys.chi_squared_per_dof].insert(chi_squared_per_dof());
  if (collection.count(keys.size)) collection[keys.size].insert(size());
  if (collection.count(keys.beta)) collection[keys.beta].insert(beta());
  if (collection.count(keys.beta_error)) collection[keys.beta_error].insert(beta_error());
  if (collection.count(keys.angle)) collection[keys.angle].insert(angle());
  if (collection.count(keys.angle_error)) collection[keys.angle_error].insert(angle_error());
}
//----------------------------------------------------------------------------------------------

//__Draw Fit Track______________________________________________________________________________
void track::draw(plot::canvas& canvas,
                 const real size,
                 const plot::color color,
                 const bool with_errors) const {
  if (fit_converged()) {
    canvas.add_line(front(), back(), size, color);
  }
}
//----------------------------------------------------------------------------------------------

//__Draw Guess Track____________________________________________________________________________
void track::draw_guess(plot::canvas& canvas,
                       const real size,
                       const plot::color color,
                       const bool with_errors) const {
  // TODO: finish
}
//----------------------------------------------------------------------------------------------

//__Track Output Stream Operator________________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const track& track) {
  static const std::string bar(80, '-');
  os << bar << "\n";

  if (track.fit_diverged()) {
    os << "* Track Status: " << util::io::bold << "DIVERGED" << util::io::reset_font << "\n";
    const auto guess = track.guess_fit();
    os << "* Guess Parameters: \n"
       << "    T0: " << guess.t0.value << "  (+/- " << guess.t0.error << ")\n"
       << "    X0: " << guess.x0.value << "  (+/- " << guess.x0.error << ")\n"
       << "    Y0: " << guess.y0.value << "  (+/- " << guess.y0.error << ")\n"
       << "    Z0: " << guess.z0.value << "  (+/- " << guess.z0.error << ")\n"
       << "    VX: " << guess.vx.value << "  (+/- " << guess.vx.error << ")\n"
       << "    VY: " << guess.vy.value << "  (+/- " << guess.vy.error << ")\n"
       << "    VZ: " << guess.vz.value << "  (+/- " << guess.vz.error << ")\n";
  } else {
    os << "* Track Status: " << util::io::bold << "CONVERGED" << util::io::reset_font << "\n";

    os.precision(7);
    os << "* Parameters: \n"
       << "    T0: " << track.t0_value() << "  (+/- " << track.t0_error() << ")\n"
       << "    X0: " << track.x0_value() << "  (+/- " << track.x0_error() << ")\n"
       << "    Y0: " << track.y0_value() << "  (+/- " << track.y0_error() << ")\n"
       << "    Z0: " << track.z0_value() << "  (+/- " << track.z0_error() << ")\n"
       << "    VX: " << track.vx_value() << "  (+/- " << track.vx_error() << ")\n"
       << "    VY: " << track.vy_value() << "  (+/- " << track.vy_error() << ")\n"
       << "    VZ: " << track.vz_value() << "  (+/- " << track.vz_error() << ")\n";

    os.precision(6);
    os << "* Event: \n";
    os << "    front: " << track.front() << "\n\n";
    const auto points = track.event();
    const auto& detectors = track.detectors();
    const auto size = points.size();
    for (size_t i = 0; i < size; ++i)
      os << "      " << detectors[i] << " " << points[i] << "\n";
    os << "\n    back:  " << track.back()  << "\n";

    os.precision(7);
    os << "* Statistics: \n"
       << "    dof:      " << track.degrees_of_freedom()             << "\n"
       << "    chi2:     " << track.chi_squared() << " = ";
    util::io::print_range(track.chi_squared_vector(), " + ", "", os) << "\n";
    os << "    chi2/dof: " << track.chi_squared_per_dof()            << "\n"
       << "    p-value:  " << track.chi_squared_p_value()            << "\n"
       << "    cov mat:  | ";
    const auto matrix = track.covariance_matrix();
    for (size_t i = 0; i < 6; ++i) {
      if (i > 0) os << "              | ";
      for (size_t j = 0; j < 6; ++j) {
        const auto cell = matrix[6*i+j];
        if (i == j) {
          os << util::io::bold << util::io::underline
             << cell << util::io::reset_font << " ";
        } else {
          os << cell << " ";
        }
      }
      os << "|\n";
    }

    os.precision(6);
    os << "* Dynamics: \n"
       << "    beta:  " << track.beta()  << "  (+/- " << track.beta_error()  << ")\n"
       << "    unit:  " << track.unit()  << "  (+/- " << track.unit_error()  << ")\n"
       << "    angle: " << track.angle() << "  (+/- " << track.angle_error() << ")\n";
  }

  return os << bar;
}
//----------------------------------------------------------------------------------------------

//__Fit all Seeds to Tracks_____________________________________________________________________
template<class EventVector,
  typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
const track_vector independent_fit_seeds(const EventVector& seeds,
                                         const Coordinate direction) {
  track_vector out;
  out.reserve(seeds.size());
  for (const auto& seed : seeds)
    out.emplace_back(seed, direction);
  return out;
}
const track_vector independent_fit_seeds(const event_vector& seeds,
                                         const Coordinate direction) {
  return independent_fit_seeds<>(seeds, direction);
}
const track_vector independent_fit_seeds(const full_event_vector& seeds,
                                         const Coordinate direction) {
  return independent_fit_seeds<>(seeds, direction);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Overlap Size between Events_________________________________________________________________
template<class Event,
  typename Point = typename Event::value_type,
  typename = std::enable_if_t<is_r4_type_v<Point>>>
std::size_t _overlap_size(const Event& first,
                          const Event& second) {
  const auto first_size = first.size();
  const auto second_size = second.size();

  if (first_size == 1UL && second_size == 1UL)
    return static_cast<std::size_t>(first.front() == second.front());

  const auto first_begin = std::cbegin(first);
  const auto first_end = std::cend(first);
  const auto second_begin = std::cbegin(second);
  const auto second_end = std::cend(second);

  std::size_t first_index{}, second_index{}, count{};
  while (first_index < first_size && second_index < second_size) {
    const auto& first_elem = first[first_index];
    const auto& second_elem = second[second_index];
    if (first_elem == second_elem) {
      ++count;
      ++first_index;
      ++second_index;
    } else {
      if (first_size - first_index >= second_size - second_index) {
        const auto search = util::algorithm::binary_find_first(second_begin + second_index,
                                                               second_end,
                                                               first_elem,
                                                               t_ordered<Point>{});
        if (search != second_end) {
          ++count;
          second_index = util::type::distance(second_begin, search);
        }
        ++first_index;
      } else {
        const auto search = util::algorithm::binary_find_first(first_begin + first_index,
                                                               first_end,
                                                               second_elem,
                                                               t_ordered<Point>{});
        if (search != first_end) {
          ++count;
          first_index = util::type::distance(first_begin, search);
        }
        ++second_index;
      }
    }
  }
  return count;
}
//----------------------------------------------------------------------------------------------

//__Time-Ordered Overlap Size between Events____________________________________________________
template<class Event,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
std::size_t _time_ordered_overlap_size(const Event& first,
                                       const Event& second) {
  return first.front().t <= second.front().t
    ? _overlap_size(first, second)
    : _overlap_size(second, first);
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Fit all Seeds to Tracks using Overlaps______________________________________________________
template<class EventVector,
  typename Event = typename EventVector::value_type,
  typename = std::enable_if_t<is_r4_type_v<typename Event::value_type>>>
const track_vector overlap_fit_seeds(const EventVector& seeds,
                                     const Coordinate direction,
                                     const std::size_t min_overlap) {
  const auto size = seeds.size();
  if (size == 0UL)
    return track_vector{};
  else if (size == 1UL)
    return independent_fit_seeds(seeds, direction);

  track_vector out;
  out.reserve(size);

  const auto sorted_seeds = util::algorithm::copy_sort_range(seeds, util::type::size_greater<Event>{});
  const auto track_buffer = independent_fit_seeds(sorted_seeds, direction);

  util::index_vector<> track_indices(size);
  /* FIXME: what to do here?
  util::algorithm::stable_sort_range(track_indices, [&](const auto left, const auto right) {
    return track_buffer[left].chi_squared_per_dof() > track_buffer[right].chi_squared_per_dof();
  });
  */

  util::bit_vector visited(size);

  std::size_t top_index{}, bottom_index = 1UL;
  while (top_index < size) {
    bottom_index = visited.first_unset(bottom_index);
    const auto top_track_index = track_indices[top_index];
    if (bottom_index == size) {
      out.push_back(track_buffer[top_track_index]);
      visited.set(top_index);
      top_index = visited.first_unset(1 + top_index);
      bottom_index = 1 + top_index;
      continue;
    }
    if (min_overlap <= _time_ordered_overlap_size(sorted_seeds[top_track_index],
                                                  sorted_seeds[track_indices[bottom_index]])) {
      visited.set(bottom_index);
    }
    ++bottom_index;
  }

  for (std::size_t i = 0; i < size; ++i)
    if (!visited[i])
      out.push_back(track_buffer[track_indices[i]]);

  out.shrink_to_fit();
  return out;
}
const track_vector overlap_fit_seeds(const event_vector& seeds,
                                     const Coordinate direction,
                                     const std::size_t min_overlap) {
  return overlap_fit_seeds<>(seeds, direction, min_overlap);
}
const track_vector overlap_fit_seeds(const full_event_vector& seeds,
                                     const Coordinate direction,
                                     const std::size_t min_overlap) {
  return overlap_fit_seeds<>(seeds, direction, min_overlap);
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

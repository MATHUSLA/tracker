/*
 * src/tracker/track.cc
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

#include <tracker/track.hh>

#include <ROOT/TMinuit.h>

#include <tracker/geometry.hh>
#include <tracker/stat.hh>
#include <tracker/units.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/error.hh>
#include <tracker/util/io.hh>
#include <tracker/util/math.hh>

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

//__Track Parameter Type________________________________________________________________________
struct _track_parameters { fit_parameter t0, x0, y0, z0, vx, vy, vz; };
//----------------------------------------------------------------------------------------------

//__Fast Guess of Initial Track Parameters______________________________________________________
_track_parameters _guess_track(const full_event& points) {
  using namespace stat::type;
  using namespace stat::error;

  const auto& first = points.front();
  const auto& last = points.back();

  const uncertain_real first_t(first.t, first.width.t);
  const uncertain_real first_x(first.x, uniform(first.width.x));
  const uncertain_real first_y(first.y, uniform(first.width.y));
  const uncertain_real first_z(first.z, uniform(first.width.z));
  const auto dt = uncertain_real(last.t, last.width.t)          - first_t;
  const auto dx = uncertain_real(last.x, uniform(last.width.x)) - first_x;
  const auto dy = uncertain_real(last.y, uniform(last.width.y)) - first_y;
  const auto dz = uncertain_real(last.z, uniform(last.width.z)) - first_z;
  const auto vx = dx / dt;
  const auto vy = dy / dt;
  const auto vz = dz / dt;

  return {{first_t, first_t.error, 0, 0},
          {first_x, first_x.error, 0, 0},
          {first_y, first_y.error, 0, 0},
          {first_z, first_z.error, 0, 0},
          {vx,      vx.error,      0, 0},
          {vy,      vy.error,      0, 0},
          {vz,      vz.error,      0, 0}};
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
void _fit_event_minuit(const full_event& points,
                       const track::fit_settings& settings,
                       _track_parameters& parameters,
                       real_vector& covariance_matrix) {
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

  switch (settings.parameter_direction) {
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

  constexpr const std::size_t dimension = 6;
  Double_t matrix[dimension][dimension];
  minuit.mnemat(&matrix[0][0], dimension);
  covariance_matrix.clear();
  covariance_matrix.reserve(dimension * dimension);
  for (std::size_t i = 0; i < dimension; ++i)
    for (std::size_t j = 0; j < dimension; ++j)
      covariance_matrix.push_back(matrix[i][j]);
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

track::track(const std::vector<hit>& points,
             const track::fit_settings& settings)
    : track(add_width(points), settings) {}

//__Track Constructor___________________________________________________________________________
track::track(const std::vector<full_hit>& points,
             const track::fit_settings& settings)
    : _full_event(points), _settings(settings) {

  auto fit_track = _guess_track(_full_event);
  _fit_event_minuit(_full_event, _settings, fit_track, _covariance);

  _t0 = std::move(fit_track.t0);
  _x0 = std::move(fit_track.x0);
  _y0 = std::move(fit_track.y0);
  _z0 = std::move(fit_track.z0);
  _vx = std::move(fit_track.vx);
  _vy = std::move(fit_track.vy);
  _vz = std::move(fit_track.vz);

  const auto& full_event_begin = _full_event.cbegin();
  const auto& full_event_end = _full_event.cend();

  std::transform(full_event_begin, full_event_end, std::back_inserter(_delta_chi2),
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
const r4_point track::operator()(const real p) const {
  // TODO: specialize for each fit parameter
  const auto dt = (p - _z0.value) / _vz.value;
  return {dt + _t0.value, std::fma(dt, _vx.value, _x0.value), std::fma(dt, _vy.value, _y0.value), p};
}
//----------------------------------------------------------------------------------------------

//__Get Position of Track at Fixed T____________________________________________________________
const r4_point track::point(const real t) const {
  const auto dt = t - _t0.value;
  return {t,
          std::fma(dt, _vx.value, _x0.value),
          std::fma(dt, _vy.value, _y0.value),
          std::fma(dt, _vz.value, _z0.value)};
}
//----------------------------------------------------------------------------------------------

//__Get Error in Position of Track at Fixed T___________________________________________________
const r4_point track::point_error(const real t) const {
  // TODO: complex
  return {};
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter from Track________________________________________________________________
const fit_parameter track::fit_of(const track::parameter p) const {
  switch (p) {
    case track::parameter::T0: return _t0;
    case track::parameter::X0: return _x0;
    case track::parameter::Y0: return _y0;
    case track::parameter::Z0: return _z0;
    case track::parameter::VX: return _vx;
    case track::parameter::VY: return _vy;
    case track::parameter::VZ: return _vz;
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

//__Relativistic Beta for the Track_____________________________________________________________
real track::beta() const {
  return util::math::hypot(_vx.value, _vy.value, _vz.value) / units::speed_of_light;
}
//----------------------------------------------------------------------------------------------

//__Error in Relativistic Beta for the Track____________________________________________________
real track::beta_error() const {
  constexpr auto c_inv = 1.0L / units::speed_of_light;
  constexpr auto c_inv2 = c_inv * c_inv;
  const real_array<9> covariance{
    _covariance[6*3+3] * c_inv2, _covariance[6*3+4] * c_inv2, _covariance[6*3+5] * c_inv2,
    _covariance[6*4+3] * c_inv2, _covariance[6*4+4] * c_inv2, _covariance[6*4+5] * c_inv2,
    _covariance[6*5+3] * c_inv2, _covariance[6*5+4] * c_inv2, _covariance[6*5+5] * c_inv2};
  const real_array<3> gradient{
    2 * _vx.value * c_inv, 2 * _vy.value * c_inv, 2 * _vz.value * c_inv};
  return std::sqrt(stat::error::propagate(gradient, covariance));
}
//----------------------------------------------------------------------------------------------

//__Unit Vector along Track_____________________________________________________________________
const r3_point track::unit() const {
  return r3_point{_vx.value, _vy.value, _vz.value} / util::math::hypot(_vx.value, _vy.value, _vz.value);
}
//----------------------------------------------------------------------------------------------

//__Error in Unit Vector along Track____________________________________________________________
const r3_point track::unit_error() const {
  const auto base = 1.0L / std::pow(util::math::square(_vx.value, _vy.value, _vz.value), 1.5L);
  const real_array<9> covariance{
    _covariance[6*3+3], _covariance[6*3+4], _covariance[6*3+5],
    _covariance[6*4+3], _covariance[6*4+4], _covariance[6*4+5],
    _covariance[6*5+3], _covariance[6*5+4], _covariance[6*5+5]};
  const real_array<3> x_gradient{
    base * util::math::square(_vy.value, _vz.value),
    base * -_vx.value * _vy.value,
    base * -_vx.value * _vz.value};
  const real_array<3> y_gradient{
    base * -_vy.value * _vz.value,
    base * util::math::square(_vz.value, _vx.value),
    base * -_vy.value * _vz.value};
  const real_array<3> z_gradient{
    base * -_vz.value * _vx.value,
    base * -_vy.value * _vz.value,
    base * util::math::square(_vx.value, _vy.value)};
  return {
    std::sqrt(stat::error::propagate(x_gradient, covariance)),
    std::sqrt(stat::error::propagate(y_gradient, covariance)),
    std::sqrt(stat::error::propagate(z_gradient, covariance))
  };
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared Test Statistic__________________________________________________________________
real track::chi_squared() const {
  return std::accumulate(_delta_chi2.cbegin(), _delta_chi2.cend(), 0.0L);
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
  switch (_settings.parameter_direction) {
    case Coordinate::T:
      if (p == track::parameter::T0 || q == track::parameter::T0) {
        return 0;
      } else {
        const auto p_index = _shift_covariance_index(p)
                           - (p == track::parameter::X0 || p == track::parameter::Y0);
        const auto q_index = _shift_covariance_index(q)
                           - (q == track::parameter::X0 || q == track::parameter::Y0);
        return _covariance[6 * p_index + q_index];
      }
    case Coordinate::X:
      if (p == track::parameter::X0 || q == track::parameter::X0) {
        return 0;
      } else {
        const auto p_index = _shift_covariance_index(p) - (p == track::parameter::Y0);
        const auto q_index = _shift_covariance_index(q) - (q == track::parameter::Y0);
        return _covariance[6 * p_index + q_index];
      }
    case Coordinate::Y:
      if (p == track::parameter::Y0 || q == track::parameter::Y0) {
        return 0;
      } else {
        return _covariance[6 * _shift_covariance_index(p) + _shift_covariance_index(q)];
      }
    case Coordinate::Z:
      if (p == track::parameter::Z0 || q == track::parameter::Z0) {
        return 0;
      } else {
        return _covariance[6 * _shift_covariance_index(p) + _shift_covariance_index(q)];
      }
  }
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

//__Get Event from Track________________________________________________________________________
const std::vector<hit> track::event() const {
  std::vector<hit> out;
  out.reserve(_full_event.size());
  util::algorithm::back_insert_transform(_full_event, out,
    [](const auto& point) { return reduce_to_r4(point); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Track Output Stream Operator________________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const track& track) {
  static const std::string bar(80, '-');
  os << bar << "\n";
  os.precision(7);
  os << "* Track Parameters: \n"
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
  const auto detectors = track.detectors();
  const auto size = points.size();
  for (size_t i = 0; i < size; ++i)
    os << "      " << detectors[i] << " " << points[i] << "\n";
  os << "\n    back:  " << track.back()  << "\n";

  os.precision(7);
  os << "* Statistics: \n"
     << "    dof:      " << track.degrees_of_freedom()               << "\n"
     << "    chi2:     " << track.chi_squared() << " = ";
  util::io::print_range(track.chi_squared_vector(), " + ", "", os)   << "\n";
  os << "    chi2/dof: " << track.chi_squared_per_dof()              << "\n"
     << "    p-value:  " << stat::chi_squared_p_value(track)         << "\n"
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
     << "    beta:  " << track.beta()  << "  (+/- " << track.beta_error() << ")\n"
     << "    unit:  " << track.unit()  << "  (+/- " << track.unit_error() << ")\n";

  return os << bar;
}
//----------------------------------------------------------------------------------------------

//__Fit all Seeds to Tracks_____________________________________________________________________
template<class EventVector,
    typename = std::enable_if_t<is_r4_type_v<typename EventVector::value_type::value_type>>>
const track_vector fit_seeds(const EventVector& seeds,
                             const track::fit_settings& settings) {
  track_vector out;
  out.reserve(seeds.size());
  for (const auto& seed : seeds)
    out.emplace_back(seed, settings);
  return out;
}
const track_vector fit_seeds(const event_vector& seeds,
                             const track::fit_settings& settings) {
  return fit_seeds<>(seeds, settings);
}
const track_vector fit_seeds(const full_event_vector& seeds,
                             const track::fit_settings& settings) {
  return fit_seeds<>(seeds, settings);
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

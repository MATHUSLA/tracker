/*
 * src/tracker/analysis/vertex.cc
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

#include <tracker/analysis/vertex.hh>

#include <tracker/core/stat.hh>
#include <tracker/core/units.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/error.hh>
#include <tracker/util/io.hh>
#include <tracker/util/math.hh>

#include "../helper/analysis.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Calculate Distance from Vertex to Track_____________________________________________________
const stat::type::uncertain_real _vertex_track_r3_distance(const real t,
                                                           const real x,
                                                           const real y,
                                                           const real z,
                                                           const track& track) {
  const auto track_point = track.at_t(t);
  const auto dx = track_point.x - x;
  const auto dy = track_point.y - y;
  const auto dz = track_point.z - z;
  const auto total_dt = t - track.t0_value();
  const auto distance = util::math::hypot(dx, dy, dz);
  const auto inverse_distance = 1.0L / distance;
  const auto dx_by_D = inverse_distance * dx;
  const auto dy_by_D = inverse_distance * dy;
  const auto dz_by_D = inverse_distance * dz;
  const real_array<6UL> gradient{
    util::math::fused_product(track.vx_value(), dx_by_D,
                              track.vy_value(), dy_by_D,
                              track.vz_value(), dz_by_D),
    dx_by_D,
    dy_by_D,
    total_dt * dx_by_D,
    total_dt * dy_by_D,
    total_dt * dz_by_D};
  return stat::type::uncertain_real(
    distance,
    stat::error::propagate(gradient, track.covariance_matrix()));
}
//----------------------------------------------------------------------------------------------

//__Calculate Squared Residual of Vertex wrt Track______________________________________________
real _vertex_squared_residual(const stat::type::uncertain_real& distance) {
  return util::math::sum_squares(distance.value / distance.error);
}
//----------------------------------------------------------------------------------------------

//__Calculate Squared Residual of Vertex wrt Track______________________________________________
real _vertex_squared_residual(const real t,
                              const real x,
                              const real y,
                              const real z,
                              const track& track) {
  return _vertex_squared_residual(_vertex_track_r3_distance(t, x, y, z, track));
}
//----------------------------------------------------------------------------------------------

//__Fast Guess of Initial Track Parameters______________________________________________________
vertex::fit_parameters _guess_vertex(const track_vector& tracks) {
  const auto size = tracks.size();

  std::vector<full_hit> track_fronts;
  track_fronts.reserve(size);
  util::algorithm::back_insert_transform(tracks, track_fronts, [](const auto& track) {
    const auto full_front = track.full_front();
    const auto front_z = full_front.t;
    const auto point = track.at_z(front_z);
    const auto error = track.error_at_z(front_z);
    return full_hit{point.t, point.x, point.y, point.z,
             r4_point{error.t, error.x, error.y, stat::error::uniform(full_front.width.z)}};
  });

  real_vector t_errors, x_errors, y_errors, z_errors;
  t_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, t_errors,
    [](const auto& front) { return front.width.t; });
  x_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, x_errors,
    [](const auto& front) { return stat::error::uniform(front.width.x); });
  y_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, y_errors,
    [](const auto& front) { return stat::error::uniform(front.width.y); });
  z_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, z_errors,
    [](const auto& front) { return stat::error::uniform(front.width.z); });

  const auto average_point = std::accumulate(track_fronts.cbegin(), track_fronts.cend(), r4_point{},
    [](const auto sum, const auto& point) { return sum + reduce_to_r4(point); })
      / static_cast<real>(size);

  return {{average_point.t, stat::error::propagate_average(t_errors), 0, 0},
          {average_point.x, stat::error::propagate_average(x_errors), 0, 0},
          {average_point.y, stat::error::propagate_average(y_errors), 0, 0},
          {average_point.z, stat::error::propagate_average(z_errors), 0, 0}};
}
//----------------------------------------------------------------------------------------------

//__Gaussian Negative Log Likelihood Calculation________________________________________________
thread_local track_vector&& _nll_fit_tracks = {};
void _gaussian_nll(Int_t&, Double_t*, Double_t& out, Double_t* x, Int_t) {
  out = std::accumulate(_nll_fit_tracks.cbegin(), _nll_fit_tracks.cend(), 0.0L,
    [&](const auto sum, const auto& track) {
      const auto distance = _vertex_track_r3_distance(x[0], x[1], x[2], x[3], track);
      return sum + std::fma(0.5L, _vertex_squared_residual(distance), std::log(distance.error));
  });
}
//----------------------------------------------------------------------------------------------

//__MINUIT Gaussian Fitter______________________________________________________________________
void _fit_tracks_minuit(const track_vector& tracks,
                        vertex::fit_parameters& parameters,
                        vertex::covariance_matrix_type& covariance_matrix) {
  auto& t = parameters.t;
  auto& x = parameters.x;
  auto& y = parameters.y;
  auto& z = parameters.z;
  _nll_fit_tracks = tracks;

  TMinuit minuit;
  helper::minuit::initialize(minuit, "T", t, "X", x, "Y", y, "Z", z);
  helper::minuit::execute(minuit, _gaussian_nll);
  helper::minuit::get_parameters(minuit, t, x, y, z);
  helper::minuit::get_covariance<vertex::free_parameter_count>(minuit, covariance_matrix);
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Vertex Constructor__________________________________________________________________________
vertex::vertex(const track_vector& tracks) : _tracks(tracks) {
  if (_tracks.size() > 1) {
    _guess = _guess_vertex(_tracks);
    _final = _guess;
    _fit_tracks_minuit(_tracks, _final, _covariance);
    util::algorithm::back_insert_transform(_tracks, _delta_chi2,
      [&](const auto& track) {
        return _vertex_squared_residual(
          _final.t.value, _final.x.value, _final.y.value, _final.z.value, track); });
  } else {
    _delta_chi2.resize(1, 0);
    _guess = {};
    _final = {};
    _covariance = {};
  }
}
//----------------------------------------------------------------------------------------------

//__Vertex Point________________________________________________________________________________
const r4_point vertex::point() const {
  return {_final.t.value, _final.x.value, _final.y.value, _final.z.value};
}
//----------------------------------------------------------------------------------------------

//__Error in Calculation of Vertex Point________________________________________________________
const r4_point vertex::point_error() const {
  return {_final.t.error, _final.x.error, _final.y.error, _final.z.error};
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter from Vertex_______________________________________________________________
const fit_parameter vertex::fit_of(const vertex::parameter p) const {
  switch (p) {
    case vertex::parameter::T: return _final.t;
    case vertex::parameter::X: return _final.x;
    case vertex::parameter::Y: return _final.y;
    case vertex::parameter::Z: return _final.z;
  }
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter Value from Vertex_________________________________________________________
real vertex::value(const vertex::parameter p) const {
  return fit_of(p).value;
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter Error from Vertex_________________________________________________________
real vertex::error(const vertex::parameter p) const {
  return fit_of(p).error;
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared Test Statistic__________________________________________________________________
real vertex::chi_squared() const {
  return std::accumulate(_delta_chi2.cbegin(), _delta_chi2.cend(), 0.0L);
}
//----------------------------------------------------------------------------------------------

//__Vertex Degrees of Freedom___________________________________________________________________
size_t vertex::degrees_of_freedom() const {
  return 4UL;
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared per Degree of Freedom___________________________________________________________
real vertex::chi_squared_per_dof() const {
  return chi_squared() / degrees_of_freedom();
}
//----------------------------------------------------------------------------------------------

//__Get Variance of a Vertex Parameter__________________________________________________________
real vertex::variance(const vertex::parameter p) const {
  return covariance(p, p);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Get Shift Index of Vertex Parameters for Covariance Matrix__________________________________
constexpr std::size_t _shift_covariance_index(const vertex::parameter p) {
  switch (p) {
    case vertex::parameter::T: return 0;
    case vertex::parameter::X: return 1;
    case vertex::parameter::Y: return 2;
    case vertex::parameter::Z: return 3;
  }
}
//----------------------------------------------------------------------------------------------
} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Get Covariance between Vertex Parameters____________________________________________________
real vertex::covariance(const vertex::parameter p,
                        const vertex::parameter q) const {
  return _covariance[4 * _shift_covariance_index(p) + _shift_covariance_index(q)];
}
//----------------------------------------------------------------------------------------------

//__Fill Plots with Vertexing Variables_________________________________________________________
void vertex::fill_plots(plot::histogram_collection& collection,
                        const vertex::plotting_keys& keys) const {
  if (collection.count(keys.t)) collection[keys.t].insert(t_value() / units::time);
  if (collection.count(keys.x)) collection[keys.x].insert(x_value() / units::length);
  if (collection.count(keys.y)) collection[keys.y].insert(y_value() / units::length);
  if (collection.count(keys.z)) collection[keys.z].insert(z_value() / units::length);
  if (collection.count(keys.t_error)) collection[keys.t_error].insert(t_error() / units::time);
  if (collection.count(keys.x_error)) collection[keys.x_error].insert(x_error() / units::length);
  if (collection.count(keys.y_error)) collection[keys.y_error].insert(y_error() / units::length);
  if (collection.count(keys.z_error)) collection[keys.z_error].insert(z_error() / units::length);
  if (collection.count(keys.chi_squared_per_dof)) collection[keys.chi_squared_per_dof].insert(chi_squared_per_dof());
  if (collection.count(keys.size)) collection[keys.size].insert(size());
}
//----------------------------------------------------------------------------------------------

//__Vertex Output Stream Operator_______________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const vertex& vertex) {
  static const std::string bar(80, '-');
  os << bar << "\n";

  os << "* Vertex:\n"
     << "    " << vertex.point() << " (+/- " << vertex.point_error() << ")\n";

  os << "* Tracks: \n";
  for (const auto& track : vertex.tracks()) {
    os << "    (" << track.t0_value() << ", "
                  << track.x0_value() << ", "
                  << track.y0_value() << ", "
                  << track.z0_value() << ", "
                  << track.vx_value() << ", "
                  << track.vy_value() << ", "
                  << track.vz_value() << ")\n";
  }

  os.precision(7);
  os << "* Statistics: \n"
     << "    dof:      " << vertex.degrees_of_freedom()               << "\n"
     << "    chi2:     " << vertex.chi_squared() << " = ";
  util::io::print_range(vertex.chi_squared_vector(), " + ", "", os)   << "\n";
  os << "    chi2/dof: " << vertex.chi_squared_per_dof()              << "\n"
     << "    p-value:  " << stat::chi_squared_p_value(vertex)         << "\n"
     << "    cov mat:  | ";
  const auto matrix = vertex.covariance_matrix();
  for (size_t i = 0; i < 4; ++i) {
    if (i > 0) os << "              | ";
    for (size_t j = 0; j < 4; ++j) {
      const auto cell = matrix[4*i+j];
      if (i == j) {
        os << util::io::bold << util::io::underline
           << cell << util::io::reset_font << " ";
      } else {
        os << cell << " ";
      }
    }
    os << "|\n";
  }

  return os << bar;
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

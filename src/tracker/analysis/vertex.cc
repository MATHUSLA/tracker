/*
 * src/tracker/vertex.cc
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

#include <tracker/vertex.hh>

#include <ROOT/TMinuit.h>

#include <tracker/stat.hh>

#include <tracker/util/error.hh>
#include <tracker/util/io.hh>
#include <tracker/util/math.hh>

#include "analysis_helper.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

const stat::type::uncertain_real _vertex_track_r3_distance(const real t,
                                                           const real x,
                                                           const real y,
                                                           const real z,
                                                           const track& track) {
  const auto track_point = track.point(t);
  const auto dx = track_point.x - x;
  const auto dy = track_point.y - y;
  const auto dz = track_point.z - z;
  const auto total_dt = track.t0_value() - t;
  const auto distance = util::math::hypot(dx, dy, dz);
  const auto inverse_distance = 1.0L / distance;
  const auto dx_by_D = inverse_distance * dx;
  const auto dy_by_D = inverse_distance * dy;
  const auto dz_by_D = inverse_distance * dz;
  const real_array<6UL> gradient{
    -util::math::fused_product(track.vx_value(), dx_by_D, track.vy_value(), dy_by_D, track.vz_value(), dz_by_D),
    dx_by_D,
    dy_by_D,
    total_dt * dx_by_D,
    total_dt * dy_by_D,
    total_dt * dz_by_D};
  const auto covariance = to_array<36UL>(track.covariance_matrix());
  return stat::type::uncertain_real(distance, stat::error::propagate(gradient, covariance));
}

//__Calculate Squared Residual of Vertex wrt Track______________________________________________
real _vertex_squared_residual(const stat::type::uncertain_real& distance) {
  return util::math::square(distance.value / distance.error);
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

//__Vertex Parameter Type_______________________________________________________________________
struct _vertex_parameters { fit_parameter t, x, y, z; };
//----------------------------------------------------------------------------------------------

//__Fast Guess of Initial Track Parameters______________________________________________________
_vertex_parameters _guess_vertex(const track_vector& tracks) {
  // TODO: fix error propagation
  const auto average_point = std::accumulate(tracks.cbegin(), tracks.cend(), r4_point{},
    [](const auto& sum, const auto& track) { return sum + track.front(); })
    / static_cast<real>(tracks.size());
  return {{average_point.t, 1, 0, 0},
          {average_point.x, 1, 0, 0},
          {average_point.y, 1, 0, 0},
          {average_point.z, 1, 0, 0}};
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
                        _vertex_parameters& parameters,
                        real_vector& covariance_matrix) {
  _nll_fit_tracks = tracks;
  auto& t = parameters.t;
  auto& x = parameters.x;
  auto& y = parameters.y;
  auto& z = parameters.z;

  TMinuit minuit;
  helper::minuit::initialize(minuit, "T", t, "X", x, "Y", y, "Z", z);
  helper::minuit::execute(minuit, _gaussian_nll);
  helper::minuit::get_parameters(minuit, t, x, y, z);
  helper::minuit::get_covariance<4UL>(minuit, covariance_matrix);
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Vertex Constructor__________________________________________________________________________
vertex::vertex(const track_vector& tracks) : _tracks(tracks) {
  auto fit_vertex = _guess_vertex(_tracks);
  _fit_tracks_minuit(_tracks, fit_vertex, _covariance);

  _t = std::move(fit_vertex.t);
  _x = std::move(fit_vertex.x);
  _y = std::move(fit_vertex.y);
  _z = std::move(fit_vertex.z);

  const auto& tracks_begin = _tracks.cbegin();
  const auto& tracks_end = _tracks.cend();

  std::transform(tracks_begin, tracks_end, std::back_inserter(_delta_chi2),
    [&](const auto& track) {
      return _vertex_squared_residual(_t.value, _x.value, _y.value, _z.value, track); });
}
//----------------------------------------------------------------------------------------------

//__Vertex Point________________________________________________________________________________
const r4_point vertex::point() const {
  return {_t.value, _x.value, _y.value, _z.value};
}
//----------------------------------------------------------------------------------------------

//__Error in Calculation of Vertex Point________________________________________________________
const r4_point vertex::point_error() const {
  return {_t.error, _x.error, _y.error, _z.error};
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter from Vertex_______________________________________________________________
const fit_parameter vertex::fit_of(const vertex::parameter p) const {
  switch (p) {
    case vertex::parameter::T: return _t;
    case vertex::parameter::X: return _x;
    case vertex::parameter::Y: return _y;
    case vertex::parameter::Z: return _z;
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
  return -4 + std::accumulate(_tracks.cbegin(), _tracks.cend(), 0UL,
    [](const auto sum, const auto& track) { return sum + track.degrees_of_freedom(); });
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

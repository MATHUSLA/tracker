/*
 * include/tracker/analysis/track.hh
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

#ifndef TRACKER__ANALYSIS__TRACK_HH
#define TRACKER__ANALYSIS__TRACK_HH
#pragma once

#include <tracker/analysis/type.hh>
#include <tracker/geometry.hh>
#include <tracker/plot.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Track Object________________________________________________________________________________
class track {
public:
  enum class parameter { T0, X0, Y0, Z0, VX, VY, VZ };
  struct fit_parameters { fit_parameter t0, x0, y0, z0, vx, vy, vz; };

  static constexpr std::size_t free_parameter_count = 6UL;
  using covariance_matrix_type = real_array<free_parameter_count * free_parameter_count>;

  track(const event& points,
        const Coordinate direction=Coordinate::Z);

  track(const full_event& points,
        const Coordinate direction=Coordinate::Z);

  track(const track& rhs) = default;
  track(track&& rhs)      = default;
  track& operator=(const track& rhs) = default;
  track& operator=(track&& rhs)      = default;

  const r4_point operator()(const real p) const;

  const r4_point point(const real t) const;
  const r4_point point_error(const real t) const;

  const r4_point at_t(const real t) const;
  const r4_point at_x(const real x) const;
  const r4_point at_y(const real y) const;
  const r4_point at_z(const real z) const;
  const r4_point at(const Coordinate c,
                    const real r) const;

  const r4_point error_at_t(const real t) const;
  const r4_point error_at_x(const real x) const;
  const r4_point error_at_y(const real y) const;
  const r4_point error_at_z(const real z) const;
  const r4_point error_at(const Coordinate c,
                          const real r) const;

  const fit_parameters guess_fit() const { return _guess; }
  const fit_parameters final_fit() const { return _final; }

  const fit_parameter t0() const { return _final.t0; }
  const fit_parameter x0() const { return _final.x0; }
  const fit_parameter y0() const { return _final.y0; }
  const fit_parameter z0() const { return _final.z0; }
  const fit_parameter vx() const { return _final.vx; }
  const fit_parameter vy() const { return _final.vy; }
  const fit_parameter vz() const { return _final.vz; }
  const fit_parameter fit_of(const parameter p) const;

  real t0_value() const { return _final.t0.value; }
  real x0_value() const { return _final.x0.value; }
  real y0_value() const { return _final.y0.value; }
  real z0_value() const { return _final.z0.value; }
  real vx_value() const { return _final.vx.value; }
  real vy_value() const { return _final.vy.value; }
  real vz_value() const { return _final.vz.value; }
  real value(const parameter p) const;

  real t0_error() const { return _final.t0.error; }
  real x0_error() const { return _final.x0.error; }
  real y0_error() const { return _final.y0.error; }
  real z0_error() const { return _final.z0.error; }
  real vx_error() const { return _final.vx.error; }
  real vy_error() const { return _final.vy.error; }
  real vz_error() const { return _final.vz.error; }
  real error(const parameter p) const;

  real beta() const;
  real beta_error() const;
  const r3_point unit() const;
  const r3_point unit_error() const;

  real angle() const;
  real angle_error() const;

  real chi_squared() const;
  size_t degrees_of_freedom() const;
  real chi_squared_per_dof() const;
  const real_vector& chi_squared_vector() const { return _delta_chi2; }

  real variance(const parameter p) const;
  real covariance(const parameter p,
                  const parameter q) const;
  const covariance_matrix_type& covariance_matrix() const { return _covariance; }

  const hit front() const;
  const hit back() const;
  const analysis::event event() const;

  const full_hit full_front() const { return _full_event.front(); }
  const full_hit full_back() const { return _full_event.back(); }
  const analysis::full_event& full_event() const { return _full_event; }

  std::size_t size() const { return _full_event.size(); }
  bool empty() const { return _full_event.empty(); }

  const geometry::structure_vector& detectors() const { return _detectors; }

  Coordinate direction() const { return _direction; }

  struct plotting_keys {
    plot::histogram::name_type t0, x0, y0, z0, vx, vy, vz,
      t0_error, x0_error, y0_error, z0_error, vx_error, vy_error, vz_error,
      chi_squared_per_dof,
      size,
      beta, beta_error,
      angle, angle_error;
  };

  void fill_plots(plot::histogram_collection& collection,
                  const plotting_keys& keys) const;

  // TODO: void draw(plot::canvas& canvas) const;

private:
  fit_parameters _guess, _final;
  std::vector<full_hit> _full_event;
  real_vector _delta_chi2;
  covariance_matrix_type _covariance;
  geometry::structure_vector _detectors;
  Coordinate _direction;
};
//----------------------------------------------------------------------------------------------

//__Track Output Stream Operator________________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const track& track);
//----------------------------------------------------------------------------------------------

//__Vector of Tracks____________________________________________________________________________
using track_vector = std::vector<track>;
//----------------------------------------------------------------------------------------------

//__Fit all Seeds to Tracks_____________________________________________________________________
const track_vector independent_fit_seeds(const event_vector& seeds,
                                         const Coordinate direction=Coordinate::Z);
const track_vector independent_fit_seeds(const full_event_vector& seeds,
                                         const Coordinate direction=Coordinate::Z);
//----------------------------------------------------------------------------------------------

//__Fit all Seeds to Tracks using Overlaps______________________________________________________
const track_vector overlap_fit_seeds(const event_vector& seeds,
                                     const Coordinate direction=Coordinate::Z,
                                     const std::size_t min_overlap=2UL);
const track_vector overlap_fit_seeds(const full_event_vector& seeds,
                                     const Coordinate direction=Coordinate::Z,
                                     const std::size_t min_overlap=2UL);
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__TRACK_HH */

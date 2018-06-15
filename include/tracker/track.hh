/*
 * include/tracker/track.hh
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

#ifndef TRACKER__TRACK_HH
#define TRACKER__TRACK_HH
#pragma once

#include <tracker/analysis.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

using namespace type;

//__Track Object________________________________________________________________________________
class track {
public:
  enum class parameter { T0, X0, Y0, Z0, VX, VY, VZ };
  struct fit_settings {
    fit_settings() {}
    Coordinate          parameter_direction = Coordinate::Z;
    std::string         command_name        = "MIGRAD";
    std::vector<double> command_parameters  = {};
    bool                graphics_on         = false;
    integer             print_level         = -1;
    double              error_def           = 0.5;
    integer             max_iterations      = 300;
  };

  track(const event& points,
        const fit_settings& settings={});

  track(const full_event& points,
        const fit_settings& settings={});

  track(const track& rhs) = default;
  track(track&& rhs)      = default;
  track& operator=(const track& rhs) = default;
  track& operator=(track&& rhs)      = default;

  const hit operator()(const real z) const;

  const fit_parameter t0() const { return _t0; }
  const fit_parameter x0() const { return _x0; }
  const fit_parameter y0() const { return _y0; }
  const fit_parameter z0() const { return _z0; }
  const fit_parameter vx() const { return _vx; }
  const fit_parameter vy() const { return _vy; }
  const fit_parameter vz() const { return _vz; }
  const fit_parameter fit_of(const parameter p) const;

  real t0_value() const { return _t0.value; }
  real x0_value() const { return _x0.value; }
  real y0_value() const { return _y0.value; }
  real z0_value() const { return _z0.value; }
  real vx_value() const { return _vx.value; }
  real vy_value() const { return _vy.value; }
  real vz_value() const { return _vz.value; }
  real value(const parameter p) const;

  real t0_error() const { return _t0.error; }
  real x0_error() const { return _x0.error; }
  real y0_error() const { return _y0.error; }
  real z0_error() const { return _z0.error; }
  real vx_error() const { return _vx.error; }
  real vy_error() const { return _vy.error; }
  real vz_error() const { return _vz.error; }
  real error(const parameter p) const;

  real beta() const;
  real beta_error() const;
  const r3_point unit() const;
  const r3_point unit_error() const;

  real chi_squared() const;
  size_t degrees_of_freedom() const;
  real chi_squared_per_dof() const;
  const real_vector& chi_squared_vector() const { return _delta_chi2; }

  real variance(const parameter p) const;
  real covariance(const parameter p,
                  const parameter q) const;
  const real_vector& covariance_matrix() const { return _covariance; }

  const hit front() const;
  const hit back() const;
  const std::vector<hit> event() const;
  const std::vector<full_hit>& full_event() const { return _full_event; }
  const fit_settings& settings() const { return _settings; }
  const std::vector<std::string>& detectors() const { return _detectors; }

private:
  fit_parameter _t0, _x0, _y0, _z0, _vx, _vy, _vz;
  std::vector<full_hit> _full_event;
  real_vector _delta_chi2, _covariance;
  std::vector<std::string> _detectors;
  fit_settings _settings;
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
const track_vector fit_seeds(const event_vector& seeds,
                             const track::fit_settings& settings={});
const track_vector fit_seeds(const full_event_vector& seeds,
                             const track::fit_settings& settings={});
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__TRACK_HH */

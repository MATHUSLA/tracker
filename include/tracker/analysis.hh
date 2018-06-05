/*
 * include/analysis.hh
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

#ifndef TRACKER__ANALYSIS_HH
#define TRACKER__ANALYSIS_HH
#pragma once

#include <tracker/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

using namespace type;

//__Event Types_________________________________________________________________________________
using hit = r4_point;
using event = std::vector<hit>;
using event_vector = std::vector<event>;
//----------------------------------------------------------------------------------------------

//__Extended Event Types________________________________________________________________________
struct full_hit { real t, x, y, z; r4_point error; };
using full_event = std::vector<full_hit>;
using full_event_vector = std::vector<full_event>;
//----------------------------------------------------------------------------------------------

//__Full Hit Stream Operator Overload___________________________________________________________
inline std::ostream& operator<<(std::ostream& os,
                                const full_hit& point) {
  return os << "[(" << point.t << ", " << point.x << ", " << point.y << ", " << point.z
            << ") +/- " << point.error << "]";
}
//----------------------------------------------------------------------------------------------

//__Full Hit Equality___________________________________________________________________________
inline constexpr bool operator==(const full_hit& left,
                                 const full_hit& right) {
  return left.t == right.t
      && left.x == right.x
      && left.y == right.y
      && left.z == right.z
      && left.error == right.error;
}
//----------------------------------------------------------------------------------------------

//__Center Events by Coordinate_________________________________________________________________
const event centralize(const event& points,
                       const Coordinate coordinate);
const full_event centralize(const full_event& points,
                            const Coordinate coordinate);
//----------------------------------------------------------------------------------------------

//__Collapse Points by R4 Interval______________________________________________________________
const event collapse(const event& points,
                     const r4_point& ds);
const full_event collapse(const full_event& points,
                          const r4_point& ds);
//----------------------------------------------------------------------------------------------

//__Event Partition Type________________________________________________________________________
struct event_partition { event_vector parts; Coordinate coordinate; real interval; };
struct full_event_partition { full_event_vector parts; Coordinate coordinate; real interval; };
//----------------------------------------------------------------------------------------------

//__Partition Points by Coordinate______________________________________________________________
const event_partition partition(const event& points,
                                const Coordinate coordinate,
                                const real interval);
const full_event_partition partition(const full_event& points,
                                     const Coordinate coordinate,
                                     const real interval);
//----------------------------------------------------------------------------------------------

//__Fast Check if Points Form a Line____________________________________________________________
bool fast_line_check(const event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2);
bool fast_line_check(const full_event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2);
bool fast_line_check(const event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2,
                     const Coordinate x3);
bool fast_line_check(const full_event& points,
                     const real threshold,
                     const Coordinate x1,
                     const Coordinate x2,
                     const Coordinate x3);
//----------------------------------------------------------------------------------------------

//__Seeding Algorithm___________________________________________________________________________
const event_vector seed(const size_t n,
                        const event_partition& partition,
                        const real line_threshold);
const full_event_vector seed(const size_t n,
                             const full_event_partition& partition,
                             const real line_threshold);
//----------------------------------------------------------------------------------------------

//__Check if Seeds can be Joined________________________________________________________________
bool seeds_compatible(const event& first,
                      const event& second,
                      const size_t difference);
bool seeds_compatible(const full_event& first,
                      const full_event& second,
                      const size_t difference);
//----------------------------------------------------------------------------------------------

//__Join Two Seeds______________________________________________________________________________
const event join(const event& first,
                 const event& second,
                 const size_t difference);
const full_event join(const full_event& first,
                      const full_event& second,
                      const size_t difference);
//----------------------------------------------------------------------------------------------

//__Seed Join___________________________________________________________________________________
const event_vector join_all(const event_vector& seeds);
const full_event_vector join_all(const full_event_vector& seeds);
//----------------------------------------------------------------------------------------------

//__Fitting Parameter Type______________________________________________________________________
struct fit_parameter { real value, error, min, max; };
//----------------------------------------------------------------------------------------------

//__Fit Settings Type with Default Values_______________________________________________________
struct fit_settings {
  std::string         command_name       = "MIGRAD";
  std::vector<double> command_parameters = {};
  bool                graphics_on        = false;
  integer             print_level        = -1;
  double              error_def          = 0.5;
  integer             max_iterations     = 250;
};
//----------------------------------------------------------------------------------------------

//__Track Object________________________________________________________________________________
class track {
public:
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

  real t0_value() const { return _t0.value; }
  real x0_value() const { return _x0.value; }
  real y0_value() const { return _y0.value; }
  real z0_value() const { return _z0.value; }
  real vx_value() const { return _vx.value; }
  real vy_value() const { return _vy.value; }
  real vz_value() const { return _vz.value; }

  real t0_error() const { return _t0.error; }
  real x0_error() const { return _x0.error; }
  real y0_error() const { return _y0.error; }
  real z0_error() const { return _z0.error; }
  real vx_error() const { return _vx.error; }
  real vy_error() const { return _vy.error; }
  real vz_error() const { return _vz.error; }

  real residual() const;
  real squared_residual() const;
  const real_vector residual_vector() const;
  const real_vector& squared_residual_vector() const { return _squared_residuals; }

  real beta() const;

  real chi_squared() const;
  size_t degrees_of_freedom() const;
  real chi_squared_per_dof() const;
  const real_vector& chi_squared_vector() const { return _delta_chi_squared; }

  const event event() const;
  const full_event& full_event() const { return _full_event; }
  const fit_settings& settings() const { return _settings; }
  const std::vector<std::string>& detectors() const { return _detectors; }

  const hit front() const;
  const hit back() const;

private:
  fit_parameter _t0, _x0, _y0, _z0, _vx, _vy, _vz;
  std::vector<full_hit> _full_event;
  real_vector _squared_residuals, _delta_chi_squared;
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
                             const fit_settings& settings={});
const track_vector fit_seeds(const full_event_vector& seeds,
                             const fit_settings& settings={});
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS_HH */

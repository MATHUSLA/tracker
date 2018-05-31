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

#include "point.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

using namespace type;

//__Event Types_________________________________________________________________________________
using event_points = r4_point_vector;
using event_vector = std::vector<event_points>;
//----------------------------------------------------------------------------------------------

//__Average Point_______________________________________________________________________________
const r4_point mean(const event_points& points);
//----------------------------------------------------------------------------------------------

//__Time Normalize Events_______________________________________________________________________
const event_points time_normalize(const event_points& event);
//----------------------------------------------------------------------------------------------

//__Collapse Points by R4 Interval______________________________________________________________
const event_points collapse(const event_points& event,
                            const r4_point& ds);
//----------------------------------------------------------------------------------------------

//__Event Partition Type________________________________________________________________________
struct event_partition { event_vector parts; Coordinate coordinate; };
//----------------------------------------------------------------------------------------------

//__Partition Points by Coordinate______________________________________________________________
const event_partition partition(const event_points& points,
                                const real interval,
                                const Coordinate coordinate=Coordinate::Z);
//----------------------------------------------------------------------------------------------

//__Center of Geometric Object for each Point___________________________________________________
const event_points find_centers(const event_points& points);
//----------------------------------------------------------------------------------------------

//__Fast Check if Points Form a Line____________________________________________________________
bool fast_line_check(const event_points& points,
                     const real threshold);
//----------------------------------------------------------------------------------------------

//__Seeding Algorithm___________________________________________________________________________
const event_vector seed(const size_t n,
                        const event_points& event,
                        const r4_point& collapse_ds,
                        const real layer_dz,
                        const real line_dr);
//----------------------------------------------------------------------------------------------

//__Check if Seeds can be Joined________________________________________________________________
bool seeds_compatible(const event_points& first,
                      const event_points& second,
                      const size_t difference);
//----------------------------------------------------------------------------------------------

//__Join Two Seeds______________________________________________________________________________
const event_points join(const event_points& first,
                        const event_points& second,
                        const size_t difference);
//----------------------------------------------------------------------------------------------

//__Seed Join___________________________________________________________________________________
const event_vector join_all(const event_vector& seeds);
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
  track(const event_points& event,
        const fit_settings& settings={});

  track(const track& rhs) = default;
  track(track&& rhs)      = default;
  track& operator=(const track& rhs) = default;
  track& operator=(track&& rhs)      = default;

  const r4_point operator()(const real z) const;

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

  const event_points& event() const { return _event; }
  const fit_settings& settings() const { return _settings; }
  const std::vector<std::string>& detectors() const { return _detectors; }

private:
  fit_parameter _t0, _x0, _y0, _z0, _vx, _vy, _vz;
  event_points _event;
  real_vector _squared_residuals, _delta_chi_squared;
  std::vector<std::string> _detectors;
  fit_settings _settings;
};
//----------------------------------------------------------------------------------------------

//__Output Stream Operator______________________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const track& t);
//----------------------------------------------------------------------------------------------

//__Vector of Tracks____________________________________________________________________________
using track_vector = std::vector<track>;
//----------------------------------------------------------------------------------------------

//__Add Track from Seed to Track Vector_________________________________________________________
track_vector& operator+=(track_vector& tracks,
                         const event_points& seed);
//----------------------------------------------------------------------------------------------

//__Fit all Seeds to Tracks_____________________________________________________________________
const track_vector fit_seeds(const event_vector& seeds,
                             const fit_settings& settings={});
//----------------------------------------------------------------------------------------------

/*
//__Perform Full Analysis of Event______________________________________________________________
const track_vector full(const event_points& event,
                        const size_t max_seed_size,
                        const r4_point& max_collapse_ds,
                        const real max_layer_dz,
                        const real max_line_dr,
                        const fit_settings& settings={});
//----------------------------------------------------------------------------------------------
*/

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS_HH */

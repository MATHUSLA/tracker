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

//__Event Partition Type________________________________________________________________________
struct event_partition { event_vector parts; Coordinate coordinate; };
//----------------------------------------------------------------------------------------------

//__Collapse Points by R4 Interval______________________________________________________________
event_points collapse(const event_points& event,
                      const r4_point& ds);
//----------------------------------------------------------------------------------------------

//__Partition Points by Coordinate______________________________________________________________
event_partition partition(const event_points& points,
                          const real interval,
                          const Coordinate coordinate=Coordinate::Z);
//----------------------------------------------------------------------------------------------

//__Seeding Algorithm___________________________________________________________________________
event_vector seed(const size_t n,
                  const event_points& event,
                  const r4_point& collapse_ds,
                  const real layer_dz,
                  const real line_dr);
//----------------------------------------------------------------------------------------------

//__Seed Search_________________________________________________________________________________
event_vector seeds_with(const r4_point& point,
                        const event_vector& seeds);
//----------------------------------------------------------------------------------------------

//__Seed Search Front___________________________________________________________________________
event_vector seeds_starting_with(const r4_point& point,
                                 const event_vector& seeds);
//----------------------------------------------------------------------------------------------

//__Seed Search Back____________________________________________________________________________
event_vector seeds_ending_with(const r4_point& point,
                               const event_vector& seeds);
//----------------------------------------------------------------------------------------------

//__Seed Merge__________________________________________________________________________________
event_vector merge(const event_vector& seeds);
//----------------------------------------------------------------------------------------------

//__Fitting Parameter Types_____________________________________________________________________
struct fit_parameter { std::string name; real value, error, min, max; };
using fit_parameter_vector = std::vector<fit_parameter>;

struct track_parameters {
  fit_parameter t0, x0, y0, z0, vx, vy, vz, mean, variance;
};
//----------------------------------------------------------------------------------------------

//__Fit Settings Type with Default Values_______________________________________________________
struct fit_settings {
  double              error_def          = 0.5;
  integer             max_iterations     = 500;
  std::string         command_name       = "MIGRAD";
  std::vector<double> command_parameters = {};
};
//----------------------------------------------------------------------------------------------

//__Perform Gaussian Fit to Events______________________________________________________________
void fit_event(const event_points& events,
               fit_parameter_vector& parameters,
               const fit_settings& settings=fit_settings{});
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS_HH */

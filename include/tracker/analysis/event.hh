/*
 * include/tracker/analysis/event.hh
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

#ifndef TRACKER__ANALYSIS__EVENT_HH
#define TRACKER__ANALYSIS__EVENT_HH
#pragma once

#include <tracker/analysis/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Calculate Number of Hits per unit Length____________________________________________________
const r4_point event_density(const event& points);
const r4_point event_density(const full_event& points);
//----------------------------------------------------------------------------------------------

//__Calculate Number of Hits per unit Volume____________________________________________________
// TODO: real event_volume_density(const event& points);
// TODO: real event_volume_density(const full_event& points);
//----------------------------------------------------------------------------------------------

//__Calculate Number of Hits per Geometric Element______________________________________________
real geometric_event_density(const event& points);
real geometric_event_density(const full_event& points);
//----------------------------------------------------------------------------------------------

//__Find The Errors Associated with a Hit from Geometry_________________________________________
const full_hit add_width(const hit& point);
const full_event add_width(const event& points);
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

//__Reset Partition by new Interval_____________________________________________________________
const event_partition repartition(const event_partition& previous,
                                  const real interval);
const full_event_partition repartition(const full_event_partition& previous,
                                       const real interval);
//----------------------------------------------------------------------------------------------

//__Reduce Event Partition to Events____________________________________________________________
const event reduce_partition(const event_partition& previous);
const full_event reduce_partition(const full_event_partition& previous);
//----------------------------------------------------------------------------------------------

//__Calculate Density of Partition______________________________________________________________
// TODO: real partition_density(const event_partition& points);
// TODO: real partition_density(const full_event_partition& points);
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__EVENT_HH */

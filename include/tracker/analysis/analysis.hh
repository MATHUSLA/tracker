/*
 * include/tracker/analysis/analysis.hh
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

#ifndef TRACKER__ANALYSIS__ANALYSIS_HH
#define TRACKER__ANALYSIS__ANALYSIS_HH
#pragma once

#include <tracker/analysis/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Calculate Number of Hits per unit Length____________________________________________________
const r4_point event_density(const event& points);
const r4_point event_density(const full_event& points);
//----------------------------------------------------------------------------------------------

//__Calculate Number of Hits per unit Volume____________________________________________________
// real event_volume_density(const event& points);
// real event_volume_density(const full_event& points);
//----------------------------------------------------------------------------------------------

//__Calculate Number of Hits per Geometric Element______________________________________________
real geometric_event_density(const event& points);
real geometric_event_density(const full_event& points);
//----------------------------------------------------------------------------------------------

//__Find The Errors Associated with a Hit from Geometry_________________________________________
const full_hit add_width(const hit& point);
const full_event add_width(const event& points);
//----------------------------------------------------------------------------------------------

//__Center Events by Coordinate_________________________________________________________________
const event centralize(const event& points,
                       const Coordinate coordinate);
const full_event centralize(const full_event& points,
                            const Coordinate coordinate);
//----------------------------------------------------------------------------------------------

//__Compress Points by R4 Interval______________________________________________________________
const event compress(const event& points,
                     const r4_point& ds);
const full_event compress(const full_event& points,
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

//__Reset Partition by new Interval_____________________________________________________________
const event_partition repartition(const event_partition& previous,
                                  const real interval);
const full_event_partition repartition(const full_event_partition& previous,
                                       const real interval);
//----------------------------------------------------------------------------------------------

//__Calculate Density of Partition______________________________________________________________
// real partition_density(const event_partition& points);
// real partition_density(const full_event_partition& points);
//----------------------------------------------------------------------------------------------

//__Fast Check if Points Form a Line____________________________________________________________
bool is_linear(const event& points,
               const real threshold,
               const Coordinate x1,
               const Coordinate x2);
bool is_linear(const full_event& points,
               const real threshold,
               const Coordinate x1,
               const Coordinate x2);
bool is_linear(const event& points,
               const real threshold,
               const Coordinate x1,
               const Coordinate x2,
               const Coordinate x3);
bool is_linear(const full_event& points,
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

//__Join Two Seeds in Sequence__________________________________________________________________
const event sequential_join(const event& first,
                            const event& second,
                            const size_t difference);
const full_event sequential_join(const full_event& first,
                                 const full_event& second,
                                 const size_t difference);
//----------------------------------------------------------------------------------------------

//__Join Two Seeds Such That One is a Subset of the Other_______________________________________
const event subset_join(const event& first,
                        const event& second);
const full_event subset_join(const full_event& first,
                             const full_event& second);
//----------------------------------------------------------------------------------------------

//__Join Two Seeds Which form a Loop____________________________________________________________
const event loop_join(const event& first,
                      const event& second);
const full_event loop_join(const full_event& first,
                           const full_event& second);
//----------------------------------------------------------------------------------------------

//__Optimally Join All Seeds by Sequence________________________________________________________
const event_vector sequential_join_all(const event_vector& seeds);
const full_event_vector sequential_join_all(const full_event_vector& seeds);
//----------------------------------------------------------------------------------------------

//__Optimally Join All Seeds by Subset__________________________________________________________
const event_vector subset_join_all(const event_vector& seeds);
const full_event_vector subset_join_all(const full_event_vector& seeds);
//----------------------------------------------------------------------------------------------

//__Optimally Join All Seeds by Loop____________________________________________________________
const event_vector loop_join_all(const event_vector& seeds);
const full_event_vector loop_join_all(const full_event_vector& seeds);
//----------------------------------------------------------------------------------------------

//__Seed Join___________________________________________________________________________________
const event_vector join_all(const event_vector& seeds);
const full_event_vector join_all(const full_event_vector& seeds);
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__ANALYSIS_HH */

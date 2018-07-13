/*
 * include/tracker/analysis/seed.hh
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

#ifndef TRACKER__ANALYSIS__SEED_HH
#define TRACKER__ANALYSIS__SEED_HH
#pragma once

#include <tracker/analysis/event.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

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

//__Check if Points are Monotonic in One Coordinate_____________________________________________
bool is_monotonic(const event& points,
                  const Coordinate c);
bool is_monotonic(const full_event& points,
                  const Coordinate c);
//----------------------------------------------------------------------------------------------

//__Cylinder Topology Type______________________________________________________________________
struct cylinder {
  cylinder(const real radius);
  bool operator()(const event& points);
  bool operator()(const full_event& points);
};
//----------------------------------------------------------------------------------------------

//__Double Cone Topology Type___________________________________________________________________
struct double_cone {
  double_cone(const real radius);
  bool operator()(const event& points);
  bool operator()(const full_event& points);
};
//----------------------------------------------------------------------------------------------

//__Seeding Algorithm___________________________________________________________________________
const event_vector seed(const std::size_t n,
                        const event_partition& partition,
                        const real line_threshold);
const full_event_vector seed(const std::size_t n,
                             const full_event_partition& partition,
                             const real line_threshold);
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__SEED_HH */

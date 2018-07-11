/*
 * demo/prototype/geometry.hh
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

#ifndef TRACKER__PROTOTYPE__GEOMETRY_HH
#define TRACKER__PROTOTYPE__GEOMETRY_HH
#pragma once

#include <tracker/analysis/type.hh>
#include <tracker/geometry.hh>

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace geometry = MATHUSLA::TRACKER::geometry;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

//__Total Geometry of the Prototype Detector____________________________________________________
const geometry::structure_vector prototype_geometry();
//----------------------------------------------------------------------------------------------

//__Hits per Total Geometry_____________________________________________________________________
type::real modified_geometry_event_density(const analysis::event& event);
//----------------------------------------------------------------------------------------------

//__Combine Pair of Hits if they Occur in Overlapping RPCs______________________________________
const geometry::box_volume combine_rpc_volume_pair(const geometry::box_volume& first,
                                                   const geometry::box_volume& second);
//----------------------------------------------------------------------------------------------

//__Check If RPC Combine Created a Valid Strip Overlap__________________________________________
bool was_combine_successful(const geometry::box_volume& combined);
//----------------------------------------------------------------------------------------------

//__Construct True Hit from RPC Hit Volumes_____________________________________________________
const analysis::hit construct_hit(const type::real top_time,
                                  const type::real bottom_time,
                                  const geometry::box_volume& combined);
//----------------------------------------------------------------------------------------------

//__Combine All Hits that Occur in Overlapping RPCs_____________________________________________
const analysis::event combine_rpc_hits(const analysis::event& points,
                                       analysis::event& combined_rpc_hits,
                                       analysis::full_event& original_rpc_hits);
//----------------------------------------------------------------------------------------------

//__Reset Seed Vector Using RPC Combination Hits________________________________________________
const analysis::full_event_vector reset_seeds(const analysis::event_vector& joined_seeds,
                                              const analysis::event& combined_rpc_hits,
                                              const analysis::full_event& original_rpc_hits);
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

#endif /* TRACKER__PROTOTYPE__GEOMETRY_HH */

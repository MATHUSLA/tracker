/*
 * include/tracker/analysis/join.hh
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

#ifndef TRACKER__ANALYSIS__JOIN_HH
#define TRACKER__ANALYSIS__JOIN_HH
#pragma once

#include <tracker/analysis/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Join Two Seeds in Sequence__________________________________________________________________
const event sequential_join(const event& first,
                            const event& second,
                            const std::size_t difference);
const full_event sequential_join(const full_event& first,
                                 const full_event& second,
                                 const std::size_t difference);
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

#endif /* TRACKER__ANALYSIS__JOIN_HH */

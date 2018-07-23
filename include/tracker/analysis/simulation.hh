/*
 * include/tracker/analysis/simulation.hh
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

#ifndef TRACKER__ANALYSIS__SIMULATION_HH
#define TRACKER__ANALYSIS__SIMULATION_HH
#pragma once

#include <tracker/analysis/event.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { namespace simulation { ////////////////////////////////////////////////////

//__Compress Points by R4 Interval______________________________________________________________
const event compress(const event& points);
const full_event compress(const full_event& points);
//----------------------------------------------------------------------------------------------

//__Time Smear Points by Detector Time Resolution_______________________________________________
const event time_smear(const event& points);
const full_event time_smear(const full_event& points);
//----------------------------------------------------------------------------------------------

//__Simulate the Detector Efficiency____________________________________________________________
const event use_efficiency(const event& points,
                           const real efficiency);
const full_event use_efficiency(const full_event& points,
                                const real efficiency);
//----------------------------------------------------------------------------------------------

//__Add Noise to Points_________________________________________________________________________
const event add_noise(const event& points);
const full_event add_noise(const full_event& points);
//----------------------------------------------------------------------------------------------

} } /* namespace analysis::simulation */ ///////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__SIMULATION_HH */

/*
 * demo/module_box/tracking.cc
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

#include <tracker/analysis.hh>

#include <tracker/core/stat.hh>
#include <tracker/core/units.hh>
#include <tracker/geometry.hh>
#include <tracker/plot.hh>
#include <tracker/reader.hh>

#include "geometry.hh"

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace geometry = MATHUSLA::TRACKER::geometry;
namespace plot     = MATHUSLA::TRACKER::plot;
namespace stat     = MATHUSLA::TRACKER::stat;
namespace reader   = MATHUSLA::TRACKER::reader;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

//__Module Box Tracking Algorithm_______________________________________________________________
int module_box_tracking(int argc,
                        char* argv[]) {
  return 1;
}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

//__Main Function: ModuleBox Tracker____________________________________________________________
int main(int argc,
         char* argv[]) {
  return MATHUSLA::module_box_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------
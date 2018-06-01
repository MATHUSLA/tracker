/*
 * demo/prototype_tracking.cc
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
#include <tracker/geometry.hh>
#include <tracker/reader.hh>
#include <tracker/plot.hh>

namespace MATHUSLA {

//__Combine Hits if they Occur in Overlapping RPCs______________________________________________
const TRACKER::analysis::event_points combine_rpc_hits(const TRACKER::analysis::event_points& event) {
  using namespace TRACKER;
  
  return event;
}
//----------------------------------------------------------------------------------------------

//__Main Function: Prototype Tracker____________________________________________________________
int main(int argc, char* argv[]) {
  using namespace TRACKER;

  const auto options = reader::parse_input(argc, argv);
  const auto detector_map = reader::root::import_detector_map(options.geometry_map_file);

  plot::init();
  geometry::open(options.geometry_file);
  for (const auto& path : reader::root::search_directory(options.root_directory)) {
    for (const auto& event : reader::root::import_events(path, options, detector_map)) {
      plot::canvas canvas(path);

      for (const auto& name : geometry::full_structure_except({"world", "Sandstone", "Marl", "Mix", "Earth"})) {
        const auto limits = geometry::limits_of(name);
        canvas.add_point(limits.center, 0.25, plot::color::MAGENTA);
      }
      canvas.add_points(event, 1.5, {90, 90, 90});




      canvas.draw();
    }
  }

  geometry::close();
  plot::end();

  return 0;
}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

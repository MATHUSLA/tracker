/*
 * src/tracker.cc
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

#include "analysis.hh"
#include "geometry.hh"
#include "plot.hh"
#include "reader.hh"

#include "util/io.hh"

//__Main Function: Tracker______________________________________________________________________
int main(int argc, char* argv[]) {
  using namespace MATHUSLA;
  using namespace MATHUSLA::TRACKER;

  const auto options = reader::parse_input(argc, argv);

  plot::init();
  geometry::open(options.geometry_file);
  const auto detector_map = reader::root::import_detector_map(options.geometry_map_file);

  const auto paths = reader::root::search_directory(options.root_directory);
  std::cout << "File Count: " << paths.size() << "\n";
  for (const auto& path : paths) {
    std::cout << path << "\n";

    const auto events = reader::root::import_events(path, options, detector_map);

    std::cout << "Processing " << events.size() << " Events\n";

    for (const auto& unsorted_event : events) {
      plot::canvas canvas(path);

      // demo for Prototype only
      for (const auto& name : geometry::full_structure_except({"world", "Sandstone", "Marl", "Mix", "Earth"})) {
        const auto limits = geometry::limits_of(name);
        canvas.add_point(limits.center, 0.25, plot::color::MAGENTA);
      }

      const auto event           = analysis::time_normalize(unsorted_event);
      const auto collapsed_event = analysis::collapse(event, options.collapse_size);
      const auto layered_event   = analysis::partition(collapsed_event, options.layer_depth, options.layer_axis).parts;

      canvas.add_points(event, 1.5, {90, 90, 90});

      util::io::print_range(event, "\n", "OLD ") << "\n\n";
      util::io::print_range(collapsed_event, "\n", "NEW ") << "\n\n";
      for (const auto& layer : layered_event)
        util::io::print_range(layer, "\n", "LAYER ") << "\n\n\n";

      const auto seeds = analysis::seed(options.seed_size, unsorted_event, options.collapse_size, options.layer_depth, options.line_width);
      std::cout << "seeds (" << seeds.size() << "):\n\n";
      for (const auto& seed : seeds)
        util::io::print_range(seed) << "\n";

      util::io::newline();

      auto joined = analysis::join_all(seeds);

      std::cout << "joined (" << joined.size() << "):\n\n";
      for (const auto& seed : joined)
        util::io::print_range(seed) << "\n";

      const auto tracks = analysis::fit_seeds(joined, {"MIGRAD", {}, false, -1, 0.5, 150});
      std::cout << "\nTRACK FITTING (" << tracks.size() << "):\n\n";
      for (const auto& track : tracks) {
        const auto& event = track.event();
        canvas.add_points(event, 0.5);
        for (const auto& point : event) {
          const auto limits = geometry::limits_of_volume(point);
          canvas.add_point(limits.center, 0.25, plot::color::BLUE);
          canvas.add_box(limits.min, limits.max, 2, plot::color::BLUE);
        }
        canvas.add_line(event.front(), event.back());
        std::cout << track << "\n";
        canvas.add_line(track(event.front().z), track(event.back().z), 1, plot::color::RED);
      }

      canvas.draw();
    }
  }

  geometry::close();
  plot::end();

  std::cout << "Done!\n";

  return 0;
}
//----------------------------------------------------------------------------------------------

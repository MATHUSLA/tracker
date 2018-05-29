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
#include "units.hh"

#include "util/command_line_parser.hh"
#include "util/error.hh"
#include "util/io.hh"
#include "util/string.hh"

//__Missing Path Exit Command___________________________________________________________________
void exit_on_missing_path(const std::string& path,
                          const std::string& name) {
  using namespace MATHUSLA::util;
  error::exit_when(!io::path_exists(path),
    "[FATAL ERROR] ", name, " Missing: The file ", path, " cannot be found.\n");
}
//----------------------------------------------------------------------------------------------

//__Main Function: Tracker______________________________________________________________________
int main(int argc, char* argv[]) {
  using namespace MATHUSLA;
  using namespace MATHUSLA::TRACKER;

  using util::cli::option;

  option help_opt   ('h', "help",     "MATHUSLA Tracking Algorithm", option::no_arguments);
  option geo_opt    ('g', "geometry", "Geometry Import",             option::required_arguments);
  option map_opt    ('m', "map",      "Detector Map",                option::required_arguments);
  option root_opt   ('d', "data",     "ROOT Data Directory",         option::required_arguments);
  option script_opt ('s', "script",   "Tracking Script",             option::required_arguments);
  option quiet_opt  ('q', "quiet",    "Quiet Mode",                  option::no_arguments);

  util::cli::parse(argv, {&help_opt, &geo_opt, &root_opt, &map_opt, &script_opt, &quiet_opt});

  util::error::exit_when(argc == 1 || !(script_opt.count || geo_opt.count || root_opt.count),
    "[FATAL ERROR] Insufficient Arguments: ",
    "Must include arguments for geometry and ROOT directory or for a tracking script. \n",
    "              Run \'./tracker --help\' for more details.\n");

  if (script_opt.count) exit_on_missing_path(script_opt.argument, "Tracking Script");

  reader::script::tracking_options options;

  if (script_opt.count) {
    options = reader::script::read(script_opt.argument);
    geo_opt.count += !options.geometry_file.empty();
    map_opt.count += !options.geometry_map_file.empty();
    root_opt.count += !options.root_directory.empty();
  } else {
    options.geometry_file = geo_opt.count ? geo_opt.argument : "";
    options.geometry_map_file = map_opt.count ? map_opt.argument : "";
    options.root_directory = root_opt.count ? root_opt.argument : "";
  }

  if (geo_opt.count) exit_on_missing_path(options.geometry_file, "Geometry File");
  if (map_opt.count) exit_on_missing_path(options.geometry_map_file, "Geometry Map");
  if (root_opt.count) exit_on_missing_path(options.root_directory, "ROOT Directory");

  units::define();
  plot::init();
  geometry::open(options.geometry_file);
  const auto& detector_map = reader::root::import_detector_map(options.geometry_map_file);

  const auto paths = reader::root::search_directory(options.root_directory);
  std::cout << "File Count: " << paths.size() << "\n";
  for (const auto& path : paths) {
    std::cout << path << "\n";

    const auto events = map_opt.count
      ? reader::root::import_events(path,
          options.root_time_key, options.root_detector_key, detector_map)
      : reader::root::import_events(path,
          options.root_time_key, options.root_x_key, options.root_y_key, options.root_z_key);

    std::cout << "Processing " << events.size() << " Events\n";

    for (const auto& unsorted_event : events) {
      plot::canvas canvas;

      // demo for Prototype only
      for (const auto& name : geometry::full_structure()) {
        if (name == "world" || name == "Sandstone" || name == "Marl" || name == "Mix" || name == "Earth") continue;
        const auto limits = geometry::limits_of(name);
        canvas.add_point(limits.center, 0.25, plot::color::MAGENTA);
      }

      const auto event           = analysis::time_normalize(unsorted_event);
      const auto collapsed_event = analysis::collapse(event, options.collapse_size);
      const auto layered_event   = analysis::partition(collapsed_event, options.layer_depth, options.layer_axis).parts;

      for (const auto& point : event)
        canvas.add_point(point, 1.5, {90, 90, 90});

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

      auto tracks = analysis::fit_seeds(joined);

      std::cout << "\nTRACK FITTING (" << tracks.size() << "):\n\n";
      for (const auto& track : tracks) {
        const auto& event = track.event();
        for (const auto& point : event) {
          canvas.add_point(point, 0.5);
          const auto limits = geometry::limits_of_volume(point);
          canvas.add_point(limits.center, 0.25, plot::color::BLUE);
          canvas.add_box(limits.min, limits.max, 2, plot::color::BLUE);
        }
        canvas.add_line(event.front(), event.back());
        std::cout << track(event.front().z) << " " << track(event.back().z) << " " << track.chi_squared() << "\n";
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

#include "analysis.hh"
#include "geometry.hh"
#include "reader.hh"
#include "units.hh"

#include "util/command_line_parser.hh"
#include "util/error.hh"
#include "util/io.hh"
#include "util/string.hh"

//__Missing File Exit Command___________________________________________________________________
void exit_on_missing_file(const MATHUSLA::util::cli::option& option,
                          const std::string& name) {
  using namespace MATHUSLA::util;
  error::exit_when(option.count && !io::path_exists(option.argument),
    "[FATAL ERROR] ", name, " Missing: ",
    "The file ", option.argument, " cannot be found.\n");
}
//----------------------------------------------------------------------------------------------

//__Find Key for Option in Tracking Option Map__________________________________________________
void set_option_from_script_map(const MATHUSLA::TRACKER::reader::script::tracking_options& script_map,
                                const std::string& key,
                                MATHUSLA::util::cli::option& option) {
  using namespace MATHUSLA::TRACKER::reader::script;
  const auto& search = script_map.find(key);
  if (search != script_map.cend()) {
    option.argument = search->second.c_str();
    option.count = 1;
  }
}
//----------------------------------------------------------------------------------------------

//__Main Function: Tracker______________________________________________________________________
int main(int argc, char* argv[]) {
  using namespace MATHUSLA;
  using namespace MATHUSLA::TRACKER;
  using util::cli::option;

  option help_opt   ('h', "help",     "MATHUSLA Tracking Algorithm", option::no_arguments);
  option geo_opt    ('g', "geometry", "Geometry Import",             option::required_arguments);
  option root_opt   ('d', "dir",      "ROOT Data Directory",         option::required_arguments);
  option map_opt    ('m', "map",      "Detector Map",                option::required_arguments);
  option script_opt ('s', "script",   "Tracking Script",             option::required_arguments);
  option quiet_opt  ('q', "quiet",    "Quiet Mode",                  option::no_arguments);

  util::cli::parse(argv, {&help_opt, &geo_opt, &root_opt, &map_opt, &script_opt, &quiet_opt});

  util::error::exit_when(argc == 1 || (!script_opt.count && !geo_opt.count && !root_opt.count),
    "[FATAL ERROR] Insufficient Arguments: ",
    "Must include arguments for geometry and ROOT directory or tracking script. \n",
    "              Run \'./tracker --help\' for more details.\n");

  exit_on_missing_file(script_opt, "Tracking Script");

  type::r4_point collapse_size;
  type::real layer_depth;
  type::real line_width;
  type::integer seed_size;

  if (script_opt.count) {
    const auto& script_map = reader::script::read(script_opt.argument);
    const auto& script_map_end = script_map.cend();

    for (const auto& entry : script_map) {
      std::cout << entry.first << ": " << entry.second << "\n";
    }

    set_option_from_script_map(script_map, "geometry-file", geo_opt);
    set_option_from_script_map(script_map, "root-data", root_opt);
    set_option_from_script_map(script_map, "geometry-map", map_opt);

    const auto& collapse_size_search = script_map.find("collapse-size");
    if (collapse_size_search != script_map_end) {
      std::vector<std::string> point;
      util::string::split(collapse_size_search->second, point, ",");
      point.insert(point.cend(), {"1", "1", "1", "1"});
      collapse_size = {
        static_cast<type::real>(std::stold(point[0])),
        static_cast<type::real>(std::stold(point[1])),
        static_cast<type::real>(std::stold(point[2])),
        static_cast<type::real>(std::stold(point[3])) };
    }

  }

  exit_on_missing_file(geo_opt, "Geometry File");
  exit_on_missing_file(root_opt, "ROOT Directory");
  exit_on_missing_file(map_opt, "Geometry Map");

  units::define();
  /*
  std::cout << "DEMO:\n";
  geometry::open(geo_opt.argument);

  auto paths = reader::root::search_directory(root_opt.argument);
  std::cout << "File Count: " << paths.size() << "\n";

  for (const auto& path : paths) {
    std::cout << path << "\n";
    auto events = reader::root::import_events(path, {{"Time", "X", "Y", "Z"}});
    std::cout << "Processing " << events.size() << " Events\n";

    for (const auto& unsorted_event : events) {
      const auto& event = t_copy_sort(unsorted_event);
      const auto& collapsed_event = analysis::collapse(event, {20, 20, 20, 20});
      const auto& layered_event = analysis::partition(collapsed_event, 500).parts;

      for (const auto& point : event) {
        std::cout << "OLD " << point  << "\n";
      }
      std::cout << "\n";
      for (const auto& point : collapsed_event) {
        std::cout << "NEW " << point  << "\n";
      }
      std::cout << "\n";

      for (const auto& layer : layered_event) {
        for (size_t j = 0; j < layer.size(); ++j) {
          const auto& point = layer[j];
          std::cout << "LAYER " << point  << "\n";
        }
        std::cout << "\n\n";
      }

      std::cout << "\n";

      auto seeds = analysis::seed(3, unsorted_event, {20, 20, 20, 20}, 500, 0.8);

      std::cout << "\n";

      analysis::merge(seeds);  // not implemented!!!

    }
  }

  std::cout << "Done!\n";
  geometry::close();
  */
  return 0;
}
//----------------------------------------------------------------------------------------------

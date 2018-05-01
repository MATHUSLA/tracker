#include "analysis.hh"
#include "geometry.hh"
#include "root_helper.hh"
#include "units.hh"
#include "util.hh"

#include "util/CommandLineParser.hh"

using Option = MATHUSLA::CommandLineOption;

auto help_opt   = new Option('h', "help",     "MATHUSLA Tracking Algorithm", Option::NoArguments);
auto geo_opt    = new Option('g', "geometry", "Geometry Import",             Option::RequiredArguments);
auto root_opt   = new Option('d', "dir",      "ROOT Directory",              Option::RequiredArguments);
auto script_opt = new Option('s', "script",   "Tracking Script",             Option::RequiredArguments);
auto quiet_opt  = new Option('q', "quiet",    "Quiet Mode",                  Option::NoArguments);

int main(int argc, char* argv[]) {
  using namespace MATHUSLA;

  CommandLineParser::parse(argv, {help_opt, geo_opt, root_opt, script_opt, quiet_opt});

  error::exit_when(argc == 1 || !geo_opt->count || !root_opt->count,
    "[FATAL ERROR] Insufficient Arguments: ",
    "Must include arguments for geometry and ROOT directory. \n",
    "              Run \'./tracker --help\' for more details.\n");

  error::exit_when(!io::path_exists(geo_opt->argument),
    "[FATAL ERROR] Geometry File Missing: ",
    "The file ", geo_opt->argument, " cannot be found.\n");

  error::exit_when(!io::path_exists(root_opt->argument),
    "[FATAL ERROR] ROOT Directory Missing: ",
    "The directory ", root_opt->argument, " cannot be found.\n");

  Units::Define();

  using namespace MATHUSLA::TRACKER;

  std::cout << "DEMO:\n";
  geometry::open(geo_opt->argument);

  auto paths = root::search_directory(root_opt->argument);
  std::cout << "File Count: " << paths.size() << "\n";

  for (const auto& path : paths) {
    std::cout << path << "\n";
    auto events = root::import_events(path, {{"Time", "X", "Y", "Z"}});
    std::cout << "Processing " << events.size() << " Events\n";

    for (const auto& unsorted_event : events) {
      const auto& event = t_copy_sort(unsorted_event);
      const auto& collapsed_event = analysis::collapse(event, {2, 2, 2, 2});
      const auto& layered_event = analysis::partition(collapsed_event, 50).parts;

      for (const auto& point : event) {
        std::cout << "OLD " << point << " " << geometry::volume(point) << "\n";
      }
      std::cout << "\n";
      for (const auto& point : collapsed_event) {
        std::cout << "NEW " << point << " " << geometry::volume(point) << "\n";
      }
      std::cout << "\n";

      for (const auto& layer : layered_event) {
        for (size_t j = 0; j < layer.size(); ++j) {
          const auto& layer_point = layer[j];
          std::cout << "LAYER " << layer_point << " " << geometry::volume(layer_point) << "\n";
        }
        std::cout << "\n\n";
      }

      std::cout << "\n";
    }
  }

  std::cout << "Done!\n";

  geometry::close();

  return 0;
}

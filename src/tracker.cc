#include "analysis.hh"
#include "geometry.hh"
#include "util/CommandLineParser.hh"
#include "util.hh"

using Option = MATHUSLA::CommandLineOption;

auto help_opt   = new Option('h', "help",     "MATHUSLA Tracking Algorithm", Option::NoArguments);
auto geo_opt    = new Option('g', "geometry", "Geometry Import",             Option::RequiredArguments);
auto root_opt   = new Option('d', "dir",      "ROOT Directory",              Option::RequiredArguments);
auto script_opt = new Option('s', "script",   "Tracking Script",             Option::RequiredArguments);
auto quiet_opt  = new Option('q', "quiet",    "Quiet Mode",                  Option::NoArguments);

int main(int argc, char* argv[]) {
  MATHUSLA::CommandLineParser::parse(argv, {
    help_opt, geo_opt, root_opt, script_opt, quiet_opt});

  MATHUSLA::Error::exit_when(argc == 1 || !geo_opt->count || !root_opt->count,
    "[FATAL ERROR] Insufficient Arguments: ",
    "Must include arguments for geometry and ROOT directory. \n",
    "              Run \'./tracker --help\' for more details.\n");

  MATHUSLA::Error::exit_when(!MATHUSLA::IO::path_exists(geo_opt->argument),
    "[FATAL ERROR] Geometry File Missing: ",
    "The file ", geo_opt->argument, " cannot be found.\n");

  MATHUSLA::Error::exit_when(!MATHUSLA::IO::path_exists(root_opt->argument),
    "[FATAL ERROR] ROOT Directory Missing: ",
    "The directory ", root_opt->argument, " cannot be found.\n");

  MATHUSLA::Units::Define();

  using namespace MATHUSLA;
  using namespace MATHUSLA::TRACKER;

  Geometry::Initialize(geo_opt->argument);

  std::cout << "DEMO:\n";

  auto paths = Analysis::search_directory(root_opt->argument);
  std::cout << "File Count: " << paths.size() << "\n";

  for (const auto& path : paths) {

    std::cout << path << "\n";
    auto events = Analysis::import_events(path, {{"Time", "X", "Y", "Z"}});
    std::cout << "Processing " << events.size() << " Events\n";

    for (const auto& unsorted_event : events) {

      const auto& event = t_copy_sort(unsorted_event);
      const auto& collapsed_event = Analysis::collapse(event, {1, 1, 1, 1});
      const auto& layered_event = Analysis::partition<>(collapsed_event, 60).parts;

      for (size_t i = 0; i < event.size(); ++i) {

        const auto& old_point = event[i];
        std::cout << "OLD (" << old_point.t << " "
                             << old_point.x << " "
                             << old_point.y << " "
                             << old_point.z << ")    ";

        if (i < collapsed_event.size()) {
          const auto& new_point = collapsed_event[i];
          std::cout << "NEW (" << new_point.t << " "
                               << new_point.x << " "
                               << new_point.y << " "
                               << new_point.z << ")    ";
        }

        if (i < layered_event.size()) {
          for (size_t j = 0; j < layered_event[i].size(); ++j) {
            const auto& layer_point = layered_event[i][j];
            std::cout << "LAYER (" << layer_point.t << " "
                                   << layer_point.x << " "
                                   << layer_point.y << " "
                                   << layer_point.z << ")  ";
          }
        }

        std::cout << "\n";
      }
      std::cout << "\n";
    }
  }

  std::cout << "Done!\n";

  return 0;
}

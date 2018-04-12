#include "Geant4/G4RunManager.hh"
#include "Geant4/G4UImanager.hh"
#include "Geant4/G4UIExecutive.hh"
#include "Geant4/G4VisExecutive.hh"
#include "Geant4/QBBC.hh"

#include "geometry.hh"
#include "util/CommandLineParser.hh"
#include "util/FileIO.hh"

using Option = MATHUSLA::CommandLineOption;

auto help_opt  = new Option('h', "help",     "MATHUSLA Particle Tracker", Option::NoArguments);
auto geo_opt   = new Option('g', "geometry", "Geometry Import",           Option::RequiredArguments);
auto root_opt  = new Option('d', "dir",      "ROOT Dir",                  Option::RequiredArguments);
auto vis_opt   = new Option('v', "vis",      "Visualization",             Option::NoArguments);
auto quiet_opt = new Option('q', "quiet",    "Quiet Mode",                Option::NoArguments);

int main(int argc, char* argv[]) {
  MATHUSLA::CommandLineParser::parse(argv, {
    help_opt, geo_opt, root_opt, vis_opt, quiet_opt});

  if (argc == 1 || !geo_opt->count || !root_opt->count) {
    std::cout << "[FATAL ERROR] Insufficient Arguments: "
              << "Must include arguments for geometry and ROOT directory. \n"
              << "Run \'./tracker --help\' for more details.\n";
    exit(EXIT_FAILURE);
  }

  auto runManager = new G4RunManager;
  auto uiManager = G4UImanager::GetUIpointer();

  G4UIExecutive* ui = nullptr;
  G4VisExecutive* vis = nullptr;
  if (vis_opt->count) {
    ui = new G4UIExecutive(argc, argv);
    vis = new G4VisExecutive("Quiet");
    vis->Initialize();
  }

  runManager->SetUserInitialization(new QBBC);

  if (!MATHUSLA::IO::path_exists(geo_opt->argument)) {
    std::cout << "[FATAL ERROR] Geometry File Missing: "
              << "The file " << geo_opt->argument << " cannot be found.\n";
    exit(EXIT_FAILURE);
  }

  runManager->SetUserInitialization(
    new MATHUSLA::TRACKER::Geometry(geo_opt->argument));

  if (!MATHUSLA::IO::path_exists(root_opt->argument)) {
    std::cout << "[FATAL ERROR] ROOT Directory Missing: "
              << "The file " << root_opt->argument << " cannot be found.\n";
    exit(EXIT_FAILURE);
  }

  uiManager->ApplyCommand("/run/initialize");
  uiManager->ApplyCommand("/run/verbose 0");

  // example
  std::cout << MATHUSLA::TRACKER::Geometry::DetectorName({0, 0, -42.8*cm}) << "\n";

  if (vis_opt->count) {
    uiManager->ApplyCommand("/vis/open OGL 700x700-0+0");
    uiManager->ApplyCommand("/vis/drawVolume");
    uiManager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi  90. 180.");
    uiManager->ApplyCommand("/vis/viewer/set/lightsThetaPhi    180.   0.");
    uiManager->ApplyCommand("/vis/viewer/zoom 1.4");
    uiManager->ApplyCommand("/vis/scene/add/axes 0 0 0 1 m");
    ui->SessionStart();
    delete ui;
    delete vis;
  }

  delete runManager;
  return 0;
}

#include "Geant4/G4RunManager.hh"
#include "Geant4/G4UImanager.hh"
#include "Geant4/G4UIExecutive.hh"
#include "Geant4/G4VisExecutive.hh"
#include "Geant4/QBBC.hh"
#include "Geant4/G4StepLimiterPhysics.hh"

#include "geometry/Construction.hh"
#include "util/CommandLineParser.hh"
#include "util/FileIO.hh"

using Option = MATHUSLA::CommandLineOption;

auto help_opt  = new Option('h', "help",     "MATHUSLA Particle Tracker", Option::NoArguments);
auto geo_opt   = new Option('g', "geometry", "Geometry Import",           Option::RequiredArguments);
auto root_opt  = new Option('d', "dir",      "ROOT Dir",                  Option::RequiredArguments);
auto quiet_opt = new Option('q', "quiet",    "Quiet Mode",                Option::NoArguments);

int main(int argc, char* argv[]) {
  MATHUSLA::CommandLineParser::parse(argv, {
    help_opt, geo_opt, root_opt, quiet_opt});

  if (argc == 1 || !geo_opt->count || !root_opt->count) {
    std::cout << "[FATAL ERROR] Insufficient Arguments: "
              << "Must include arguments for geometry and ROOT directory. \n"
              << "Run \'./tracker --help\' for more details.\n";
    exit(EXIT_FAILURE);
  }

  auto ui = new G4UIExecutive(argc, argv);

  auto runManager = new G4RunManager;
  auto uiManager = G4UImanager::GetUIpointer();

  auto visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  runManager->SetUserInitialization(new QBBC);

  if (!MATHUSLA::IO::path_exists(geo_opt->argument)) {
    std::cout << "[FATAL ERROR] Geometry File Missing: "
              << "The file " << geo_opt->argument << " cannot be found.\n";
    exit(EXIT_FAILURE);
  }

  runManager->SetUserInitialization(
    new MATHUSLA::TRACKER::Construction(geo_opt->argument));

  if (!MATHUSLA::IO::path_exists(root_opt->argument)) {
    std::cout << "[FATAL ERROR] ROOT Directory Missing: "
              << "The file " << root_opt->argument << " cannot be found.\n";
    exit(EXIT_FAILURE);
  }

  // temporary for debugging
  uiManager->ApplyCommand("/run/initialize");
  uiManager->ApplyCommand("/vis/open OGL 700x700-0+0");
  uiManager->ApplyCommand("/vis/drawVolume");
  uiManager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi  90. 180.");
  uiManager->ApplyCommand("/vis/viewer/set/lightsThetaPhi    180.   0.");
  uiManager->ApplyCommand("/vis/viewer/zoom 1.4");
  uiManager->ApplyCommand("/vis/scene/add/axes 0 0 0 1 m");

  if (ui) {
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
  return 0;
}

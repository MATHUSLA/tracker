#include "geometry/Reader.hh"

#include <string_view>

#include "Geant4/G4GDMLParser.hh"

namespace MATHUSLA { namespace TRACKER {

G4VPhysicalVolume* Reader::Load(const std::string& path) {
  static G4ThreadLocal G4GDMLParser _parser;
  _parser.Clear();
  _parser.Read(path, false);
  return _parser.GetWorldVolume();
}

} } /* namespace MATHUSLA::TRACKER */

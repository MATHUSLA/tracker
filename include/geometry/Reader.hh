#ifndef GEOMETRY_READER_HH
#define GEOMETRY_READER_HH
#pragma once

#include "Geant4/G4VPhysicalVolume.hh"

namespace MATHUSLA { namespace TRACKER {

namespace Reader {

G4VPhysicalVolume* Load(const std::string& path);

}

} } /* namespace MATHUSLA */

#endif /* GEOMETRY_READER_HH */

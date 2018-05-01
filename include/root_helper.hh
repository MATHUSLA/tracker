#ifndef TRACKER__ROOT_HELPER_HH
#define TRACKER__ROOT_HELPER_HH
#pragma once

#include <unordered_map>

#include "ROOT/TFile.h"
#include "ROOT/TKey.h"

#include "event.hh"

namespace MATHUSLA { namespace TRACKER {

namespace root { ///////////////////////////////////////////////////////////////////////////////

//__ROOT Directory Search_______________________________________________________________________
std::vector<std::string> search_directory(const std::string& path);
//----------------------------------------------------------------------------------------------

//__ROOT Import Types___________________________________________________________________________
using point_keys    = std::array<std::string, 4>;
using detector_keys = std::array<std::string, 2>;
using detector_map  = std::unordered_map<type::integer, std::string>;
//----------------------------------------------------------------------------------------------

//__ROOT Detector Map Import____________________________________________________________________
detector_map import_detector_map(const std::string& path);
//----------------------------------------------------------------------------------------------

//__ROOT Event Import___________________________________________________________________________
event_vector import_events(const std::string& path,
                           const point_keys& keys);
event_vector import_events(const std::string& path,
                           const detector_keys& keys,
                           const detector_map& map);
//----------------------------------------------------------------------------------------------

//__ROOT File Key Traversal_____________________________________________________________________
template<class BinaryFunction>
inline BinaryFunction traverse_file(const std::string& path, BinaryFunction f) {
  auto file = TFile::Open(path.c_str());
  TIter next(file->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = static_cast<TKey*>(next())))
    f(file, key);
  delete key;
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ROOT_HELPER_HH */

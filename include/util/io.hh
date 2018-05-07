#ifndef UTIL__IO_HH
#define UTIL__IO_HH
#pragma once

#include <algorithm>
#include <iostream>
#include <sys/stat.h>

#if defined(_WIN32)
#include <windows.h>
#endif

namespace MATHUSLA {

namespace util { ///////////////////////////////////////////////////////////////////////////////

namespace io { /////////////////////////////////////////////////////////////////////////////////

//__Print Range of Printable Elements___________________________________________________________
template<class Range>
inline std::ostream& print_range(const Range& range, const std::string& spacer=" ") {
  std::for_each(range.cbegin(), --range.cend(),
    [&](const auto& element) { std::cout << element << spacer; });
  if (range.size() >= 1)
    std::cout << range.back();
  return std::cout;
}
//----------------------------------------------------------------------------------------------

//__Create Directory____________________________________________________________________________
inline int create_directory(const std::string& dir) {
  #if defined(_WIN32)
    return _mkdir(dir.c_str());
  #else
    return mkdir(dir.c_str(), 0733);
  #endif
}
//----------------------------------------------------------------------------------------------

//__Check if Path Exists________________________________________________________________________
inline bool path_exists(const std::string& path) {
  struct stat info;
  return !stat(path.c_str(), &info);
}
//----------------------------------------------------------------------------------------------

} /* namespace io */ ///////////////////////////////////////////////////////////////////////////

} /* namespace util */ /////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__IO_HH */

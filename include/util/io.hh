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
std::ostream& print_range(const Range& range, const std::string& spacer=" ", std::ostream& os=std::cout) {
  const auto& begin = range.cbegin();
  const auto& end = range.cend() - 1;
  if (end - begin >= 0) {
    std::for_each(begin, end, [&](const auto& element) { os << element << spacer; });
    os << *end;
  }
  return os;
}
//----------------------------------------------------------------------------------------------

//__Print Newline Characters____________________________________________________________________
inline std::ostream& newline(const size_t count=1, std::ostream& os=std::cout) {
  return os << std::string(count, '\n');
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

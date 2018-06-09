/*
 * include/tracker/util/io.hh
 *
 * Copyright 2018 Brandon Gomes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef UTIL__IO_HH
#define UTIL__IO_HH
#pragma once

#include <algorithm>
#include <iostream>

namespace MATHUSLA {

namespace util { ///////////////////////////////////////////////////////////////////////////////

namespace io { /////////////////////////////////////////////////////////////////////////////////

//__Hidden Includes To Keep Global Namespace Clean______________________________________________
#include <sys/stat.h>
#if defined(_WIN32)
#include <windows.h>
#endif
//----------------------------------------------------------------------------------------------

//__Print Range of Printable Elements___________________________________________________________
template<class Range>
std::ostream& print_range(const Range& range,
                          const std::string& spacer=" ",
                          const std::string& prefix="",
                          std::ostream& os=std::cout) {
  const auto& begin = range.cbegin();
  const auto& end = range.cend() - 1;
  if (end - begin >= 0) {
    std::for_each(begin, end, [&](const auto& element) { os << prefix << element << spacer; });
    os << prefix << *end;
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

//__Bold, Italicise, and Underline Text_________________________________________________________
inline std::ostream& bold(std::ostream& os=std::cout) {
  return os << "\e[1m";
}
inline std::ostream& italics(std::ostream& os=std::cout) {
  return os << "\e[3m";
}
inline std::ostream& underline(std::ostream& os=std::cout) {
  return os << "\e[4m";
}
inline std::ostream& reset_font(std::ostream& os=std::cout) {
  return os << "\e[0m";
}
//----------------------------------------------------------------------------------------------

} /* namespace io */ ///////////////////////////////////////////////////////////////////////////

} /* namespace util */ /////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__IO_HH */

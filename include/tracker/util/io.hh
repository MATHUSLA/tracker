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
#include <iomanip>

namespace MATHUSLA {

namespace util { namespace io { ////////////////////////////////////////////////////////////////

//__Enclosed Includes To Keep Global Namespace Clean____________________________________________
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
  const auto begin = std::cbegin(range);
  const auto end = std::cend(range) - 1;
  if (end - begin >= 0) {
    std::for_each(begin, end, [&](const auto& element) { os << prefix << element << spacer; });
    os << prefix << *end;
  }
  return os;
}
//----------------------------------------------------------------------------------------------

//__Repeat String_______________________________________________________________________________
inline std::ostream& repeat(const std::size_t count,
                            const std::string& string,
                            std::ostream& os=std::cout) {
  for (std::size_t i{}; i < count; ++i)
    os << string;
  return os;
}
inline std::ostream& repeat(const std::size_t count,
                            const char character,
                            std::ostream& os=std::cout) {
  return os << std::string(count, character);
}
//----------------------------------------------------------------------------------------------

//__Print Newline Characters____________________________________________________________________
inline std::ostream& newline(const std::size_t count=1,
                             std::ostream& os=std::cout) {
  return repeat(count, '\n', os);
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

//__Escape Character____________________________________________________________________________
inline std::ostream& escape(std::ostream& os=std::cout) {
  return os << "\x1B";
}
//----------------------------------------------------------------------------------------------

//__Bold, Italicise, and Underline Text_________________________________________________________
inline std::ostream& bold(std::ostream& os=std::cout) {
  return escape(os) << "[1m";
}
inline std::ostream& italics(std::ostream& os=std::cout) {
  return escape(os) << "[3m";
}
inline std::ostream& underline(std::ostream& os=std::cout) {
  return escape(os) << "[4m";
}
inline std::ostream& reset_font(std::ostream& os=std::cout) {
  return escape(os) << "[0m";
}
//----------------------------------------------------------------------------------------------

//__Set Failbit of Stream_______________________________________________________________________
inline std::ostream& set_failbit(std::ostream& os=std::cout) {
  os.setstate(std::ios_base::failbit);
  return os;
}
//----------------------------------------------------------------------------------------------

//__Clear Stream________________________________________________________________________________
inline std::ostream& clear(std::ostream& os=std::cout) {
  os.clear();
  return os;
}
//----------------------------------------------------------------------------------------------

//__Swap Buffer of Stream_______________________________________________________________________
inline std::ostream& swap_buffer(std::ostream& os,
                                 std::streambuf* buffer,
                                 std::streambuf* previous) {
  previous = os.rdbuf(buffer);
  return os;
}
//----------------------------------------------------------------------------------------------

//__Remove Buffers from Stream__________________________________________________________________
inline void remove_buffer(std::ostream& os) {
  swap_buffer(os, nullptr, nullptr);
}
template<class ...Stream>
void remove_buffer(std::ostream& os,
                   Stream& ...oss) {
  remove_buffer(os);
  remove_buffer(oss...);
}
//----------------------------------------------------------------------------------------------

} } /* namespace util::io */ ///////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__IO_HH */

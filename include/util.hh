#ifndef TRACKER__UTIL_HH
#define TRACKER__UTIL_HH
#pragma once

#include <cstdio>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

#if defined(_WIN32)
#include <windows.h>
#endif

namespace MATHUSLA {

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Forward Arguments to std::cout______________________________________________________________
template<class T> void _cout_forward(T&& arg) {
  std::cout << std::forward<T>(arg);
}
template<class T, class... Args> void _cout_forward(T&& arg, Args&&... args) {
  _cout_forward(arg);
  _cout_forward(args...);
}
//----------------------------------------------------------------------------------------------

} /* annonymous namespace */ ///////////////////////////////////////////////////////////////////

namespace Error { //////////////////////////////////////////////////////////////////////////////

//__Boolean Exit Convenience Function___________________________________________________________
template<class... Args>
void exit_when(bool value, int code, Args&&... msgs) {
  if (value) {
    _cout_forward(msgs...);
    exit(code);
  }
}
template<class... Args>
void exit_when(bool value, Args&&... msgs) {
  exit_when(value, EXIT_FAILURE, msgs...);
}
//----------------------------------------------------------------------------------------------

} /* namespace Error */ ////////////////////////////////////////////////////////////////////////

namespace IO { /////////////////////////////////////////////////////////////////////////////////

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

//__Remove File_________________________________________________________________________________
inline bool remove_file(const std::string& path) {
  return !remove(path.c_str());
}
//----------------------------------------------------------------------------------------------

//__Rename File_________________________________________________________________________________
inline bool rename_file(const std::string& path, const std::string& new_path) {
  return !rename(path.c_str(), new_path.c_str());
}
//----------------------------------------------------------------------------------------------

} /* namespace IO */ ///////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* TRACKER__UTIL_HH */

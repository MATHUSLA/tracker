#ifndef UTIL__ERROR_HH
#define UTIL__ERROR_HH
#pragma once

#include <iostream>

namespace MATHUSLA {

namespace util { ///////////////////////////////////////////////////////////////////////////////

namespace error { //////////////////////////////////////////////////////////////////////////////

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

} /* namespace error */ ////////////////////////////////////////////////////////////////////////

} /* namespace util */ /////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__ERROR_HH */
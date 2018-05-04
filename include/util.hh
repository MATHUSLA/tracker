#ifndef TRACKER__UTIL_HH
#define TRACKER__UTIL_HH
#pragma once

#include <algorithm>
#include <cstdio>
#include <functional>
#include <iostream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#if defined(_WIN32)
#include <windows.h>
#endif

namespace MATHUSLA {

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

namespace combinatorics { //////////////////////////////////////////////////////////////////////

namespace detail { /////////////////////////////////////////////////////////////////////////////
//__Boolean-Like Object_________________________________________________________________________
struct hidden_bool {
  bool data;
  hidden_bool() : data(0) {}
  hidden_bool(bool data) : data(data) {}

  hidden_bool(const hidden_bool& rhs) = default;
  hidden_bool(hidden_bool&& rhs)      = default;
  hidden_bool& operator=(const hidden_bool& rhs) = default;
  hidden_bool& operator=(hidden_bool&& rhs)      = default;

  operator bool() const { return data; }
};
//----------------------------------------------------------------------------------------------
} /* namespace detail */ ///////////////////////////////////////////////////////////////////////

//__Dynamic Bit Vector__________________________________________________________________________
class bit_vector : public std::vector<detail::hidden_bool> {
public:
  bit_vector(std::size_t count, std::size_t size) {
    resize(size, 0);
    for (auto i = count > size ? 0 : size - count; i < size; ++i) operator[](i) = true;
  }

  std::size_t count() const { return std::count(begin(), end(), true); }

  bool next_permutation() { return std::next_permutation(begin(), end()); }
};
//----------------------------------------------------------------------------------------------

//__Bit Vector Printer__________________________________________________________________________
inline std::ostream& operator<<(std::ostream& os, const bit_vector& bits) {
  for (const auto& bit : bits) os << bit;
  return os;
}
//----------------------------------------------------------------------------------------------

//__Vector of Bit Vectors_______________________________________________________________________
using bit_vector_sequence = std::vector<bit_vector>;
//----------------------------------------------------------------------------------------------

//__Generate Bit Vector Sequences In-place______________________________________________________
inline bit_vector_sequence generate_bit_sequences(std::vector<std::pair<std::size_t, std::size_t>> setup) {
  const auto&& size = setup.size();
  if (size == 0) return {};
  bit_vector_sequence out;
  out.reserve(size);
  for (const auto& pair : setup)
    out.push_back(bit_vector(pair.first, pair.second));
  return out;
}
//----------------------------------------------------------------------------------------------

//__First Order Bit Permutation Sequencer_______________________________________________________
template<class UnaryFunction>
inline UnaryFunction order1_permutations(std::size_t count, std::size_t total, UnaryFunction f) {
  bit_vector chooser(count, total);
  do { f(chooser); } while (chooser.next_permutation());
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Second Order Bit Permutation Sequencer______________________________________________________
template<class UnaryFunction>
inline UnaryFunction order2_permutations(std::size_t count, bit_vector_sequence& vectors, UnaryFunction f) {
  const auto&& size = vectors.size();
  bit_vector chooser(count, size);
  std::size_t index = 0;
  do {
    do {
      f(chooser);
      for (index = 0; index < size; ++index)
        if (chooser[index] && vectors[index].next_permutation()) break;
    } while (index != size);
  } while (chooser.next_permutation());
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

} /* namespace combinatorics */ ////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* TRACKER__UTIL_HH */

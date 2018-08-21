/*
 * include/tracker/util/bit_vector.hh
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

#ifndef UTIL__BIT_VECTOR_HH
#define UTIL__BIT_VECTOR_HH
#pragma once

#include <algorithm>
#include <ostream>
#include <vector>

namespace MATHUSLA {

namespace util { ///////////////////////////////////////////////////////////////////////////////

namespace detail { /////////////////////////////////////////////////////////////////////////////
//__Boolean-Like Object_________________________________________________________________________
struct hidden_bool {
  bool data;
  hidden_bool() : data(false) {}
  hidden_bool(const bool data) : data(data) {}
  hidden_bool(const hidden_bool& rhs)            = default;
  hidden_bool(hidden_bool&& rhs)                 = default;
  hidden_bool& operator=(const hidden_bool& rhs) = default;
  hidden_bool& operator=(hidden_bool&& rhs)      = default;
  operator bool() const { return data; }
};
//----------------------------------------------------------------------------------------------
} /* namespace detail */ ///////////////////////////////////////////////////////////////////////

//__Boolean-Like Bit____________________________________________________________________________
using bit = detail::hidden_bool;
//----------------------------------------------------------------------------------------------

//__Dynamic Bit Vector__________________________________________________________________________
class bit_vector : public std::vector<bit> {
public:
  explicit bit_vector(const std::size_t count,
                      const std::size_t size) {
    resize(size, false);
    for (std::size_t i = count > size ? 0UL : size - count; i < size; ++i)
      (*this)[i] = true;
  }

  explicit bit_vector(const std::size_t size) : bit_vector(0UL, size) {}

  std::size_t count() const noexcept {
    return std::count(cbegin(), cend(), true);
  }

  std::size_t first_set(const std::size_t start=0UL) const noexcept {
    const auto s = size();
    for (std::size_t i = start; i < s; ++i)
      if ((*this)[i]) return i;
    return s;
  }

  std::size_t last_set(const std::size_t start=0UL) const noexcept {
    const auto s = size();
    if (start >= s)
      return s;
    for (std::size_t i = s - start - 1UL; i != static_cast<std::size_t>(-1); --i)
      if ((*this)[i]) return i;
    return s;
  }

  std::size_t first_unset(const std::size_t start=0UL) const noexcept {
    const auto s = size();
    for (std::size_t i = start; i < s; ++i)
      if (!(*this)[i]) return i;
    return s;
  }

  std::size_t last_unset(const std::size_t start=0UL) const noexcept {
    const auto s = size();
    if (start >= s)
      return s;
    for (std::size_t i = s - start - 1UL; i != static_cast<std::size_t>(-1); --i)
      if (!(*this)[i]) return i;
    return s;
  }

  bit_vector& set() noexcept {
    const auto s = size();
    for (std::size_t i{}; i < s; ++i)
      (*this)[i] = true;
    return *this;
  }

  bit_vector& set(const std::size_t index) {
    (*this)[index] = true;
    return *this;
  }

  bit_vector& reset() noexcept {
    const auto s = size();
    for (std::size_t i{}; i < s; ++i)
      (*this)[i] = false;
    return *this;
  }

  bit_vector& reset(const std::size_t index) {
    (*this)[index] = false;
    return *this;
  }

  std::size_t extend(const std::size_t count=1) {
    for (std::size_t i{}; i < count; ++i)
      push_back(false);
    return size();
  }

  bool next_permutation() noexcept {
    return std::next_permutation(begin(), end());
  }

  template<class Range, class BinaryFunction>
  BinaryFunction traverse(const Range& range,
                          BinaryFunction f) const noexcept {
    auto begin = std::cbegin(range);
    const auto s = std::min(size(), static_cast<std::size_t>(std::cend(range) - begin));
    for (std::size_t i{}; i < s; ++i)
      f((*this)[i], *begin++);
    return std::move(f);
  }

  template<class InputRange, class OutputRange, class BinaryFunction>
  BinaryFunction conditional_push_back(const InputRange& in,
                                       OutputRange& out,
                                       BinaryFunction f) const noexcept {
    auto back_inserter = std::back_inserter(out);
    traverse(in, [&](const auto bit, const auto& item) { if (f(bit, item)) *back_inserter++ = item; });
    return std::move(f);
  }

  template<class InputRange, class OutputRange>
  void set_conditional_push_back(const InputRange& in,
                                 OutputRange& out) const noexcept {
    conditional_push_back(in, out, [](const auto bit, const auto&) { return bit; });
  }

  template<class InputRange, class OutputRange>
  void unset_conditional_push_back(const InputRange& in,
                                   OutputRange& out) const noexcept {
    conditional_push_back(in, out, [](const auto bit, const auto&) { return !bit; });
  }
};
//----------------------------------------------------------------------------------------------

//__Bit Vector Printer__________________________________________________________________________
inline std::ostream& operator<<(std::ostream& os,
                                const bit_vector& bits) {
  for (const auto& b : bits) os << b;
  return os;
}
//----------------------------------------------------------------------------------------------

//__Vector of Bit Vectors_______________________________________________________________________
using bit_vector_sequence = std::vector<bit_vector>;
//----------------------------------------------------------------------------------------------

//__Generate Bit Vector Sequences In-place______________________________________________________
inline bit_vector_sequence generate_bit_sequences(const std::vector<std::pair<std::size_t, std::size_t>>& setup) {
  const auto size = setup.size();
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
UnaryFunction order1_permutations(const std::size_t count,
                                  const std::size_t total,
                                  UnaryFunction f) {
  bit_vector chooser(count, total);
  do { f(chooser); } while (chooser.next_permutation());
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Second Order Bit Permutation Sequencer______________________________________________________
template<class UnaryFunction>
UnaryFunction order2_permutations(const std::size_t count,
                                  bit_vector_sequence& vectors,
                                  UnaryFunction f) {
  const auto size = vectors.size();
  bit_vector chooser(count, size);
  std::size_t index;
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

} /* namespace util */ /////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__BIT_VECTOR_HH */

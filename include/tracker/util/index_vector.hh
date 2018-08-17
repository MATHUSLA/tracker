/*
 * include/tracker/util/index_vector.hh
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

#ifndef UTIL__INDEX_VECTOR_HH
#define UTIL__INDEX_VECTOR_HH
#pragma once

#include <algorithm>
#include <cstdint>
#include <ostream>
#include <vector>
#include <tuple>

#include <tracker/util/type.hh>

namespace MATHUSLA {

namespace util { ///////////////////////////////////////////////////////////////////////////////

namespace detail { /////////////////////////////////////////////////////////////////////////////

//__Constexpr Log2______________________________________________________________________________
constexpr std::size_t log2(const std::size_t n) {
  return n <   2 ? 0 :
        (n <   4 ? 1 :
        (n <   8 ? 2 :
        (n <  16 ? 3 :
        (n <  32 ? 4 :
        (n <  64 ? 5 :
        (n < 128 ? 6 :
        (n < 256 ? 7 :
        (n < 512 ? 8 : 9))))))));
}
//----------------------------------------------------------------------------------------------

//__Allowed Unsigned Integer Types______________________________________________________________
using uint_types = std::tuple<std::uint_least8_t,
                              std::uint_least8_t,
                              std::uint_least8_t,
                              std::uint_fast8_t,
                              std::uint_fast16_t,
                              std::uint_fast32_t,
                              std::uint_fast64_t>;
//----------------------------------------------------------------------------------------------

//__Unsigned Integer with N Bits________________________________________________________________
template<std::size_t N>
using uint_size = typename std::tuple_element<log2(N), uint_types>::type;
//----------------------------------------------------------------------------------------------

} /* namespace detail */ ///////////////////////////////////////////////////////////////////////

//__Default Index Vector Size Type______________________________________________________________
constexpr std::size_t default_index_vector_size = 32UL;
//----------------------------------------------------------------------------------------------

//__Index Vector Type___________________________________________________________________________
template<std::size_t N=default_index_vector_size>
class index_vector : public std::vector<detail::uint_size<N>> {
public:
  using value_type = detail::uint_size<N>;
  using value_pair_type = std::pair<value_type, value_type>;

  index_vector() = default;

  explicit index_vector(const std::size_t size) {
    this->reserve(size);
    for (std::size_t i = 0; i < size; ++i)
      this->push_back(i);
  }

  const value_type width() const noexcept {
    return this->back() - this->front();
  }

  const value_pair_type range() const noexcept {
    return std::make_pair(this->front(), this->back());
  }

  std::size_t extend(const std::size_t count=1) {
    for (std::size_t i = 0; i < count; ++i)
      this->push_back(this->size());
    return this->size();
  }

  template<class Range, typename Value = typename Range::value_type>
  void joint_push_back(Range& range,
                       const Value& value) {
    *std::back_inserter(range) = value;
    this->push_back(util::type::size(range) - 1);
  }

  template<class Range, typename Value = typename Range::value_type>
  void joint_push_back(Range& range,
                       Value&& value) {
    *std::back_inserter(range) = std::move(value);
    this->push_back(util::type::size(range) - 1);
  }

  template<class Range, class BinaryFunction>
  BinaryFunction traverse(const Range& range,
                          BinaryFunction f) const noexcept {
    auto begin = std::cbegin(range);
    const auto s = std::min(this->size(), util::type::distance(begin, std::cend(range)));
    for (std::size_t i = 0; i < s; ++i) {
      const auto index = (*this)[i];
      f(index, *(begin + index));
    }
    return std::move(f);
  }

  template<class InputRange, class OutputRange>
  void conditional_push_back(const InputRange& in,
                             OutputRange& out) const {
    auto out_inserter = std::back_inserter(out);
    traverse(in, [&](const auto, const auto& value) { *out_inserter++ = value; });
  }
};
//----------------------------------------------------------------------------------------------

//__Index Vector Stream Overload________________________________________________________________
template<std::size_t N=default_index_vector_size>
inline std::ostream& operator<<(std::ostream& os,
                                const index_vector<N>& indices) {
  const auto size = indices.size();
  if (size == 0UL)
    return os;
  for (std::size_t i = 0; i < size - 1; ++i)
    os << indices[i] << " ";
  return os << indices.back();
}
//----------------------------------------------------------------------------------------------

//__Vector of Index Vectors_____________________________________________________________________
template<std::size_t N=default_index_vector_size>
using index_vector_sequence = std::vector<index_vector<N>>;
//----------------------------------------------------------------------------------------------

} /* namespace util */ /////////////////////////////////////////////////////////////////////////

} /* namespace MATHUSLA */

#endif /* UTIL__INDEX_VECTOR_HH */

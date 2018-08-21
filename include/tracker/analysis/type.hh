/*
 * include/tracker/analysis/type.hh
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

#ifndef TRACKER__ANALYSIS__TYPE_HH
#define TRACKER__ANALYSIS__TYPE_HH
#pragma once

#include <tracker/core/type.hh>
#include <tracker/core/units.hh>

#include <tracker/util/functional.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

using namespace type;

//__Event Types_________________________________________________________________________________
using hit = r4_point;
using event = std::vector<hit>;
using event_vector = std::vector<event>;
//----------------------------------------------------------------------------------------------

//__Extended Event Types________________________________________________________________________
struct full_hit { real t, x, y, z; r4_point width; };
using full_event = std::vector<full_hit>;
using full_event_vector = std::vector<full_event>;
//----------------------------------------------------------------------------------------------

//__Fitting Parameter Type______________________________________________________________________
struct fit_parameter { real value, error, min, max; };
//----------------------------------------------------------------------------------------------

//__Full Hit Stream Operator Overload___________________________________________________________
inline std::ostream& operator<<(std::ostream& os,
                                const full_hit& point) {
  return os << "[(" << point.t / units::time   << ", "
                    << point.x / units::length << ", "
                    << point.y / units::length << ", "
                    << point.z / units::length
            << ") +/- " << units::scale_r4_length(point.width) << "]";
}
//----------------------------------------------------------------------------------------------

//__Full Hit Equality___________________________________________________________________________
constexpr bool operator==(const full_hit& left,
                          const full_hit& right) {
  return left.t == right.t && left.x == right.x && left.y == right.y && left.z == right.z
      && left.width == right.width;
}
//----------------------------------------------------------------------------------------------

//__Fitting Parameter Equality__________________________________________________________________
constexpr bool operator==(const fit_parameter& left,
                          const fit_parameter& right) {
  return left.value == right.value && left.error == right.error
    && left.min == right.min && left.max == right.max;
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

namespace std { ////////////////////////////////////////////////////////////////////////////////

//__Analysis Event Hash Function________________________________________________________________
template<>
struct hash<MATHUSLA::TRACKER::analysis::event> {
  std::size_t operator()(const MATHUSLA::TRACKER::analysis::event& in) const {
    return MATHUSLA::util::functional::hash_combine_range(in);
  }
};
//----------------------------------------------------------------------------------------------

//__Analysis Event Vector Hash Function_________________________________________________________
template<>
struct hash<MATHUSLA::TRACKER::analysis::event_vector> {
  std::size_t operator()(const MATHUSLA::TRACKER::analysis::event_vector& in) const {
    return MATHUSLA::util::functional::hash_combine_range(in);
  }
};
//----------------------------------------------------------------------------------------------

//__Analysis Full Hit Hash Function_____________________________________________________________
template<>
struct hash<MATHUSLA::TRACKER::analysis::full_hit> {
  std::size_t operator()(const MATHUSLA::TRACKER::analysis::full_hit& in) const {
    return MATHUSLA::util::functional::hash_combine(
      in.t, in.x, in.y, in.z, in.width.t, in.width.x, in.width.y, in.width.z);
  }
};
//----------------------------------------------------------------------------------------------

//__Analysis Full Event Hash Function___________________________________________________________
template<>
struct hash<MATHUSLA::TRACKER::analysis::full_event> {
  std::size_t operator()(const MATHUSLA::TRACKER::analysis::full_event& in) const {
    return MATHUSLA::util::functional::hash_combine_range(in);
  }
};
//----------------------------------------------------------------------------------------------

//__Analysis Full Event Vector Hash Function____________________________________________________
template<>
struct hash<MATHUSLA::TRACKER::analysis::full_event_vector> {
  std::size_t operator()(const MATHUSLA::TRACKER::analysis::full_event_vector& in) const {
    return MATHUSLA::util::functional::hash_combine_range(in);
  }
};
//----------------------------------------------------------------------------------------------

//__Analysis Full Event Vector Hash Function____________________________________________________
template<>
struct hash<MATHUSLA::TRACKER::analysis::fit_parameter> {
  std::size_t operator()(const MATHUSLA::TRACKER::analysis::fit_parameter& in) const {
    return MATHUSLA::util::functional::hash_combine(in.value, in.error, in.min, in.max);
  }
};
//----------------------------------------------------------------------------------------------

} /* namespace std */ //////////////////////////////////////////////////////////////////////////

#endif /* TRACKER__ANALYSIS__TYPE_HH */

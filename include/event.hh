#ifndef TRACKER__EVENT_HH
#define TRACKER__EVENT_HH
#pragma once

#include "point.hh"

namespace MATHUSLA { namespace TRACKER {

using namespace type;

//__Event Types_________________________________________________________________________________
using event_points = r4_point_vector;
using event_vector = std::vector<event_points>;
//----------------------------------------------------------------------------------------------

//__Event Tuple Types___________________________________________________________________________
template<std::size_t N> using event_tuple = r4_point_array<N>;
using event_pair      = event_tuple<2>;
using event_triple    = event_tuple<3>;
using event_quadruple = event_tuple<4>;
template<std::size_t N> using event_tuple_vector = std::vector<event_tuple<N>>;
//----------------------------------------------------------------------------------------------

//__Event Partition Type________________________________________________________________________
struct event_partition { event_vector parts; Coordinate coordinate; };
//----------------------------------------------------------------------------------------------

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__EVENT_HH */

/*
 * src/tracker/vertex.cc
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

#include <tracker/vertex.hh>

#include <ROOT/TMinuit.h>

#include <tracker/stat.hh>

#include <tracker/util/io.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Vertex Constructor__________________________________________________________________________
vertex::vertex(const track_vector& tracks) : _tracks(tracks) {

}
//----------------------------------------------------------------------------------------------

//__Error in Calculation of Vertex Point________________________________________________________
const r4_point vertex::point_error() const {
  return {};
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared Test Statistic__________________________________________________________________
real vertex::chi_squared() const {
  return std::accumulate(_delta_chi2.cbegin(), _delta_chi2.cend(), 0.0L);
}
//----------------------------------------------------------------------------------------------

//__Vertex Degrees of Freedom___________________________________________________________________
size_t vertex::degrees_of_freedom() const {
  return -4 + std::accumulate(_tracks.cbegin(), _tracks.cend(), 0UL,
    [](const auto sum, const auto& track) { return sum + track.degrees_of_freedom(); });
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared per Degree of Freedom___________________________________________________________
real vertex::chi_squared_per_dof() const {
  return chi_squared() / degrees_of_freedom();
}
//----------------------------------------------------------------------------------------------

//__Vertex Output Stream Operator_______________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const vertex& vertex) {
  static const std::string bar(80, '-');
  os << bar << "\n";

  os << "* Vertex:\n"
     << "    " << vertex.point() << " (+/- " << vertex.point_error() << ")\n";

  os << "* Tracks: \n";
  for (const auto& track : vertex.tracks()) {
    os << "    (" << track.t0_value() << ", "
                  << track.x0_value() << ", "
                  << track.y0_value() << ", "
                  << track.z0_value() << ", "
                  << track.vx_value() << ", "
                  << track.vy_value() << ", "
                  << track.vz_value() << ")\n";
  }

  os.precision(7);
  os << "* Statistics: \n"
     << "    dof:      " << vertex.degrees_of_freedom()               << "\n"
     << "    chi2:     " << vertex.chi_squared() << " = ";
  util::io::print_range(vertex.chi_squared_vector(), " + ", "", os)   << "\n";
  os << "    chi2/dof: " << vertex.chi_squared_per_dof()              << "\n"
     << "    p-value:  " << stat::chi_squared_p_value(vertex)         << "\n";

  return os << bar;
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

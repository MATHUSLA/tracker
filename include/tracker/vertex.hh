/*
 * include/tracker/vertex.hh
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

#ifndef TRACKER__VERTEX_HH
#define TRACKER__VERTEX_HH
#pragma once

#include <tracker/track.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

using namespace type;

//__Vertex Object_______________________________________________________________________________
class vertex {
public:
  enum class parameter { T, X, Y, Z };
  struct fit_parameters { fit_parameter t, x, y, z; };

  static constexpr std::size_t free_parameter_count = 4UL;
  using covariance_matrix_type = real_array<free_parameter_count * free_parameter_count>;

  vertex(const track_vector& tracks);

  vertex(const vertex& rhs) = default;
  vertex(vertex&& rhs)      = default;
  vertex& operator=(const vertex& rhs) = default;
  vertex& operator=(vertex&& rhs)      = default;

  const r4_point point() const;
  const r4_point point_error() const;

  const fit_parameters guess_fit() const { return _guess; }
  const fit_parameters final_fit() const { return _final; }

  const fit_parameter t() const { return _final.t; }
  const fit_parameter x() const { return _final.x; }
  const fit_parameter y() const { return _final.y; }
  const fit_parameter z() const { return _final.z; }
  const fit_parameter fit_of(const parameter p) const;

  real t_value() const { return _final.t.value; }
  real x_value() const { return _final.x.value; }
  real y_value() const { return _final.y.value; }
  real z_value() const { return _final.z.value; }
  real value(const parameter p) const;

  real t_error() const { return _final.t.error; }
  real x_error() const { return _final.x.error; }
  real y_error() const { return _final.y.error; }
  real z_error() const { return _final.z.error; }
  real error(const parameter p) const;

  real chi_squared() const;
  size_t degrees_of_freedom() const;
  real chi_squared_per_dof() const;
  const real_vector& chi_squared_vector() const { return _delta_chi2; }

  real variance(const parameter p) const;
  real covariance(const parameter p,
                  const parameter q) const;
  const covariance_matrix_type& covariance_matrix() const { return _covariance; }

  const track_vector tracks() const { return _tracks; };
  size_t count() const { return _tracks.size(); }

private:
  fit_parameters _guess, _final;
  track_vector _tracks;
  real_vector _delta_chi2;
  covariance_matrix_type _covariance;
};
//----------------------------------------------------------------------------------------------

//__Vertex Output Stream Operator_______________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const vertex& vertex);
//----------------------------------------------------------------------------------------------

//__Vector of Vertices__________________________________________________________________________
using vertex_vector = std::vector<vertex>;
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__VERTEX_HH */

/*
 * include/tracker/analysis/vertex.hh
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

#ifndef TRACKER__ANALYSIS__VERTEX_HH
#define TRACKER__ANALYSIS__VERTEX_HH
#pragma once

#include <tracker/analysis/track.hh>
#include <tracker/plot.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Vertex Object_______________________________________________________________________________
class vertex {
public:
  class tree;
  enum class parameter { T, X, Y, Z };
  struct fit_parameters { fit_parameter t, x, y, z; };

  static constexpr std::size_t free_parameter_count = 4UL;
  using covariance_matrix_type = real_array<free_parameter_count * free_parameter_count>;

  vertex(const track_vector& tracks);
  vertex(track_vector&& tracks);

  vertex(const vertex& rhs) = default;
  vertex(vertex&& rhs)      = default;
  vertex& operator=(const vertex& rhs) = default;
  vertex& operator=(vertex&& rhs)      = default;

  const r4_point point() const;
  const r4_point point_error() const;

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

  const fit_parameters guess_fit() const { return _guess; }
  const fit_parameters final_fit() const { return _final; }

  bool fit_diverged() const noexcept;
  bool fit_converged() const noexcept { return !fit_diverged(); }

  real_vector distances() const;
  real_vector distance_errors() const;

  real chi_squared() const;
  std::size_t degrees_of_freedom() const;
  real chi_squared_per_dof() const;
  real chi_squared_p_value() const;
  const real_vector& chi_squared_vector() const { return _delta_chi2; }

  real variance(const parameter p) const;
  real covariance(const parameter p,
                  const parameter q) const;
  const covariance_matrix_type& covariance_matrix() const { return _covariance; }

  const track_vector tracks() const { return _tracks; };
  std::size_t size() const { return _tracks.size(); }
  bool empty() const { return _tracks.size() <= 1; }

  track_vector::iterator begin() { return _tracks.begin(); }
  track_vector::const_iterator begin() const { return _tracks.cbegin(); }
  track_vector::const_iterator cbegin() const { return _tracks.cbegin(); }
  track_vector::iterator end() { return _tracks.end(); }
  track_vector::const_iterator end() const { return _tracks.cend(); }
  track_vector::const_iterator cend() const { return _tracks.cend(); }

  std::size_t reset(const track_vector& tracks);

  std::size_t insert(const track& track);
  std::size_t insert(const track_vector& tracks);

  std::size_t remove(const std::size_t index);
  std::size_t remove(const std::vector<std::size_t>& indices);

  std::size_t prune_on_chi_squared(const real max_chi_squared);

  struct plotting_keys {
    plot::histogram::name_type t, x, y, z,
      t_error, x_error, y_error, z_error,
      distance, distance_error,
      chi_squared_per_dof,
      size;
  };

  void fill_plots(plot::histogram_collection& collection,
                  const plotting_keys& keys) const;

  void draw(plot::canvas& canvas,
            const real size,
            const plot::color color,
            const bool with_errors=false) const;

  void draw_guess(plot::canvas& canvas,
                  const real size,
                  const plot::color color,
                  const bool with_errors=false) const;

  bool operator==(const vertex& other) const noexcept { return _tracks == other._tracks; }
  bool operator!=(const vertex& other) const noexcept { return !(*this == other);        }

protected:
  fit_parameters _guess, _final;
  track_vector _tracks;
  real_vector _delta_chi2;
  covariance_matrix_type _covariance;
};
//----------------------------------------------------------------------------------------------

//__Track Fitting Parameter Equality____________________________________________________________
constexpr bool operator==(const vertex::fit_parameters& left,
                          const vertex::fit_parameters& right) {
  return left.t == right.t
      && left.x == right.x
      && left.y == right.y
      && left.z == right.z;
}
constexpr bool operator!=(const vertex::fit_parameters& left,
                          const vertex::fit_parameters& right) {
  return !(left == right);
}
//----------------------------------------------------------------------------------------------

//__Track Data Tree Specialization______________________________________________________________
class vertex::tree : public analysis::tree {
public:
  using branch_value_type = std::vector<double>;
  using branch_type = branch<branch_value_type>;

  tree(const std::string& name);
  tree(const std::string& name,
       const std::string& title);

  branch_type t, x, y, z,
              t_error, x_error, y_error, z_error,
              chi_squared, chi_squared_per_dof, chi_squared_p_value,
              size;

  void insert(const vertex& vertex);
  void clear();
  void reserve(std::size_t capacity);

private:
  std::vector<std::reference_wrapper<branch_type>> _branches;
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

#endif /* TRACKER__ANALYSIS__VERTEX_HH */

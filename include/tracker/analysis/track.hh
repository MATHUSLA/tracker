/*
 * include/tracker/analysis/track.hh
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

#ifndef TRACKER__ANALYSIS__TRACK_HH
#define TRACKER__ANALYSIS__TRACK_HH
#pragma once

#include <tracker/analysis/type.hh>
#include <tracker/analysis/tree.hh>
#include <tracker/geometry.hh>
#include <tracker/plot.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Track Object________________________________________________________________________________
class track {
public:
  class tree;
  // TODO: class graph;
  enum class parameter { T0, X0, Y0, Z0, VX, VY, VZ };
  struct fit_parameters { fit_parameter t0, x0, y0, z0, vx, vy, vz; };

  static constexpr std::size_t free_parameter_count = 6UL;
  using covariance_matrix_type = real_array<free_parameter_count * free_parameter_count>;

  using container_type = analysis::full_event;
  using value_type = typename container_type::value_type;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;
  using reverse_iterator = typename container_type::reverse_iterator;
  using const_reverse_iterator = typename container_type::const_reverse_iterator;

  track();
  track(const event& points,
        const Coordinate direction=Coordinate::Z);
  track(const full_event& points,
        const Coordinate direction=Coordinate::Z);
  track(const track& rhs) = default;
  track(track&& rhs) noexcept = default;
  track& operator=(const track& rhs) = default;
  track& operator=(track&& rhs) noexcept = default;

  const r4_point operator()(const real p) const;

  const r4_point origin() const;
  const r3_point ray() const;
  const r4_point origin_error() const;
  const r3_point ray_error() const;

  const r4_point point(const real t) const;
  const r4_point point_error(const real t) const;

  const r4_point at_t(const real t) const;
  const r4_point at_x(const real x) const;
  const r4_point at_y(const real y) const;
  const r4_point at_z(const real z) const;
  const r4_point at(const Coordinate c,
                    const real r) const;

  const r4_point error_at_t(const real t) const;
  const r4_point error_at_x(const real x) const;
  const r4_point error_at_y(const real y) const;
  const r4_point error_at_z(const real z) const;
  const r4_point error_at(const Coordinate c,
                          const real r) const;

  const fit_parameter t0() const { return _final.t0; }
  const fit_parameter x0() const { return _final.x0; }
  const fit_parameter y0() const { return _final.y0; }
  const fit_parameter z0() const { return _final.z0; }
  const fit_parameter vx() const { return _final.vx; }
  const fit_parameter vy() const { return _final.vy; }
  const fit_parameter vz() const { return _final.vz; }
  const fit_parameter fit_of(const parameter p) const;

  real t0_value() const { return _final.t0.value; }
  real x0_value() const { return _final.x0.value; }
  real y0_value() const { return _final.y0.value; }
  real z0_value() const { return _final.z0.value; }
  real vx_value() const { return _final.vx.value; }
  real vy_value() const { return _final.vy.value; }
  real vz_value() const { return _final.vz.value; }
  real value(const parameter p) const;

  real t0_error() const { return _final.t0.error; }
  real x0_error() const { return _final.x0.error; }
  real y0_error() const { return _final.y0.error; }
  real z0_error() const { return _final.z0.error; }
  real vx_error() const { return _final.vx.error; }
  real vy_error() const { return _final.vy.error; }
  real vz_error() const { return _final.vz.error; }
  real error(const parameter p) const;

  const fit_parameters guess_fit() const { return _guess; }
  const fit_parameters final_fit() const { return _final; }

  bool fit_diverged() const noexcept;
  bool fit_converged() const noexcept { return !fit_diverged(); }

  real beta() const;
  real beta_error() const;
  const r3_point unit() const;
  const r3_point unit_error() const;

  real angle() const;
  real angle_error() const;

  real chi_squared() const;
  std::size_t degrees_of_freedom() const;
  real chi_squared_per_dof() const;
  real chi_squared_p_value() const;
  const real_vector& chi_squared_vector() const { return _delta_chi2; }

  real variance(const parameter p) const;
  real covariance(const parameter p,
                  const parameter q) const;
  const covariance_matrix_type& covariance_matrix() const { return _covariance; }

  const hit front() const;
  const hit back() const;
  const analysis::event event() const;

  const full_hit full_front() const noexcept { return _full_event.front(); }
  const full_hit full_back() const noexcept { return _full_event.back(); }
  const analysis::full_event& full_event() const noexcept { return _full_event; }

  template<class Geometry=void>
  const geometry::structure_vector detectors() const {
    geometry::structure_vector out;
    out.reserve(size());
    util::algorithm::back_insert_transform(_full_event, out,
      [](const auto& point) { return geometry::custom::volume<Geometry>(reduce_to_r3(point)); });
    return out;
  }

  Coordinate direction() const noexcept { return _direction; }

  std::size_t size() const noexcept { return _full_event.size(); }
  bool empty() const noexcept { return _full_event.empty(); }

  iterator       begin()        noexcept { return _full_event.begin();  }
  const_iterator begin()  const noexcept { return _full_event.cbegin(); }
  iterator       end()          noexcept { return _full_event.end();    }
  const_iterator end()    const noexcept { return _full_event.cend();   }
  const_iterator cbegin() const noexcept { return _full_event.cbegin(); }
  const_iterator cend()   const noexcept { return _full_event.cend();   }

  reverse_iterator       rbegin()        noexcept { return _full_event.rbegin();  }
  const_reverse_iterator rbegin()  const noexcept { return _full_event.crbegin(); }
  reverse_iterator       rend()          noexcept { return _full_event.rend();    }
  const_reverse_iterator rend()    const noexcept { return _full_event.crend();   }
  const_reverse_iterator crbegin() const noexcept { return _full_event.crbegin(); }
  const_reverse_iterator crend()   const noexcept { return _full_event.crend();   }

  std::size_t reset(const analysis::event& points);
  std::size_t reset(const analysis::full_event& points);

  std::size_t reset(const analysis::event& points,
                    const Coordinate direction);
  std::size_t reset(const analysis::full_event& points,
                    const Coordinate direction);

  std::size_t insert(const hit& point);
  std::size_t insert(const analysis::event& points);
  std::size_t insert(const full_hit& point);
  std::size_t insert(const analysis::full_event& points);

  std::size_t remove(const std::size_t index);
  std::size_t remove(const std::vector<std::size_t>& indices);

  std::size_t prune_on_chi_squared(const real max_chi_squared);

  void clear();

  void reparameterize(const Coordinate direction);

  struct plotting_keys {
    plot::histogram::name_type t0, x0, y0, z0, vx, vy, vz,
      t0_error, x0_error, y0_error, z0_error, vx_error, vy_error, vz_error,
      chi_squared, chi_squared_per_dof, chi_squared_p_value,
      size,
      beta, beta_error,
      angle, angle_error;
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

  bool operator==(const track& other) const noexcept {
    return _direction == other._direction && _full_event == other._full_event;
  }
  bool operator!=(const track& other) const noexcept {
    return !(*this == other);
  }

  std::size_t hash() const;

protected:
  fit_parameters _guess, _final;
  analysis::full_event _full_event;
  real_vector _delta_chi2;
  covariance_matrix_type _covariance;
  Coordinate _direction;
};
//----------------------------------------------------------------------------------------------

//__Track Fitting Parameter Equality____________________________________________________________
constexpr bool operator==(const track::fit_parameters& left,
                          const track::fit_parameters& right) {
  return left.t0 == right.t0
      && left.x0 == right.x0
      && left.y0 == right.y0
      && left.z0 == right.z0
      && left.vx == right.vx
      && left.vy == right.vy
      && left.vz == right.vz;
}
constexpr bool operator!=(const track::fit_parameters& left,
                          const track::fit_parameters& right) {
  return !(left == right);
}
//----------------------------------------------------------------------------------------------

//__Track Output Stream Operator________________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const track& track);
//----------------------------------------------------------------------------------------------

//__Vector of Tracks____________________________________________________________________________
using track_vector = std::vector<track>;
//----------------------------------------------------------------------------------------------

//__Track Data Tree Specialization______________________________________________________________
class track::tree : public analysis::tree {
public:
  using real_branch_value_type = std::vector<double>;
  using real_branch_type = branch<real_branch_value_type>;

  tree(const std::string& name);
  tree(const std::string& name,
       const std::string& title);

  real_branch_type t0, x0, y0, z0, vx, vy, vz,
                   t0_error, x0_error, y0_error, z0_error, vx_error, vy_error, vz_error,
                   chi_squared, chi_squared_per_dof, chi_squared_p_value,
                   size, beta, beta_error, angle, angle_error,
                   event_t, event_x, event_y, event_z;

  branch<std::vector<std::string>> event_detector;

  branch<std::vector<unsigned long>> hash;

  void insert(const track& track);
  void clear();
  void reserve(std::size_t capacity);

  template<class UnaryPredicate>
  UnaryPredicate fill_if(const track_vector& tracks,
                         UnaryPredicate f) {
    clear();
    reserve(tracks.size());
    for (const auto& track : tracks) {
      if (f(track))
        insert(track);
    }
    analysis::tree::fill();
    return std::move(f);
  }

  void fill(const track_vector& tracks={}) {
    fill_if(tracks, [](auto) { return true; });
  }

private:
  branch<unsigned long long> _count;
  std::vector<std::reference_wrapper<real_branch_type>> _vector_branches;
};
//----------------------------------------------------------------------------------------------

//__Fit all Seeds to Tracks_____________________________________________________________________
const track_vector independent_fit_seeds(const event_vector& seeds,
                                         const Coordinate direction=Coordinate::Z);
const track_vector independent_fit_seeds(const full_event_vector& seeds,
                                         const Coordinate direction=Coordinate::Z);
//----------------------------------------------------------------------------------------------

//__Fit all Seeds to Tracks using Overlaps______________________________________________________
const track_vector overlap_fit_seeds(const event_vector& seeds,
                                     const Coordinate direction=Coordinate::Z,
                                     const std::size_t min_overlap=2UL);
const track_vector overlap_fit_seeds(const full_event_vector& seeds,
                                     const Coordinate direction=Coordinate::Z,
                                     const std::size_t min_overlap=2UL);
//----------------------------------------------------------------------------------------------

//__Remove Overlaps from Tracks_________________________________________________________________
const track_vector overlap_fit_tracks(const track_vector& tracks,
                                      const std::size_t min_overlap=2UL);
//----------------------------------------------------------------------------------------------

//__Collect Points Which are Untracked__________________________________________________________
const event non_tracked_points(const event& points,
                               const track_vector& tracks,
                               const bool ignore_diverged=false);
const full_event non_tracked_points(const full_event& points,
                                    const track_vector& tracks,
                                    const bool ignore_diverged=false);
//----------------------------------------------------------------------------------------------

// TODO: implement
// class track::graph {};

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */



#endif /* TRACKER__ANALYSIS__TRACK_HH */

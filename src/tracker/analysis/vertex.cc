/*
 * src/tracker/analysis/vertex.cc
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

#include <tracker/analysis/vertex.hh>

#include <tracker/core/stat.hh>
#include <tracker/core/units.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/error.hh>
#include <tracker/util/io.hh>
#include <tracker/util/math.hh>

#include "../helper/analysis.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Calculate Pairwise Closest Approach for Tracks in R3________________________________________
r3_point _pairwise_track_r3_closest_approach(const track& first,
                                             const track& second) {
  const auto x0 = reduce_to_r3(first.origin());
  const auto x1 = reduce_to_r3(second.origin());
  const auto v0 = first.ray();
  const auto v1 = second.ray();
  const auto v0v0 = v0 * v0;
  const auto v1v1 = v1 * v1;
  const auto v0v1 = v0 * v1;
  const auto v1v0 = v1 * v0;
  const auto n1 = v0 * v1v1 - v1 * v1v0;
  const auto n0 = v1 * v0v0 - v0 * v0v1;
  const auto diff = x1 - x0;
  // const auto c0 = x0 + (diff * n1 / (v0v0 * v1v1 - v1v0 * v1v0)) * v0;
  // const auto c1 = x1 - (diff * n0 / (v1v1 * v0v0 - v0v1 * v0v1)) * v1;
  return 0.5L * ((x1 - (diff * n0 / (v1v1 * v0v0 - v0v1 * v0v1)) * v1)
               + (x0 + (diff * n1 / (v0v0 * v1v1 - v1v0 * v1v0)) * v0));
}
//----------------------------------------------------------------------------------------------

//__Calculate Distance from Vertex to Track at Fixed Time_______________________________________
real _vertex_track_r3_distance(const real t,
                               const real x,
                               const real y,
                               const real z,
                               const track& track) {
  const auto track_point = track.at_t(t);
  return util::math::hypot(track_point.x - x, track_point.y - y, track_point.z - z);
}
//----------------------------------------------------------------------------------------------

//__Calculate Distance with Error from Vertex to Track at Fixed Time____________________________
const stat::type::uncertain_real _vertex_track_r3_distance_with_error(const real t,
                                                                      const real x,
                                                                      const real y,
                                                                      const real z,
                                                                      const track& track) {
  const auto track_point = track.at_t(t);
  const auto dx = track_point.x - x;
  const auto dy = track_point.y - y;
  const auto dz = track_point.z - z;
  const auto total_dt = t - track.t0_value();
  const auto distance = util::math::hypot(dx, dy, dz);
  const auto inverse_distance = 1.0L / distance;
  const auto dx_by_D = dx * inverse_distance;
  const auto dy_by_D = dy * inverse_distance;
  const auto dz_by_D = dz * inverse_distance;
  const real_array<6UL> gradient{
    -util::math::fused_product(track.vx_value(), dx_by_D,
                               track.vy_value(), dy_by_D,
                               track.vz_value(), dz_by_D),
    dx_by_D,
    dy_by_D,
    total_dt * dx_by_D,
    total_dt * dy_by_D,
    total_dt * dz_by_D};
  return stat::type::uncertain_real(
    distance,
    stat::error::propagate(gradient, track.covariance_matrix()));
}
//----------------------------------------------------------------------------------------------

//__Calculate Squared Residual of Vertex wrt Track______________________________________________
real _vertex_squared_residual(const stat::type::uncertain_real& distance) {
  return util::math::sum_squares(distance.value / distance.error);
}
//----------------------------------------------------------------------------------------------

//__Calculate Squared Residual of Vertex wrt Track______________________________________________
real _vertex_squared_residual(const real t,
                              const real x,
                              const real y,
                              const real z,
                              const track& track) {
  return _vertex_squared_residual(_vertex_track_r3_distance_with_error(t, x, y, z, track));
}
//----------------------------------------------------------------------------------------------

//__Fast Guess of Initial Track Parameters______________________________________________________
vertex::fit_parameters _guess_vertex(const track_vector& tracks) {
  const auto size = tracks.size();

  std::vector<full_hit> track_fronts;
  track_fronts.reserve(size);
  util::algorithm::back_insert_transform(tracks, track_fronts, [](const auto& track) {
    const auto front_t = track.t0_value();
    const auto point = track.at_t(front_t);
    const auto error = track.error_at_t(front_t);
    return full_hit{point.t, point.x, point.y, point.z,
             r4_point{track.t0_error(), error.x, error.y, error.z}};
  });

  real_vector t_errors, x_errors, y_errors, z_errors;
  t_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, t_errors,
    [](const auto& front) { return front.width.t; });
  x_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, x_errors,
    [](const auto& front) { return stat::error::uniform(front.width.x); });
  y_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, y_errors,
    [](const auto& front) { return stat::error::uniform(front.width.y); });
  z_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, z_errors,
    [](const auto& front) { return stat::error::uniform(front.width.z); });

  const auto average_point = std::accumulate(track_fronts.cbegin(), track_fronts.cend(), r4_point{},
    [](const auto sum, const auto& point) { return sum + reduce_to_r4(point); })
      / static_cast<real>(size);

  return {{average_point.t, stat::error::propagate_average(t_errors), 0, 0},
          {average_point.x, stat::error::propagate_average(x_errors), 0, 0},
          {average_point.y, stat::error::propagate_average(y_errors), 0, 0},
          {average_point.z, stat::error::propagate_average(z_errors), 0, 0}};
}
//----------------------------------------------------------------------------------------------

//__Gaussian Negative Log Likelihood Calculation________________________________________________
thread_local track_vector&& _nll_fit_tracks = {};
void _gaussian_nll(Int_t&, Double_t*, Double_t& out, Double_t* x, Int_t) {
  out = std::accumulate(_nll_fit_tracks.cbegin(), _nll_fit_tracks.cend(), 0.0L,
    [&](const auto sum, const auto& track) {
      const auto distance = _vertex_track_r3_distance_with_error(x[0], x[1], x[2], x[3], track);
      // std::cout << "distance: " << distance.value << " error: " << distance.error << "\n";
      return sum + std::fma(0.5L, _vertex_squared_residual(distance), std::log(distance.error));
    });
  // std::cout << "(" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ")\n";
  // std::cout << "NLL: " << out << "\n\n";
}
//----------------------------------------------------------------------------------------------

//__MINUIT Gaussian Fitter______________________________________________________________________
bool _fit_tracks_minuit(const track_vector& tracks,
                        vertex::fit_parameters& parameters,
                        vertex::covariance_matrix_type& covariance_matrix) {
  using namespace helper::minuit;
  auto& t = parameters.t;
  auto& x = parameters.x;
  auto& y = parameters.y;
  auto& z = parameters.z;
  _nll_fit_tracks = tracks;
  for (const auto& track : tracks)
    if (track.empty() || track.fit_diverged())
      return false;

  TMinuit minuit;
  initialize(minuit, "T", t, "X", x, "Y", y, "Z", z);

  execute(minuit, _gaussian_nll);

  get_parameters(minuit, t, x, y, z);
  get_covariance<vertex::free_parameter_count>(minuit, covariance_matrix);

  return true;
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Vertex Default Constructor__________________________________________________________________
vertex::vertex() {
  clear();
}
//----------------------------------------------------------------------------------------------

//__Vertex Constructor__________________________________________________________________________
vertex::vertex(const track_vector& tracks) {
  reset(tracks);
}
//----------------------------------------------------------------------------------------------

//__Vertex Constructor__________________________________________________________________________
vertex::vertex(track_vector&& tracks) {
  reset(std::move(tracks));
}
//----------------------------------------------------------------------------------------------

//__Vertex Point________________________________________________________________________________
const r4_point vertex::point() const {
  return r4_point{t_value(), x_value(), y_value(), z_value()};
}
//----------------------------------------------------------------------------------------------

//__Error in Calculation of Vertex Point________________________________________________________
const r4_point vertex::point_error() const {
  return r4_point{t_error(), x_error(), y_error(), z_error()};
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter from Vertex_______________________________________________________________
const fit_parameter vertex::fit_of(const vertex::parameter p) const {
  switch (p) {
    case vertex::parameter::T: return _final.t;
    case vertex::parameter::X: return _final.x;
    case vertex::parameter::Y: return _final.y;
    case vertex::parameter::Z: return _final.z;
  }
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter Value from Vertex_________________________________________________________
real vertex::value(const vertex::parameter p) const {
  return fit_of(p).value;
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter Error from Vertex_________________________________________________________
real vertex::error(const vertex::parameter p) const {
  return fit_of(p).error;
}
//----------------------------------------------------------------------------------------------

//__Check if Track Fit Converged________________________________________________________________
bool vertex::fit_diverged() const noexcept {
  return _final == vertex::fit_parameters{};
}
//----------------------------------------------------------------------------------------------

//__Get Distance from Each Track at Time of Vertex______________________________________________
real_vector vertex::distances() const {
  real_vector out;
  out.reserve(_tracks.size());
  util::algorithm::back_insert_transform(_tracks, out, [&](const auto& track) {
    return _vertex_track_r3_distance(
      _final.t.value, _final.x.value, _final.y.value, _final.z.value, track); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Get Error in Distance of Each Track at Time of Vertex_______________________________________
real_vector vertex::distance_errors() const {
  real_vector out;
  out.reserve(_tracks.size());
  util::algorithm::back_insert_transform(_tracks, out, [&](const auto& track) {
    return _vertex_track_r3_distance_with_error(
      _final.t.value, _final.x.value, _final.y.value, _final.z.value, track).error; });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared Test Statistic__________________________________________________________________
real vertex::chi_squared() const {
  return std::accumulate(_delta_chi2.cbegin(), _delta_chi2.cend(), 0.0L);
}
//----------------------------------------------------------------------------------------------

//__Vertex Degrees of Freedom___________________________________________________________________
std::size_t vertex::degrees_of_freedom() const {
  return 4UL;
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared per Degree of Freedom___________________________________________________________
real vertex::chi_squared_per_dof() const {
  return chi_squared() / degrees_of_freedom();
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared P-Value_________________________________________________________________________
real vertex::chi_squared_p_value() const {
  return stat::chi_squared_p_value(chi_squared(), degrees_of_freedom());
}
//----------------------------------------------------------------------------------------------

//__Get Variance of a Vertex Parameter__________________________________________________________
real vertex::variance(const vertex::parameter p) const {
  return covariance(p, p);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Get Shift Index of Vertex Parameters for Covariance Matrix__________________________________
constexpr std::size_t _shift_covariance_index(const vertex::parameter p) {
  switch (p) {
    case vertex::parameter::T: return 0;
    case vertex::parameter::X: return 1;
    case vertex::parameter::Y: return 2;
    case vertex::parameter::Z: return 3;
  }
}
//----------------------------------------------------------------------------------------------
} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Get Covariance between Vertex Parameters____________________________________________________
real vertex::covariance(const vertex::parameter p,
                        const vertex::parameter q) const {
  return _covariance[4 * _shift_covariance_index(p) + _shift_covariance_index(q)];
}
//----------------------------------------------------------------------------------------------

//__Reset Vertex with Given Tracks______________________________________________________________
std::size_t vertex::reset(const track_vector& tracks) {
  _tracks = tracks;
  const auto new_size = size();

  _delta_chi2.clear();
  _delta_chi2.reserve(new_size);

  if (new_size > 1) {
    _guess = _guess_vertex(_tracks);
    _final = _guess;
    if (_fit_tracks_minuit(_tracks, _final, _covariance)) {
      util::algorithm::back_insert_transform(_tracks, _delta_chi2,
        [&](const auto& track) {
          return _vertex_squared_residual(
            _final.t.value,
            _final.x.value,
            _final.y.value,
            _final.z.value,
            track);
        });
      return new_size;
    }
  } else {
    _guess = {};
  }

  _final = {};
  _covariance = {};
  _delta_chi2.resize(new_size, 0);
  return size();
}
//----------------------------------------------------------------------------------------------

//__Insert Track into Vertex and Refit__________________________________________________________
std::size_t vertex::insert(const track& track) {
  if (std::find(_tracks.cbegin(), _tracks.cend(), track) != _tracks.cend()) {
    _tracks.push_back(track);
    _tracks.shrink_to_fit();
    return reset(_tracks);
  }
  return size();
}
//----------------------------------------------------------------------------------------------

//__Insert Tracks into Vertex and Refit_________________________________________________________
std::size_t vertex::insert(const track_vector& tracks) {
  _tracks.reserve(size() + tracks.size());
  const auto begin = _tracks.cbegin();
  std::copy_if(tracks.cbegin(), tracks.cend(), std::back_inserter(_tracks),
    [&](const auto& t) { return std::find(begin, _tracks.cend(), t) != _tracks.cend(); });
  _tracks.shrink_to_fit();
  return reset(_tracks);
}
//----------------------------------------------------------------------------------------------

//__Remove Track from Vertex and Refit__________________________________________________________
std::size_t vertex::remove(const std::size_t index) {
  // TODO: improve efficiency
  const auto s = size();
  if (index >= s)
    return s;

  track_vector saved_tracks;
  saved_tracks.reserve(s - 1);
  const auto begin = _tracks.cbegin();
  const auto end = _tracks.cend();
  saved_tracks.insert(saved_tracks.cend(), begin, begin + index);
  saved_tracks.insert(saved_tracks.cend(), begin + index + 1, end);
  return reset(saved_tracks);
}
//----------------------------------------------------------------------------------------------

//__Remove Tracks from Vertex and Refit_________________________________________________________
std::size_t vertex::remove(const std::vector<std::size_t>& indices) {
  // TODO: improve efficiency
  const auto sorted = util::algorithm::copy_sort_range(indices);
  const auto s = size();
  track_vector saved_tracks;
  saved_tracks.reserve(s);
  for (std::size_t hit_index{}, removal_index{}; hit_index < s; ++hit_index) {
    if (sorted[removal_index] == hit_index) {
      ++removal_index;
    } else {
      saved_tracks.push_back(std::move(_tracks[hit_index]));
    }
  }
  saved_tracks.shrink_to_fit();
  return reset(saved_tracks);
}
//----------------------------------------------------------------------------------------------

//__Remove Track from Vertex Below Maximum Chi-Squared and Refit________________________________
std::size_t vertex::prune_on_chi_squared(const real max_chi_squared) {
  const auto s = size();
  std::vector<std::size_t> indices;
  indices.reserve(s);
  for (std::size_t i{}; i < s; ++i) {
    if (chi_squared_vector()[i] > max_chi_squared)
      indices.push_back(i);
  }
  return remove(indices);
}
//----------------------------------------------------------------------------------------------

//__Clear Tracks from Vertex____________________________________________________________________
void vertex::clear() {
  reset({});
}
//----------------------------------------------------------------------------------------------

//__Fill Plots with Vertexing Variables_________________________________________________________
void vertex::fill_plots(plot::histogram_collection& collection,
                        const vertex::plotting_keys& keys) const {
  if (collection.count(keys.t)) collection[keys.t].insert(t_value() / units::time);
  if (collection.count(keys.x)) collection[keys.x].insert(x_value() / units::length);
  if (collection.count(keys.y)) collection[keys.y].insert(y_value() / units::length);
  if (collection.count(keys.z)) collection[keys.z].insert(z_value() / units::length);
  if (collection.count(keys.t_error)) collection[keys.t_error].insert(t_error() / units::time);
  if (collection.count(keys.x_error)) collection[keys.x_error].insert(x_error() / units::length);
  if (collection.count(keys.y_error)) collection[keys.y_error].insert(y_error() / units::length);
  if (collection.count(keys.z_error)) collection[keys.z_error].insert(z_error() / units::length);

  if (collection.count(keys.distance)) {
    auto& distance_histogram = collection[keys.distance];
    for (const auto& distance : distances())
      distance_histogram.insert(distance / units::length);
  }
  if (collection.count(keys.distance_error)) {
    auto& distance_error_histogram = collection[keys.distance_error];
    for (const auto& distance_error : distance_errors())
      distance_error_histogram.insert(distance_error / units::length);
  }

  if (collection.count(keys.chi_squared)) collection[keys.chi_squared].insert(chi_squared());
  if (collection.count(keys.chi_squared_per_dof)) collection[keys.chi_squared_per_dof].insert(chi_squared_per_dof());
  if (collection.count(keys.chi_squared_p_value)) collection[keys.chi_squared_p_value].insert(chi_squared_p_value());
  if (collection.count(keys.size)) collection[keys.size].insert(size());
}
//----------------------------------------------------------------------------------------------

//__Draw Fit Vertex_____________________________________________________________________________
void vertex::draw(plot::canvas& canvas,
                  const real size,
                  const plot::color color,
                  const bool with_errors) const {
  // TODO: decide what to do with convergence
  if (fit_converged()) {
    canvas.add_point(point(), size, color);
    if (with_errors)
      canvas.add_box(point(), point_error().x, point_error().y, point_error().z, size, color);
  }
}
//----------------------------------------------------------------------------------------------

//__Draw Guess Vertex___________________________________________________________________________
void vertex::draw_guess(plot::canvas& canvas,
                        const real size,
                        const plot::color color,
                        const bool with_errors) const {
  // TODO: add with_errors argument
  canvas.add_point(r3_point{guess_fit().x.value, guess_fit().y.value, guess_fit().z.value}, size, color);
}
//----------------------------------------------------------------------------------------------

//__Hash Implementation_________________________________________________________________________
std::size_t vertex::hash() const {
  return util::functional::hash_combine_range(_tracks, 0x13c9ee801bULL);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Print Vertex Parameters with Units__________________________________________________________
std::ostream& _print_vertex_parameters(std::ostream& os,
                                       const vertex::fit_parameters& parameters,
                                       std::size_t prefix_count) {
  return os
    << std::string(prefix_count, ' ')
      << "T: " << std::setw(10) << parameters.t.value / units::time
               << "  (+/- " << std::setw(10) << parameters.t.error / units::time     << ")  "
               << units::time_string << "\n"
    << std::string(prefix_count, ' ')
      << "X: " << std::setw(10) << parameters.x.value / units::length
               << "  (+/- " << std::setw(10) << parameters.x.error / units::length   << ")  "
               << units::length_string << "\n"
    << std::string(prefix_count, ' ')
      << "Y: " << std::setw(10) << parameters.y.value / units::length
               << "  (+/- " << std::setw(10) << parameters.y.error / units::length   << ")  "
               << units::length_string << "\n"
    << std::string(prefix_count, ' ')
      << "Z: " << std::setw(10) << parameters.z.value / units::length
               << "  (+/- " << std::setw(10) << parameters.z.error / units::length   << ")  "
               << units::length_string << "\n";
}
//----------------------------------------------------------------------------------------------
} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Vertex Output Stream Operator_______________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const vertex& vertex) {
  static const std::string bar(80, '-');
  os << bar << "\n";
  os.precision(6);

  if (vertex.fit_diverged()) {

    os << "* Vertex Status: " << util::io::bold << "DIVERGED" << util::io::reset_font << "\n"
       << "* Guess Parameters:\n";
    _print_vertex_parameters(os, vertex.guess_fit(), 4);

    os << "* Tracks: \n";
    for (const auto& track : vertex.tracks()) {
      os << "    (" << track.t0_value() / units::time     << ", "
                    << track.x0_value() / units::length   << ", "
                    << track.y0_value() / units::length   << ", "
                    << track.z0_value() / units::length   << ", "
                    << track.vx_value() / units::velocity << ", "
                    << track.vy_value() / units::velocity << ", "
                    << track.vz_value() / units::velocity << ")\n";
    }

  } else {

    os << "* Vertex Status: " << util::io::bold << "CONVERGED" << util::io::reset_font << "\n"
       << "* Parameters:\n";
    _print_vertex_parameters(os, vertex.final_fit(), 4);

    os << "* Tracks: \n";
    const auto size = vertex.size();
    const auto& tracks = vertex.tracks();
    const auto& distances = vertex.distances();
    const auto& errors = vertex.distance_errors();
    for (std::size_t i{}; i < size; ++i) {
      const auto& track = tracks[i];
      os << "    " << distances[i] / units::length
                   << "  (+/- " << errors[i] / units::length << ") "
                   << units::length_string << "\n      from ("
                   << track.t0_value() / units::time     << ", "
                   << track.x0_value() / units::length   << ", "
                   << track.y0_value() / units::length   << ", "
                   << track.z0_value() / units::length   << ", "
                   << track.vx_value() / units::velocity << ", "
                   << track.vy_value() / units::velocity << ", "
                   << track.vz_value() / units::velocity << ")\n";
    }

    os << "* Statistics: \n"
       << "    dof:      " << vertex.degrees_of_freedom()             << "\n"
       << "    chi2:     " << vertex.chi_squared() << " = ";
    util::io::print_range(vertex.chi_squared_vector(), " + ", "", os) << "\n";
    os << "    chi2/dof: " << vertex.chi_squared_per_dof()            << "\n"
       << "    p-value:  " << vertex.chi_squared_p_value()            << "\n"
       << "    cov mat:  | ";
    const auto matrix = vertex.covariance_matrix();
    os << std::right;
    for (size_t i = 0; i < 4; ++i) {
      if (i > 0) os << "              | ";
      for (size_t j = 0; j < 4; ++j) {
        const auto cell = matrix[4*i+j];
        real cell_unit{1.0L};
        if (i == 0) cell_unit *= units::time;
        else cell_unit *= units::length;
        if (j == 0) cell_unit *= units::time;
        else cell_unit *= units::length;

        if (i == j) {
          os << util::io::bold << std::setw(14)
             << cell / cell_unit << util::io::reset_font << " ";
        } else {
          os << std::setw(14) << cell / cell_unit << " ";
        }
      }
      os << "|\n";
    }
  }

  return os << bar;
}
//----------------------------------------------------------------------------------------------

//__Vertex Data Tree Constructor________________________________________________________________
vertex::tree::tree(const std::string& name)
    : tree(name, name) {}
//----------------------------------------------------------------------------------------------

//__Vertex Data Tree Constructor________________________________________________________________
vertex::tree::tree(const std::string& name,
                   const std::string& title)
    : analysis::tree(name, title),
      t(emplace_branch<real_branch_value_type>("t")),
      x(emplace_branch<real_branch_value_type>("x")),
      y(emplace_branch<real_branch_value_type>("y")),
      z(emplace_branch<real_branch_value_type>("z")),
      t_error(emplace_branch<real_branch_value_type>("t_error")),
      x_error(emplace_branch<real_branch_value_type>("x_error")),
      y_error(emplace_branch<real_branch_value_type>("y_error")),
      z_error(emplace_branch<real_branch_value_type>("z_error")),
      chi_squared(emplace_branch<real_branch_value_type>("chi_squared")),
      chi_squared_per_dof(emplace_branch<real_branch_value_type>("chi_squared_per_dof")),
      chi_squared_p_value(emplace_branch<real_branch_value_type>("chi_squared_p_value")),
      size(emplace_branch<real_branch_value_type>("size")),
      track_hash(emplace_branch<decltype(track_hash)::value_type>("track_hash")),
      hash(emplace_branch<decltype(hash)::value_type>("hash")),
      _count(emplace_branch<decltype(_count)::value_type>("N")),
      _vector_branches({t, x, y, z,
                        t_error, x_error, y_error, z_error,
                        chi_squared, chi_squared_per_dof, chi_squared_p_value,
                        size}) {}
//----------------------------------------------------------------------------------------------

//__Track Data Tree Insertion___________________________________________________________________
void vertex::tree::insert(const vertex& vertex) {
  t.get().push_back(vertex.t_value() / units::time);
  x.get().push_back(vertex.x_value() / units::length);
  y.get().push_back(vertex.y_value() / units::length);
  z.get().push_back(vertex.z_value() / units::length);
  t_error.get().push_back(vertex.t_error() / units::time);
  x_error.get().push_back(vertex.x_error() / units::length);
  y_error.get().push_back(vertex.y_error() / units::length);
  z_error.get().push_back(vertex.z_error() / units::length);
  chi_squared.get().push_back(vertex.chi_squared());
  chi_squared_per_dof.get().push_back(vertex.chi_squared_per_dof());
  chi_squared_p_value.get().push_back(vertex.chi_squared_p_value());
  size.get().push_back(vertex.size());
  for (const auto& track : vertex.tracks())
    track_hash.get().push_back(track.hash());
  hash.get().push_back(vertex.hash());
  ++_count;
}
//----------------------------------------------------------------------------------------------

//__Clear Vertex Data Tree______________________________________________________________________
void vertex::tree::clear() {
  _count = 0UL;
  for (auto& entry : _vector_branches)
    entry.get().get().clear();
  track_hash.get().clear();
  hash.get().clear();
}
//----------------------------------------------------------------------------------------------

//__Reserve Space for Vertex Data Tree__________________________________________________________
void vertex::tree::reserve(std::size_t capacity) {
  for (auto& entry : _vector_branches)
    entry.get().get().reserve(capacity);
  hash.get().reserve(capacity);
}
//----------------------------------------------------------------------------------------------

//__Pairwise Fit Tracks to Vertices_____________________________________________________________
const vertex_vector pairwise_fit_tracks(const track_vector& tracks) {
  const auto size = tracks.size();
  if (size == 0UL)
    return vertex_vector{};
  if (size <= 2UL)
    return vertex_vector{vertex{tracks}};

  vertex_vector out;
  out.reserve(size);

  util::bit_vector join_list(size);

  vertex current_vertex{};
  track_vector pair;
  pair.reserve(2UL);

  std::size_t top_index{}, bottom_index{1UL};
  while (top_index < size) {
    pair.push_back(tracks[top_index]);
    bottom_index = join_list.first_unset(bottom_index);

    real best_fit{-1.0L};
    std::size_t best_bottom_index{bottom_index};
    while (bottom_index < size) {
      pair.push_back(tracks[bottom_index]);
      current_vertex.reset(pair);

      if (current_vertex.fit_converged()) {
        if (best_fit < 0.0L) {
          best_fit = current_vertex.chi_squared_per_dof();
        } else {
          const auto current_fit = current_vertex.chi_squared_per_dof();
          if (current_fit < best_fit) {
            best_fit = current_fit;
            best_bottom_index = bottom_index;
          }
        }
      }

      bottom_index = join_list.first_unset(1UL + bottom_index);
      pair.pop_back();
    }

    if (best_bottom_index < size) {
      pair.push_back(tracks[best_bottom_index]);
      out.emplace_back(pair);
    }

    join_list.set(top_index);
    top_index = join_list.first_unset(1UL + top_index);
    bottom_index = 1UL + top_index;
    pair.clear();
  }

  out.shrink_to_fit();
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

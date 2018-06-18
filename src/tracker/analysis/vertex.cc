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

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Vertex Parameter Type_______________________________________________________________________
struct _vertex_parameters { fit_parameter t, x, y, z; };
//----------------------------------------------------------------------------------------------

//__Fast Guess of Initial Track Parameters______________________________________________________
_vertex_parameters _guess_vertex(const track_vector& tracks) {
  return {};
}
//----------------------------------------------------------------------------------------------

//__Gaussian Negative Log Likelihood Calculation________________________________________________
thread_local track_vector&& _nll_fit_tracks = {};
void _gaussian_nll(Int_t&, Double_t*, Double_t& out, Double_t* x, Int_t) {
  out = 0;
}
//----------------------------------------------------------------------------------------------

//__MINUIT Gaussian Fitter______________________________________________________________________
void _fit_tracks_minuit(const track_vector& tracks,
                        _vertex_parameters& parameters,
                        real_vector& covariance_matrix) {
  TMinuit minuit;
  //minuit.SetGraphicsMode(settings.graphics_on);
  //minuit.SetPrintLevel(settings.print_level);
  //minuit.SetErrorDef(settings.error_def);
  //minuit.SetMaxIterations(settings.max_iterations);

  minuit.Command("SET STR 2");

  auto& t = parameters.t;
  minuit.DefineParameter(0, "T", t.value, t.error, t.min, t.max);
  auto& x = parameters.x;
  minuit.DefineParameter(1, "X", x.value, x.error, x.min, x.max);
  auto& y = parameters.y;
  minuit.DefineParameter(2, "Y", y.value, y.error, y.min, y.max);
  auto& z = parameters.z;
  minuit.DefineParameter(3, "Z", z.value, z.error, z.min, z.max);

  _nll_fit_tracks = tracks;
  minuit.SetFCN(_gaussian_nll);

  /*
  Int_t error_flag;
  auto command_parameters = settings.command_parameters;
  minuit.mnexcm(
    settings.command_name.c_str(),
    command_parameters.data(),
    command_parameters.size(),
    error_flag);

  switch (error_flag) {
    case 1:
    case 2:
    case 3: util::error::exit("[FATAL ERROR] Unknown MINUIT Command \"", settings.command_name,
                              "\". Exited with Error Code ", error_flag, ".\n");
    //case 4: util::error::exit("[FATAL ERROR] MINUIT Exited Abnormally ",
    //                          "with Error Code ", error_flag, ".\n");
    default: break;
  }
  */

  Double_t value, error;
  minuit.GetParameter(0, value, error);
  t.value = value;
  t.error = error;
  minuit.GetParameter(1, value, error);
  x.value = value;
  x.error = error;
  minuit.GetParameter(2, value, error);
  y.value = value;
  y.error = error;
  minuit.GetParameter(3, value, error);
  z.value = value;
  z.error = error;

  constexpr const std::size_t dimension = 4;
  Double_t matrix[dimension][dimension];
  minuit.mnemat(&matrix[0][0], dimension);
  covariance_matrix.clear();
  covariance_matrix.reserve(dimension * dimension);
  for (std::size_t i = 0; i < dimension; ++i)
    for (std::size_t j = 0; j < dimension; ++j)
      covariance_matrix.push_back(matrix[i][j]);
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Vertex Constructor__________________________________________________________________________
vertex::vertex(const track_vector& tracks) : _tracks(tracks) {
  auto fit_vertex = _guess_vertex(_tracks);
  _fit_tracks_minuit(_tracks, fit_vertex, _covariance);

  _t = std::move(fit_vertex.t);
  _x = std::move(fit_vertex.x);
  _y = std::move(fit_vertex.y);
  _z = std::move(fit_vertex.z);

}
//----------------------------------------------------------------------------------------------

//__Vertex Point________________________________________________________________________________
const r4_point vertex::point() const {
  return {_t.value, _x.value, _y.value, _z.value};
}
//----------------------------------------------------------------------------------------------

//__Error in Calculation of Vertex Point________________________________________________________
const r4_point vertex::point_error() const {
  return {_t.error, _x.error, _y.error, _z.error};
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter from Vertex_______________________________________________________________
const fit_parameter vertex::fit_of(const vertex::parameter p) const {
  switch (p) {
    case vertex::parameter::T: return _t;
    case vertex::parameter::X: return _x;
    case vertex::parameter::Y: return _y;
    case vertex::parameter::Z: return _z;
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
     << "    p-value:  " << stat::chi_squared_p_value(vertex)         << "\n"
     << "    cov mat:  | ";
  const auto matrix = vertex.covariance_matrix();
  for (size_t i = 0; i < 4; ++i) {
    if (i > 0) os << "              | ";
    for (size_t j = 0; j < 4; ++j) {
      const auto cell = matrix[4*i+j];
      if (i == j) {
        os << util::io::bold << util::io::underline
           << cell << util::io::reset_font << " ";
      } else {
        os << cell << " ";
      }
    }
    os << "|\n";
  }

  return os << bar;
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

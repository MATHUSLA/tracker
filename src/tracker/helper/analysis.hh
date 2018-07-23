/*
 * src/tracker/helper/analysis.hh
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

#ifndef TRACKER__HELPER__ANALYSIS_HH
#define TRACKER__HELPER__ANALYSIS_HH
#pragma once

#include <ROOT/TMinuit.h>

#include <tracker/analysis.hh>

#include <tracker/util/error.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { namespace helper { ////////////////////////////////////////////////////////

namespace minuit { /////////////////////////////////////////////////////////////////////////////

namespace settings { ///////////////////////////////////////////////////////////////////////////

//__MINUIT Minimization Details_________________________________________________________________
static const std::string&        command        = "MIGRAD";
static const std::vector<double> parameters     = {};
static const bool                graphics       = false;
static const integer             print_level    = -1;
static const double              error_def      = 0.5;
static const integer             max_iterations = 1000;
static const std::string&        strategy       = "SET STR 2";
//----------------------------------------------------------------------------------------------

} /* namespace settings */ /////////////////////////////////////////////////////////////////////

namespace error { //////////////////////////////////////////////////////////////////////////////

//__MINUIT Minimization Error Codes_____________________________________________________________
// TODO: populate with better names resembling the proper error code
static const int unknown_1 = 1;
static const int unknown_2 = 2;
static const int unknown_3 = 3;
static const int diverged  = 4;
//----------------------------------------------------------------------------------------------

/* TODO: use enum instead of above
enum error {
  diverged
};
 */

} /* namespace error */ ////////////////////////////////////////////////////////////////////////

//__Function to Minimize________________________________________________________________________
using minimizing_function = void(Int_t&, Double_t*, Double_t&, Double_t*, Int_t);
//----------------------------------------------------------------------------------------------

//__Prepare MINUIT Environment__________________________________________________________________
inline TMinuit& initialize(TMinuit& minimizer) {
  minimizer.SetGraphicsMode(settings::graphics);
  minimizer.SetPrintLevel(settings::print_level);
  minimizer.SetErrorDef(settings::error_def);
  minimizer.SetMaxIterations(settings::max_iterations);
  minimizer.Command(settings::strategy.c_str());
  return minimizer;
}
//----------------------------------------------------------------------------------------------

//__Set Parameters to Minimization______________________________________________________________
inline TMinuit& set_parameters(TMinuit& minimizer,
                               const std::size_t starting_index,
                               const std::string& name,
                               const fit_parameter& parameter) {
  minimizer.DefineParameter(starting_index, name.c_str(),
    parameter.value, parameter.error, parameter.min, parameter.max);
  return minimizer;
}
template<class ...Args>
TMinuit& set_parameters(TMinuit& minimizer,
                        const std::size_t starting_index,
                        const std::string& name,
                        const fit_parameter& parameter,
                        const Args& ...rest) {
  set_parameters(minimizer, starting_index, name, parameter);
  set_parameters(minimizer, starting_index + 1, rest...);
  return minimizer;
}
template<class ...Args>
TMinuit& set_parameters(TMinuit& minimizer,
                        const std::string& name,
                        const fit_parameter& parameter,
                        const Args& ...rest) {
  return set_parameters(minimizer, 0UL, name, parameter, rest...);
}
//----------------------------------------------------------------------------------------------

//__Prepare MINUIT Environment and Set Parameters_______________________________________________
template<class ...Args>
TMinuit& initialize(TMinuit& minimizer,
                    const Args& ...parameters) {
  return set_parameters(initialize(minimizer), parameters...);
}
//----------------------------------------------------------------------------------------------

//__Execute Minimization________________________________________________________________________
inline int execute(TMinuit& minimizer,
                   minimizing_function& f,
                   bool exit_on_error=false) {
  Int_t error_flag;
  auto parameters_copy = settings::parameters;
  minimizer.SetFCN(f);
  minimizer.mnexcm(
    settings::command.c_str(),
    parameters_copy.data(),
    parameters_copy.size(),
    error_flag);

  if (exit_on_error) {
    switch (error_flag) {
      case 0: break;
      case 1:
      case 2:
      case 3: util::error::exit(error_flag,
                                "[FATAL ERROR] Unknown MINUIT Command \"", settings::command,
                                "\". Exited with Error Code ", error_flag, ".\n");
      case 4: util::error::exit(error_flag,
                                "[FATAL ERROR] MINUIT Exited Abnormally ",
                                "with Error Code ", error_flag, ".\n");
      default: util::error::exit(error_flag,
                                 "[FATAL ERROR] Unknown MINUIT Error Code: ", error_flag, ".\n");
    }
  }
  return error_flag;
}
//----------------------------------------------------------------------------------------------

//__Get Parameters from Minimization____________________________________________________________
inline void get_parameters(const TMinuit& minimizer,
                           const std::size_t starting_index,
                           fit_parameter& parameter) {
  Double_t value, error;
  minimizer.GetParameter(starting_index, value, error);
  parameter.value = value;
  parameter.error = error;
}
template<class ...Args>
void get_parameters(const TMinuit& minimizer,
                    const std::size_t starting_index,
                    fit_parameter& parameter,
                    Args& ...rest) {
  get_parameters(minimizer, starting_index, parameter);
  get_parameters(minimizer, starting_index + 1, rest...);
}
template<class ...Args>
void get_parameters(const TMinuit& minimizer,
                    fit_parameter& parameter,
                    Args& ...rest) {
  get_parameters(minimizer, 0, parameter, rest...);
}
//----------------------------------------------------------------------------------------------

//__Get Covariance Matrix From Minimization_____________________________________________________
template<std::size_t N>
void get_covariance(TMinuit& minimizer,
                    real_array<N * N>& covariance_matrix) {
  Double_t matrix[N][N];
  minimizer.mnemat(&matrix[0][0], N);
  for (std::size_t i = 0; i < N; ++i)
    for (std::size_t j = 0; j < N; ++j)
      covariance_matrix[N*i+j] = matrix[i][j];
}
//----------------------------------------------------------------------------------------------

} /* namespace minuit */ ///////////////////////////////////////////////////////////////////////

} } /* namespace analysis::helper */ ///////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__HELPER__ANALYSIS_HH */

/*
 * src/tracker/core/stat.cc
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

#include <tracker/core/stat.hh>

#include <ROOT/Math/ProbFunc.h>

namespace MATHUSLA { namespace TRACKER {

namespace stat { ///////////////////////////////////////////////////////////////////////////////

//__Calculate P-Value from Chi^2________________________________________________________________
real chi_squared_p_value(const real chi2,
                         const std::size_t dof) {
  return ROOT::Math::chisquared_cdf_c(chi2, dof);
}
//----------------------------------------------------------------------------------------------

} /* namespace stat */ /////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

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

namespace random { /////////////////////////////////////////////////////////////////////////////

//__Min Element of Distribution_________________________________________________________________
real generator::min() const {
  real out;
  detail::distribution_apply(_type, _distribution, [&](const auto& dist) { out = dist.min(); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Max Element of Distribution_________________________________________________________________
real generator::max() const {
  real out;
  detail::distribution_apply(_type, _distribution, [&](const auto& dist) { out = dist.max(); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Calculate Next Random Number________________________________________________________________
real generator::operator()() {
  real out;
  detail::distribution_apply(_type, _distribution, [&](auto& dist) { out = dist(_engine); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Seed Engine_________________________________________________________________________________
void generator::seed(const std::uint_least32_t seed) {
  _engine.seed(seed);
}
//----------------------------------------------------------------------------------------------

//__Seed Engine_________________________________________________________________________________
void generator::seed(std::seed_seq& seq) {
  _engine.seed(seq);
}
//----------------------------------------------------------------------------------------------

//__Reset Distribution__________________________________________________________________________
void generator::reset() {
  detail::distribution_apply(_type, _distribution, [](auto& dist) { dist.reset(); });
}
//----------------------------------------------------------------------------------------------

//__Discard Elements from Engine________________________________________________________________
void generator::discard(const std::size_t z) {
  _engine.discard(z);
}
//----------------------------------------------------------------------------------------------

//__Random Number Generator Equality____________________________________________________________
bool generator::operator==(const generator& other) const {
  return _engine == other._engine && _type == other._type;
  // TODO: add distribution equality
  // TODO: add parameter equality
}
//----------------------------------------------------------------------------------------------

} /* namespace random */ ///////////////////////////////////////////////////////////////////////

} /* namespace stat */ /////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

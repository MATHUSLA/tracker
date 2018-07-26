/*
 * include/tracker/core/stat.hh
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

#ifndef TRACKER__CORE__STAT_HH
#define TRACKER__CORE__STAT_HH
#pragma once

#include <random>

#include <tracker/core/type.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/math.hh>
#include <tracker/util/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace stat { ///////////////////////////////////////////////////////////////////////////////

using namespace type;

namespace type { ///////////////////////////////////////////////////////////////////////////////

//__Chi^2 Type Traits___________________________________________________________________________
template<class C, typename = void>
struct has_chi_squared_member : std::false_type {};
template<class C>
struct has_chi_squared_member<C, decltype(C::chi_squared, void())> : std::true_type {};
template<class C>
constexpr bool has_chi_squared_member_v = has_chi_squared_member<C>::value;

template<class C, typename = void>
struct has_chi_squared_method : std::false_type {};
template<class C>
struct has_chi_squared_method<C, decltype(&C::chi_squared, void())> : std::true_type {};
template<class C>
constexpr bool has_chi_squared_method_v = has_chi_squared_method<C>::value;
//----------------------------------------------------------------------------------------------

//__Degree of Freedom Type Traits_______________________________________________________________
template<class C, typename = void>
struct has_degrees_of_freedom_member : std::false_type {};
template<class C>
struct has_degrees_of_freedom_member<C, decltype(C::degrees_of_freedom, void())> : std::true_type {};
template<class C>
constexpr bool has_degrees_of_freedom_member_v = has_degrees_of_freedom_member<C>::value;

template<class C, typename = void>
struct has_degrees_of_freedom_method : std::false_type {};
template<class C>
struct has_degrees_of_freedom_method<C, decltype(&C::degrees_of_freedom, void())> : std::true_type {};
template<class C>
constexpr bool has_degrees_of_freedom_method_v = has_degrees_of_freedom_method<C>::value;
//----------------------------------------------------------------------------------------------

//__Chi^2 and DOF Type Reflection_______________________________________________________________
template<class C,
  bool = has_chi_squared_member_v<C>
      && has_degrees_of_freedom_member_v<C>>
struct has_chi2_and_dof_members : std::true_type {};
template<class C>
struct has_chi2_and_dof_members<C, false> : std::false_type {};
template<class C>
constexpr bool has_chi2_and_dof_members_v = has_chi2_and_dof_members<C>::value;

template<class C,
  bool = has_chi_squared_method_v<C>
      && has_degrees_of_freedom_method_v<C>>
struct has_chi2_and_dof_methods : std::true_type {};
template<class C>
struct has_chi2_and_dof_methods<C, false> : std::false_type {};
template<class C>
constexpr bool has_chi2_and_dof_methods_v = has_chi2_and_dof_methods<C>::value;

template<class C,
  bool =  has_chi2_and_dof_methods_v<C>
      ||  has_chi2_and_dof_members_v<C>
      || (has_chi_squared_member_v<C> && has_degrees_of_freedom_method_v<C>)
      || (has_chi_squared_method_v<C> && has_degrees_of_freedom_member_v<C>)>
struct is_chi2_dof_type : std::true_type {};
template<class C>
struct is_chi2_dof_type<C, false> : std::false_type {};
template<class C>
constexpr bool is_chi2_dof_type_v = is_chi2_dof_type<C>::value;
//----------------------------------------------------------------------------------------------

} /* namespace type */ /////////////////////////////////////////////////////////////////////////

//__Free Function Alternatives to Memeber Functions_____________________________________________
template<class T>
real chi_squared(const T& t);
template<class T>
real degrees_of_freedom(const T& t);
//----------------------------------------------------------------------------------------------

//__Calculate P-Value from Chi^2________________________________________________________________
real chi_squared_p_value(const real chi2,
                         const std::size_t dof);
//----------------------------------------------------------------------------------------------

//__Calculate P-Value from Chi^2________________________________________________________________
template<class T>
std::enable_if_t<!type::is_chi2_dof_type_v<T>, real>
chi_squared_p_value(const T& t) {
  return chi_squared_p_value(chi_squared(t), degrees_of_freedom(t));
}
template<class T>
std::enable_if_t<type::has_chi2_and_dof_methods_v<T>, real>
chi_squared_p_value(const T& t) {
  return chi_squared_p_value(t.chi_squared(), t.degrees_of_freedom());
}
template<class T>
std::enable_if_t<type::has_chi2_and_dof_members_v<T>, real>
chi_squared_p_value(const T& t) {
  return chi_squared_p_value(t.chi_squared, t.degrees_of_freedom);
}
//----------------------------------------------------------------------------------------------

//__Perform Chi^2/DOF Cut on Range______________________________________________________________
template<class Range>
std::enable_if_t<!type::is_chi2_dof_type_v<typename Range::value_type>, Range>&
chi2_per_dof_cut(const Range& range,
                 const real min,
                 const real max,
                 Range& out) {
  util::algorithm::back_insert_copy_if(range, out, [&](const auto& value) {
    return util::algorithm::between(chi_squared(value) / degrees_of_freedom(value), min, max); });
  return out;
}
template<class Range>
std::enable_if_t<type::has_chi2_and_dof_methods_v<typename Range::value_type>, Range>&
chi2_per_dof_cut(const Range& range,
                 const real min,
                 const real max,
                 Range& out) {
  util::algorithm::back_insert_copy_if(range, out, [&](const auto& value) {
    return util::algorithm::between(value.chi_squared() / value.degrees_of_freedom(), min, max); });
  return out;
}
template<class Range>
std::enable_if_t<type::has_chi2_and_dof_members_v<typename Range::value_type>, Range>&
chi2_per_dof_cut(const Range& range,
                 const real min,
                 const real max,
                 Range& out) {
  util::algorithm::back_insert_copy_if(range, out, [&](const auto& value) {
    return util::algorithm::between(value.chi_squared / value.degrees_of_freedom, min, max); });
  return out;
}
//----------------------------------------------------------------------------------------------

namespace error { //////////////////////////////////////////////////////////////////////////////

//__Propagate Error_____________________________________________________________________________
template<std::size_t N>
real propagate(const real_array<N>& gradient,
               const real_array<N*N>& covariance) {
  return std::sqrt(weighted_norm(gradient, covariance));
}
template<std::size_t N>
real propagate(const real_vector& gradient,
               const real_vector& covariance) {
  return propagate(to_array<N>(gradient), to_array<N*N>(covariance));
}
//----------------------------------------------------------------------------------------------

//__Propagate Error in Sum of Independent Errors________________________________________________
template<class ...Args>
constexpr real propagate_sum(const real error,
                             const Args... rest) {
  return util::math::hypot(error, rest...);
}
//----------------------------------------------------------------------------------------------

//__Propagate Error in Sum of Independent Errors________________________________________________
template<class Range>
constexpr real propagate_sum(const Range& range) {
  return util::math::range_hypot(range);
}
//----------------------------------------------------------------------------------------------

//__Propagate Error in Average of Independent Errors____________________________________________
template<class ...Args>
constexpr real propagate_average(const real error,
                                 const Args... rest) {
  return propagate_sum(error, rest...)
    / std::sqrt(1.0L + static_cast<real>(util::type::count(rest...)));
}
//----------------------------------------------------------------------------------------------

//__Propagate Error in Average of Independent Errors____________________________________________
template<class Range>
constexpr real propagate_average(const Range& range) {
  return propagate_sum(range)
    / std::sqrt(static_cast<real>(util::type::size(range)));
}
//----------------------------------------------------------------------------------------------

namespace detail { /////////////////////////////////////////////////////////////////////////////

//__Take Product of Event Arguments_____________________________________________________________
constexpr real even_argument_product(const real x,
                                     const real) {
  return x;
}
template<class ...Args>
constexpr real even_argument_product(const real x,
                                     const real,
                                     const Args ...args) {
  return x * even_argument_product(args...);
}
//----------------------------------------------------------------------------------------------

//__Fused Product of Ratio of Consecutive Arguments_____________________________________________
constexpr real ratio_fused_product(const real x,
                                   const real y) {
  return util::math::sum_squares(y / x);
}
template<class ...Args>
constexpr real ratio_fused_product(const real x,
                                   const real y,
                                   const Args ...args) {
  return std::fma(y / x, y / x, ratio_fused_product(args...));
}
//----------------------------------------------------------------------------------------------

} /* namespace detail */ ///////////////////////////////////////////////////////////////////////

//__Propagate Error in Product of Independent Errors____________________________________________
template<class ...Args>
constexpr real propagate_product(const real value,
                                 const real error,
                                 const Args... rest) {
  return util::math::abs(detail::even_argument_product(value, error, rest...))
    * std::sqrt(detail::ratio_fused_product(value, error, rest...));
}
//----------------------------------------------------------------------------------------------

//__Propagate Error in Product of Independent Errors____________________________________________
template<class Range>
constexpr real propagate_product(const Range& range) {
  real product = 1.0L, ratio_sum_squares = 0.0L;
  const auto end = std::cend(range);
  for (auto it = std::cbegin(range); it != end; it += 2) {
    product *= *it;
    const auto ratio = *(it+1) / *it;
    ratio_sum_squares = std::fma(ratio, ratio, ratio_sum_squares);
  }
  return util::math::abs(product) * std::sqrt(ratio_sum_squares);
}
//----------------------------------------------------------------------------------------------

//__Error of Uniform Random Variable____________________________________________________________
inline real uniform(const real width) {
  static const real inv_sqrt12 = 1.0L / std::sqrt(12.0L);
  return util::math::abs(width * inv_sqrt12);
}
//----------------------------------------------------------------------------------------------

//__Error of Uniform Random Variable____________________________________________________________
inline real uniform(const real a,
                    const real b) {
  return uniform(b - a);
}
//----------------------------------------------------------------------------------------------

} /* namespace error */ ////////////////////////////////////////////////////////////////////////

namespace type { ///////////////////////////////////////////////////////////////////////////////

//__Self-Propagating Uncertainty Type___________________________________________________________
template<class T>
struct uncertain {
  T value, error;

  uncertain(T v, T e) : value(v), error(e) {}

  static uncertain from_sum(T value1,
                            T error1,
                            T value2,
                            T error2);

  static uncertain from_difference(T value1,
                                   T error1,
                                   T value2,
                                   T error2);

  static uncertain from_product(T value1,
                                T error1,
                                T value2,
                                T error2);

  static uncertain from_quotient(T value1,
                                 T error1,
                                 T value2,
                                 T error2);

  static uncertain from_uniform(T v, T width);
  static uncertain from_uniform(T v, T a, T b);

  uncertain(const uncertain& other) = default;
  uncertain(uncertain&& other) = default;
  uncertain& operator=(const uncertain& other) = default;
  uncertain& operator=(uncertain&& other) = default;

  operator T() const { return value; }
};
//----------------------------------------------------------------------------------------------

//__Self-Propagating Real Number Alias__________________________________________________________
using uncertain_real = uncertain<real>;
//----------------------------------------------------------------------------------------------

//__Sum of Uncertainty Types____________________________________________________________________
template<class T>
constexpr uncertain<T> operator+=(uncertain<T>& left,
                                 const uncertain<T>& right) {
  left.error = stat::error::propagate_sum(left.error, right.error);
  left.value += right.value;
  return left;
}
template<class T>
constexpr uncertain<T> operator+(uncertain<T> left,
                                 const uncertain<T>& right) {
  return left += right;
}
//----------------------------------------------------------------------------------------------

//__Difference of Uncertainty Types_____________________________________________________________
template<class T>
constexpr uncertain<T> operator-=(uncertain<T>& left,
                                 const uncertain<T>& right) {
  left.error = stat::error::propagate_sum(left.error, right.error);
  left.value -= right.value;
  return left;
}
template<class T>
constexpr uncertain<T> operator-(uncertain<T> left,
                                 const uncertain<T>& right) {
  return left -= right;
}
//----------------------------------------------------------------------------------------------

//__Product of Uncertainty Types________________________________________________________________
template<class T>
constexpr uncertain<T> operator*=(uncertain<T>& left,
                                 const uncertain<T>& right) {
  left.error = stat::error::propagate_product(left.value, left.error,
                                              right.value, right.error);
  left.value *= right.value;
  return left;
}
template<class T>
constexpr uncertain<T> operator*(uncertain<T> left,
                                 const uncertain<T>& right) {
  return left *= right;
}
//----------------------------------------------------------------------------------------------

//__Quotient of Uncertainty Types_______________________________________________________________
template<class T>
constexpr uncertain<T> operator/=(uncertain<T>& left,
                                 const uncertain<T>& right) {
  left.error = stat::error::propagate_product(left.value, left.error,
                                              1.0L / right.value, right.error);
  left.value /= right.value;
  return left;
}
template<class T>
constexpr uncertain<T> operator/(uncertain<T> left,
                                 const uncertain<T>& right) {
  return left /= right;
}
//----------------------------------------------------------------------------------------------

//__Specialty Factory Constructors______________________________________________________________
template<class T>
uncertain<T> uncertain<T>::from_sum(T value1,
                                    T error1,
                                    T value2,
                                    T error2) {
  return uncertain<T>(value1, error1) + uncertain<T>(value2, error2);
}
template<class T>
uncertain<T> uncertain<T>::from_difference(T value1,
                                           T error1,
                                           T value2,
                                           T error2) {
  return uncertain<T>(value1, error1) - uncertain<T>(value2, error2);
}
template<class T>
uncertain<T> uncertain<T>::from_product(T value1,
                                        T error1,
                                        T value2,
                                        T error2) {
  return uncertain<T>(value1, error1) * uncertain<T>(value2, error2);
}
template<class T>
uncertain<T> uncertain<T>::from_quotient(T value1,
                                         T error1,
                                         T value2,
                                         T error2) {
  return uncertain<T>(value1, error1) / uncertain<T>(value2, error2);
}
template<class T>
uncertain<T> uncertain<T>::from_uniform(T value,
                                        T width) {
  return uncertain<T>(value, stat::error::uniform(width));
}
template<class T>
uncertain<T> uncertain<T>::from_uniform(T value,
                                        T a,
                                        T b) {
  return uncertain<T>(value, stat::error::uniform(a, b));
}
//----------------------------------------------------------------------------------------------

} /* namespace type */ /////////////////////////////////////////////////////////////////////////

namespace random { /////////////////////////////////////////////////////////////////////////////

//__Random Distribution Types___________________________________________________________________
enum class Distribution {
  UniformInt,
  UniformReal,
  Bernoulli,
  Binomial,
  NegativeBinomial,
  Geometric,
  Poisson,
  Exponential,
  Gamma,
  Weibull,
  ExtremeValue,
  Normal,
  LogNormal,
  ChiSquared,
  Cauchy,
  FisherF,
  StudentT,
  Discrete,
  PiecewiseConstant,
  PiecewiseLinear,
};
//----------------------------------------------------------------------------------------------

//__Random Distribution Parameter Base Type_____________________________________________________
template<Distribution D>
struct distribution_parameters { static constexpr auto type = D; };
//----------------------------------------------------------------------------------------------

//__Random Distribution Parameter Derived Type__________________________________________________
struct uniform_int       : public distribution_parameters<Distribution::UniformInt> {
  integer a, b;
  uniform_int(integer left, integer right)         : a(left), b(right) {}
};
struct uniform_real      : public distribution_parameters<Distribution::UniformReal> {
  real a, b;
  uniform_real(real left, real right)              : a(left), b(right) {}
};
struct bernoulli         : public distribution_parameters<Distribution::Bernoulli> {
  real p;
  bernoulli(real prob)                             : p(prob) {}
};
struct binomial          : public distribution_parameters<Distribution::Binomial> {
  integer t; real p;
  binomial(std::size_t trials, real prob)          : t(trials), p(prob) {}
};
struct negative_binomial : public distribution_parameters<Distribution::NegativeBinomial> {
  integer k; real p;
  negative_binomial(std::size_t trials, real prob) : k(trials), p(prob) {}
};
struct geometric         : public distribution_parameters<Distribution::Geometric> {
  real p;
  geometric(real prob)                             : p(prob) {}
};
struct poisson           : public distribution_parameters<Distribution::Poisson> {
  real mean;
  poisson(real m)                                  : mean(m) {}
};
struct exponential       : public distribution_parameters<Distribution::Exponential> {
  real lambda;
  exponential(real l)                              : lambda(l) {}
};
struct gamma             : public distribution_parameters<Distribution::Gamma> {
  real alpha, beta;
  gamma(real a, real b)                            : alpha(a), beta(b) {}
};
struct weibull           : public distribution_parameters<Distribution::Weibull> {
  real a, b;
  weibull(real pa, real pb)                        : a(pa), b(pb) {}
};
struct extreme_value     : public distribution_parameters<Distribution::ExtremeValue> {
  real a, b;
  extreme_value(real pa, real pb)                  : a(pa), b(pb) {}
};
struct normal            : public distribution_parameters<Distribution::Normal> {
  real mean, stddev;
  normal(real m, real s)                           : mean(m), stddev(s) {}
};
struct lognormal         : public distribution_parameters<Distribution::LogNormal> {
  real mean, stddev;
  lognormal(real m, real s)                        : mean(m), stddev(s) {}
};
struct chi_squared       : public distribution_parameters<Distribution::ChiSquared> {
  real n;
  chi_squared(real dof)                            : n(dof) {}
};
struct cauchy            : public distribution_parameters<Distribution::Cauchy> {
  real a, b;
  cauchy(real pa, real pb)                         : a(pa), b(pb) {}
};
struct fisher_f          : public distribution_parameters<Distribution::FisherF> {
  real m, n;
  fisher_f(real pm, real pn)                       : m(pm), n(pn) {}
};
struct student_t         : public distribution_parameters<Distribution::StudentT> {
  real n;
  student_t(real dof)                              : n(dof) {}
};
// TODO: finish
// struct discrete : public distribution_parameters<Distribution::Discrete> {
//   real_vector probabilities;
// };
struct piecewise_constant : public distribution_parameters<Distribution::PiecewiseConstant> {
  real_vector intervals, densities;
  piecewise_constant(real_vector i, real_vector d) : intervals(i), densities(d) {}
};
// TODO: finish
// struct piecewise_linear : public distribution_parameters<Distribution::PiecewiseLinear> {
//   real_vector intervals, densities;
// };
//----------------------------------------------------------------------------------------------

namespace detail { /////////////////////////////////////////////////////////////////////////////

//__Random Distribution Generator Type Union____________________________________________________
union distribution_union {
  distribution_union() {}
  std::uniform_int_distribution<integer>        uniform_int;
  std::uniform_real_distribution<real>          uniform_real;
  std::bernoulli_distribution                   bernoulli;
  std::binomial_distribution<integer>           binomial;
  std::negative_binomial_distribution<integer>  negative_binomial;
  std::geometric_distribution<integer>          geometric;
  std::poisson_distribution<integer>            poisson;
  std::exponential_distribution<real>           exponential;
  std::gamma_distribution<real>                 gamma;
  std::weibull_distribution<real>               weibull;
  std::extreme_value_distribution<real>         extreme_value;
  std::normal_distribution<real>                normal;
  std::lognormal_distribution<real>             lognormal;
  std::chi_squared_distribution<real>           chi_squared;
  std::cauchy_distribution<real>                cauchy;
  std::fisher_f_distribution<real>              fisher_f;
  std::student_t_distribution<real>             student_t;
  std::discrete_distribution<integer>           discrete;
  std::piecewise_constant_distribution<real>    piecewise_constant;
  std::piecewise_linear_distribution<real>      piecewise_linear;
  ~distribution_union() {}
};
//----------------------------------------------------------------------------------------------

//__Apply Function to Proper Distribution Type__________________________________________________
template<class UnaryFunction>
UnaryFunction distribution_apply(const Distribution type,
                                 distribution_union& dist,
                                 UnaryFunction f) {
  switch (type) {
    case Distribution::UniformInt:        f(dist.uniform_int);        break;
    case Distribution::UniformReal:       f(dist.uniform_real);       break;
    case Distribution::Bernoulli:         f(dist.bernoulli);          break;
    case Distribution::Binomial:          f(dist.binomial);           break;
    case Distribution::NegativeBinomial:  f(dist.negative_binomial);  break;
    case Distribution::Geometric:         f(dist.geometric);          break;
    case Distribution::Poisson:           f(dist.poisson);            break;
    case Distribution::Exponential:       f(dist.exponential);        break;
    case Distribution::Gamma:             f(dist.gamma);              break;
    case Distribution::Weibull:           f(dist.weibull);            break;
    case Distribution::ExtremeValue:      f(dist.extreme_value);      break;
    case Distribution::Normal:            f(dist.normal);             break;
    case Distribution::LogNormal:         f(dist.lognormal);          break;
    case Distribution::ChiSquared:        f(dist.chi_squared);        break;
    case Distribution::Cauchy:            f(dist.cauchy);             break;
    case Distribution::FisherF:           f(dist.fisher_f);           break;
    case Distribution::StudentT:          f(dist.student_t);          break;
    case Distribution::Discrete:          f(dist.discrete);           break;
    case Distribution::PiecewiseConstant: f(dist.piecewise_constant); break;
    case Distribution::PiecewiseLinear:   f(dist.piecewise_linear);   break;
    default:                                                          break;
  }
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Apply Function to Proper Distribution Type without Modification_____________________________
template<class UnaryFunction>
UnaryFunction distribution_apply(const Distribution type,
                                 const distribution_union& dist,
                                 UnaryFunction f) {
  switch (type) {
    case Distribution::UniformInt:        f(dist.uniform_int);        break;
    case Distribution::UniformReal:       f(dist.uniform_real);       break;
    case Distribution::Bernoulli:         f(dist.bernoulli);          break;
    case Distribution::Binomial:          f(dist.binomial);           break;
    case Distribution::NegativeBinomial:  f(dist.negative_binomial);  break;
    case Distribution::Geometric:         f(dist.geometric);          break;
    case Distribution::Poisson:           f(dist.poisson);            break;
    case Distribution::Exponential:       f(dist.exponential);        break;
    case Distribution::Gamma:             f(dist.gamma);              break;
    case Distribution::Weibull:           f(dist.weibull);            break;
    case Distribution::ExtremeValue:      f(dist.extreme_value);      break;
    case Distribution::Normal:            f(dist.normal);             break;
    case Distribution::LogNormal:         f(dist.lognormal);          break;
    case Distribution::ChiSquared:        f(dist.chi_squared);        break;
    case Distribution::Cauchy:            f(dist.cauchy);             break;
    case Distribution::FisherF:           f(dist.fisher_f);           break;
    case Distribution::StudentT:          f(dist.student_t);          break;
    case Distribution::Discrete:          f(dist.discrete);           break;
    case Distribution::PiecewiseConstant: f(dist.piecewise_constant); break;
    case Distribution::PiecewiseLinear:   f(dist.piecewise_linear);   break;
    default:                                                          break;
  }
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

template<class Dist, class ...Args>
void reconstruct_distribution(Dist& dist,
                              Args&& ...args) {
  new (&dist) Dist(std::forward<Args>(args)...);
}
//----------------------------------------------------------------------------------------------

inline void destruct_distribution(const Distribution type,
                                  distribution_union& dist) {
  distribution_apply(type, dist, [](auto& d) {
    using old_type = std::remove_reference_t<decltype(d)>;
    d.~old_type();
  });
}
//----------------------------------------------------------------------------------------------

template<class Dist, class ...Args>
void set_distribution(const Distribution old_type,
                      distribution_union& dist_union,
                      Dist& new_dist,
                      Args&& ...args) {
  destruct_distribution(old_type, dist_union);
  reconstruct_distribution(new_dist, std::forward<Args>(args)...);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameter Vector___________________________________________________________
template<class ParameterType>
void set_distribution_parameters(ParameterType&& parameters,
                                 distribution_union& distribution,
                                 const Distribution previous_type);
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Uniform Integer Distribution________________________________
template<>
inline void set_distribution_parameters<uniform_int>(uniform_int&& parameters,
                                                     distribution_union& distribution,
                                                     const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.uniform_int,
    parameters.a, parameters.b);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Uniform Real Distribution___________________________________
template<>
inline void set_distribution_parameters<uniform_real>(uniform_real&& parameters,
                                                      distribution_union& distribution,
                                                      const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.uniform_real,
    parameters.a, parameters.b);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Bernoulli Distribution______________________________________
template<>
inline void set_distribution_parameters<bernoulli>(bernoulli&& parameters,
                                                   distribution_union& distribution,
                                                   const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.bernoulli,
    parameters.p);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Binomial Distribution_______________________________________
template<>
inline void set_distribution_parameters<binomial>(binomial&& parameters,
                                                  distribution_union& distribution,
                                                  const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.binomial,
    parameters.t, parameters.p);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Negative Binomial Distribution______________________________
template<>
inline void set_distribution_parameters<negative_binomial>(negative_binomial&& parameters,
                                                           distribution_union& distribution,
                                                           const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.negative_binomial,
    parameters.k, parameters.p);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Geometric Distribution______________________________________
template<>
inline void set_distribution_parameters<geometric>(geometric&& parameters,
                                                   distribution_union& distribution,
                                                   const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.geometric,
    parameters.p);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Poisson Distribution________________________________________
template<>
inline void set_distribution_parameters<poisson>(poisson&& parameters,
                                                 distribution_union& distribution,
                                                 const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.poisson,
    parameters.mean);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Exponential Distribution____________________________________
template<>
inline void set_distribution_parameters<exponential>(exponential&& parameters,
                                                     distribution_union& distribution,
                                                     const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.exponential,
    parameters.lambda);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Gamma Distribution__________________________________________
template<>
inline void set_distribution_parameters<gamma>(gamma&& parameters,
                                               distribution_union& distribution,
                                               const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.gamma,
    parameters.alpha, parameters.beta);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Weibull Distribution________________________________________
template<>
inline void set_distribution_parameters<weibull>(weibull&& parameters,
                                                 distribution_union& distribution,
                                                 const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.weibull,
    parameters.a, parameters.b);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Extreme Value Distribution__________________________________
template<>
inline void set_distribution_parameters<extreme_value>(extreme_value&& parameters,
                                                       distribution_union& distribution,
                                                       const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.extreme_value,
    parameters.a, parameters.b);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Normal Distribution_________________________________________
template<>
inline void set_distribution_parameters<normal>(normal&& parameters,
                                                distribution_union& distribution,
                                                const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.normal,
    parameters.mean, parameters.stddev);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for LogNormal Distribution______________________________________
template<>
inline void set_distribution_parameters<lognormal>(lognormal&& parameters,
                                                   distribution_union& distribution,
                                                   const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.lognormal,
    parameters.mean, parameters.stddev);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Chi Squared Distribution____________________________________
template<>
inline void set_distribution_parameters<chi_squared>(chi_squared&& parameters,
                                                     distribution_union& distribution,
                                                     const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.chi_squared,
    parameters.n);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Cauchy Distribution_________________________________________
template<>
inline void set_distribution_parameters<cauchy>(cauchy&& parameters,
                                                distribution_union& distribution,
                                                const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.cauchy,
    parameters.a, parameters.b);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Fisher F Distribution_______________________________________
template<>
inline void set_distribution_parameters<fisher_f>(fisher_f&& parameters,
                                                  distribution_union& distribution,
                                                  const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.fisher_f,
    parameters.m, parameters.n);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Student T Distribution______________________________________
template<>
inline void set_distribution_parameters<student_t>(student_t&& parameters,
                                                   distribution_union& distribution,
                                                   const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.student_t,
    parameters.n);
}
//----------------------------------------------------------------------------------------------

//__Set Distribution Parameters for Student T Distribution______________________________________
template<>
inline void set_distribution_parameters<piecewise_constant>(piecewise_constant&& parameters,
                                                            distribution_union& distribution,
                                                            const Distribution previous_type) {
  set_distribution(previous_type, distribution, distribution.piecewise_constant,
    parameters.intervals.begin(), parameters.intervals.end(), parameters.densities.begin());
}
//----------------------------------------------------------------------------------------------

//__Build Distribution Parameter Type from Distribution Parameters______________________________
template<class ParameterType>
const ParameterType get_distribution_parameters(const distribution_union& distribution);
//----------------------------------------------------------------------------------------------

//__Build Uniform Integer Distribution Parameter Type from Distribution Parameters______________
template<>
inline const uniform_int get_distribution_parameters<uniform_int>(const distribution_union& distribution) {
  const auto& dist = distribution.uniform_int;
  return uniform_int(dist.a(), dist.b());
}
//----------------------------------------------------------------------------------------------

//__Build Uniform Real Distribution Parameter Type from Distribution Parameters_________________
template<>
inline const uniform_real get_distribution_parameters<uniform_real>(const distribution_union& distribution) {
  const auto& dist = distribution.uniform_real;
  return uniform_real(dist.a(), dist.b());
}
//----------------------------------------------------------------------------------------------

//__Build Normal Distribution Parameter Type from Distribution Parameters_______________________
template<>
inline const poisson get_distribution_parameters<poisson>(const distribution_union& distribution) {
  const auto& dist = distribution.poisson;
  return poisson(dist.mean());
}
//----------------------------------------------------------------------------------------------

// TODO: finish get_parameters

//__Build Normal Distribution Parameter Type from Distribution Parameters_______________________
template<>
inline const normal get_distribution_parameters<normal>(const distribution_union& distribution) {
  const auto& dist = distribution.normal;
  return normal(dist.mean(), dist.stddev());
}
//----------------------------------------------------------------------------------------------

// TODO: finish get_parameters

} /* namespace detail */ ///////////////////////////////////////////////////////////////////////

//__Random Number Generator_____________________________________________________________________
class generator {
public:
  template<class Type>
  generator(Type&& dist) {
    // FIXME: better seeding or remove
    seed(time(nullptr));
    distribution(std::forward<Type>(dist));
  }
  generator() = default;
  ~generator() = default;

  template<class Type>
  const Type distribution() const {
    return detail::get_distribution_parameters<Type>(_distribution);
  }

  template<class Type>
  void distribution(Type&& dist) {
    detail::set_distribution_parameters(std::forward<Type>(dist), _distribution, _type);
    _type = dist.type;
    reset();
  }

  real min() const;
  real max() const;

  real operator()();
  operator real() { return (*this)(); };
  real next() { return (*this)(); }

  void seed(const std::uint_least32_t seed);
  void seed(std::seed_seq& seq);
  void reset();
  void discard(const std::size_t z);

  bool operator==(const generator& other) const;
  bool operator!=(const generator& other) const { return !(*this == other); }

  template<class Container, class FillFunction>
  FillFunction fill_container(const std::size_t count,
                              Container& c,
                              FillFunction f) {
    for (std::size_t i{}; i < count; ++i)
      f(c, next());
    return std::move(f);
  }

  template<class Container, class FillFunction>
  FillFunction fill_container_until(const std::size_t count,
                                    Container& c,
                                    FillFunction f) {
    for (std::size_t i{}; i < count; ++i)
      if (!f(c, next())) break;
    return std::move(f);
  }

  template<class Container, class FillFunction>
  FillFunction fill_container_until(Container& c,
                                    FillFunction f) {
    while (f(c, next())) {}
    return std::move(f);
  }

private:
  std::mt19937 _engine;
  detail::distribution_union _distribution;
  Distribution _type;
};
//----------------------------------------------------------------------------------------------

} /* namespace random */ ///////////////////////////////////////////////////////////////////////

} /* namespace stat */ /////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__CORE__STAT_HH */

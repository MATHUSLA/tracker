#include "analysis.hh"

#include <deque>
#include <numeric>

#include "ROOT/TMinuit.h"

#include "geometry.hh"
#include "util/combinatorics.hh"

#include "util/io.hh" // TODO: REMOVE

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Collapse Points by R4 Interval______________________________________________________________
event_points collapse(const event_points& event,
                      const r4_point& ds) {
  if (event.empty())
    return event;

  event_points out;

  const auto& sorted_event = t_copy_sort(event);
  const auto&& size = sorted_event.size();

  using size_type = event_points::size_type;

  size_type index = 0;
  std::deque<size_type> marked_indicies{};
  while (index < size) {
    size_type collected = 1, missed_index = 0;
    const auto& point = sorted_event[index];
    const auto&& time_interval = point.t + ds.t;
    auto sum_t = point.t,
         sum_x = point.x,
         sum_y = point.y,
         sum_z = point.z;

    auto skipped = false;
    while (++index < size) {

      while (!marked_indicies.empty() && index++ == marked_indicies.front())
        marked_indicies.pop_front();

      const auto& next = sorted_event[index];
      if (next.t > time_interval)
        break;

      if (within_dr(point, next, ds)) {
        ++collected;
        sum_t += next.t;
        sum_x += next.x;
        sum_y += next.y;
        sum_z += next.z;
        if (skipped)
          marked_indicies.push_back(index);
      } else if (!skipped) {
        skipped = true;
        missed_index = index;
      }
    }

    if (skipped)
      index = missed_index;

    out.push_back({sum_t / collected, sum_x / collected, sum_y / collected, sum_z / collected});
  }

  return out;
}
//----------------------------------------------------------------------------------------------

//__Partition Points by Coordinate______________________________________________________________
event_partition partition(const event_points& points,
                          const real interval,
                          const Coordinate coordinate) {
  event_partition out{{}, coordinate};
  if (points.empty())
    return out;

  auto& parts = out.parts;

  const auto& sorted_points = coordinate_copy_sort(points, coordinate);
  const auto&& size = sorted_points.size();

  event_points::size_type count = 0;
  auto point_iter = sorted_points.cbegin();
  while (count < size) {
    const auto& point = *point_iter;
    event_points current_layer{point};
    ++count;

    while (count < size) {
      const auto& next = *(++point_iter);
      if ((coordinate == Coordinate::T && (next.t > point.t + interval)) ||
          (coordinate == Coordinate::X && (next.x > point.x + interval)) ||
          (coordinate == Coordinate::Y && (next.y > point.y + interval)) ||
          (coordinate == Coordinate::Z && (next.z > point.z + interval))) {
        break;
      }
      current_layer.push_back(next);
      ++count;
    }

    parts.push_back(t_sort(current_layer));
  }

  return out;
}
//----------------------------------------------------------------------------------------------

//__Seeding Algorithm___________________________________________________________________________
event_vector seed(const size_t n,
                  const event_points& event,
                  const r4_point& collapse_ds,
                  const real layer_dz,
                  const real line_dr) {
  if (n <= 2)
    return {};

  const auto& points = collapse(event, collapse_ds);
  const auto&& size = points.size();

  if (size <= n)
    return { points };

  event_vector out;
  out.reserve(std::pow(size, n) / std::pow(n/2.718, n));  // work on this limit

  const auto& layers = partition(points, layer_dz).parts;
  const auto&& layer_count = layers.size();
  if (layer_count < n)  // FIXME: unsure what to do here
    return {};

  util::combinatorics::bit_vector_sequence layer_sequence;
  for (const auto& layer : layers)
    layer_sequence.emplace_back(1, layer.size());

  util::combinatorics::order2_permutations(n, layer_sequence, [&](const auto& chooser) {
    event_points tuple;
    tuple.reserve(n);

    for (size_t i = 0; i < layer_count; ++i) {
      if (chooser[i]) {
        const auto& layer = layers[i];
        const auto& bits = layer_sequence[i];
        size_t j = 0;
        std::copy_if(layer.cbegin(), layer.cend(), std::back_inserter(tuple),
          [&](const auto& _) { return bits[j++]; });
      }
    }

    t_sort(tuple);

    const auto& line_begin = tuple.front();
    const auto& line_end = tuple.back();
    const auto&& max = std::accumulate(++tuple.cbegin(), --tuple.cend(), line_dr,
      [&](const auto& max, const auto& point) {
        return std::max(max, type::point_line_distance(point, line_begin, line_end)); });

    if (max <= line_dr)
      out.push_back(tuple);
  });

  return out;
}
//----------------------------------------------------------------------------------------------

//__Seed Search_________________________________________________________________________________
event_vector seeds_with(const r4_point& point,
                        const event_vector& seeds) {
  event_vector out;
  out.reserve(seeds.size());
  std::copy_if(seeds.cbegin(), seeds.cend(), std::back_inserter(out), [&](const auto& seed) {
    const auto& end = seed.cend();
    return end != std::find(seed.cbegin(), end, point);
  });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Seed Search Front___________________________________________________________________________
event_vector seeds_starting_with(const r4_point& point,
                                 const event_vector& seeds) {
  event_vector out;
  out.reserve(seeds.size());
  std::copy_if(seeds.cbegin(), seeds.cend(), std::back_inserter(out),
    [&](const auto& seed) { return point == seed.front(); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Seed Search Back____________________________________________________________________________
event_vector seeds_ending_with(const r4_point& point,
                               const event_vector& seeds) {
  event_vector out;
  out.reserve(seeds.size());
  std::copy_if(seeds.cbegin(), seeds.cend(), std::back_inserter(out),
    [&](const auto& seed) { return point == seed.back(); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Seed Merge__________________________________________________________________________________
event_vector merge(const event_vector& seeds) {
  event_vector out{};
  out.reserve(seeds.size());

  // TODO: implement

  for (const auto& seed : seeds) {
    for (const auto& point : seed) {
      const auto& ss = seeds_starting_with(point, seeds);
      std::for_each(ss.cbegin(), ss.cend(),
        [](const auto& seed) { util::io::print_range(seed, " :: ") << "\n"; });
    }
  }

  return out;
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////
#define FP_FAST_FMA
#define FP_FAST_FMAF
#define FP_FAST_FMAL
//__Gaussian Negative Log Likelihood Calculation________________________________________________
thread_local event_points&& _event = {};
void _gaussian_nll(Int_t&, Double_t*, Double_t& out, Double_t* parameters, Int_t) {
  const auto& t0 = parameters[0];
  const auto& x0 = parameters[1];
  const auto& y0 = parameters[2];
  const auto& z0 = parameters[3];
  const auto& vx = parameters[4];
  const auto& vy = parameters[5];
  const auto& vz = parameters[6];

  out = 0.5 * std::accumulate(_event.cbegin(), _event.cend(), 0, [&](const auto& sum, const auto& point) {
    const auto& limits = geometry::limits_of(geometry::volume(point));
    const auto& center = limits.center;
    const auto& min = limits.min;
    const auto& max = limits.max;
    const auto&& dz = (center.z - z0) / vz;
    const auto&& t_res = (dz + t0) / (2 /* nanoseconds */);
    const auto&& x_res = (std::fma(vx, dz, x0) - center.x) / (max.x - min.x);
    const auto&& y_res = (std::fma(vy, dz, y0) - center.y) / (max.y - min.y);
    return sum + (t_res*t_res + 12*x_res*x_res + 12*y_res*y_res);
  });
}
//----------------------------------------------------------------------------------------------
} /* annonymous namespace */ ///////////////////////////////////////////////////////////////////

//__MINUIT Gaussian Fitter______________________________________________________________________
void fit_event(const event_points& event,
               fit_parameter_vector& parameters,
               const fit_settings& settings) {
  TMinuit minuit;
  minuit.SetPrintLevel(-1);
  minuit.SetErrorDef(settings.error_def);

  const auto&& size = parameters.size();
  for (size_t i = 0; i < size; ++i) {
    const auto& param = parameters[i];
    minuit.DefineParameter(i,
      param.name.c_str(), param.value, param.error, param.min, param.max);
  }

  _event = event;
  minuit.SetFCN(_gaussian_nll);

  Int_t error_flag;
  auto command_parameters = settings.command_parameters;
  minuit.mnexcm(
    settings.command_name.c_str(),
    command_parameters.data(),
    command_parameters.size(),
    error_flag);

  for (size_t i = 0; i < size; ++i) {
    Double_t value, error;
    minuit.GetParameter(i, value, error);
    auto& param = parameters[i];
    param.value = value;
    param.error = error;
  }
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

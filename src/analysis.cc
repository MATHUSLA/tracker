#include "analysis.hh"

#include <deque>

#include "ROOT/TMinuit.h"

#include "util.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Collapse Points by R4 Interval______________________________________________________________
event_points collapse(const event_points& event,
                      const r4_point& ds) {
  if (event.empty())
    return event;

  event_points out;

  const auto& sorted_events = t_copy_sort(event);
  const auto&& event_size = sorted_events.size();

  event_points::size_type index = 0;
  std::deque<decltype(index)> marked_indicies{};
  while (index < event_size) {
    event_points::size_type collected = 1, missed_index = 0;
    const auto& point = sorted_events[index];
    const auto&& time_interval = point.t + ds.t;
    auto sum_t = point.t,
         sum_x = point.x,
         sum_y = point.y,
         sum_z = point.z;

    auto skipped = false;
    while (++index < event_size) {
      while (!marked_indicies.empty() && index == marked_indicies.front()) {
        ++index;
        marked_indicies.pop_front();
      }

      const auto& next = sorted_events[index];
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

    t_sort(current_layer);
    parts.push_back(current_layer);
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

  if (layer_count < n)  // FIXME: what to do about this? maybe recurse to a lower level? (seed(n-1))
    return {};

  bit_vector_sequence layer_sequence;
  for (const auto& layer : layers)
    layer_sequence.push_back(bit_vector(1, layer.size()));

  combinatorics::order2_permutations(n, layer_sequence, [&](const auto& chooser) {
    event_points tuple;
    tuple.reserve(n);

    for (size_t i = 0; i < chooser.size(); ++i) {
      if (chooser[i]) {
        const auto& layer_bits = layer_sequence[i];
        const auto&& layer_size = layer_bits.size();
        for (size_t j = 0; j < layer_size; ++j)
          if (layer_bits[j]) tuple.push_back(layers[i][j]);
      }
    }

    t_sort(tuple);

    if (n > 2) {
      const auto& line_begin = tuple.front();
      const auto& line_end = tuple.back();
      real max = -1;
      for (size_t i = 1; i < n - 1; ++i)
        max = std::max(max, type::point_line_distance(tuple[i], line_begin, line_end));
      std::cout << max << ": ";
      if (max == -1 || max > line_dr) { std::cout << "\n"; return; }
    }

    for (const auto& t : tuple)
      std::cout << t << " ";
    std::cout << "\n";

    out.push_back(tuple);
  });

  return out;
}
//----------------------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Gaussian Negative Log Likelihood Calculation________________________________________________
static thread_local event_vector&& _events = {};
static void _gaussian_nll(Int_t&, Double_t*, Double_t& value, Double_t* parameters, Int_t) {
  value = 0;
  parameters = nullptr;
}
//----------------------------------------------------------------------------------------------
} /* annonymous namespace */ ///////////////////////////////////////////////////////////////////

//__MINUIT Gaussian Fitter______________________________________________________________________
void fit_events(const event_vector& events,
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

  _events = events;
  minuit.SetFCN(_gaussian_nll);

  int error_flag;
  auto command_parameters = settings.command_parameters;
  minuit.mnexcm(
    settings.command_name.c_str(),
    command_parameters.data(),
    command_parameters.size(),
    error_flag);

  for (size_t i = 0; i < size; ++i) {
    real value, error;
    minuit.GetParameter(i, value, error);
    auto& param = parameters[i];
    param.value = value;
    param.error = error;
  }
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

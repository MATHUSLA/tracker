#ifndef TRACKER__ANALYSIS_HH
#define TRACKER__ANALYSIS_HH
#pragma once

#include "event.hh"
#include "util.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

using namespace type;

//__Collapse Points by R4 Interval______________________________________________________________
event_points collapse(const event_points& event,
                      const r4_point& ds);
//----------------------------------------------------------------------------------------------

//__Partition Points by Coordinate______________________________________________________________
event_partition partition(const event_points& points,
                          const real interval,
                          const Coordinate coordinate=Coordinate::Z);
//----------------------------------------------------------------------------------------------

//__SeedN Algorithm_____________________________________________________________________________
template<std::size_t N>
event_tuple_vector<N> seed(const event_points& event,
                           const r4_point& collapse_ds,
                           const real layer_dz) {
  const auto& points = collapse(event, collapse_ds);
  const auto&& size = points.size();

  if (size <= N)
    return {to_array<N>(points)};

  event_tuple_vector<N> out;
  out.reserve(std::pow(size, N) / std::pow(N/2.718, N));  // work on this limit

  const auto& layers = partition(points, layer_dz).parts;
  const auto&& layer_count = layers.size();

  if (layer_count < N)  // FIXME: what to do about this? maybe recurse to a lower level? (seed<N-1>)
    return {};

  bit_vector_sequence layer_sequence;
  for (const auto& layer : layers)
    layer_sequence.push_back(bit_vector(1, layer.size()));

  combinatorics::order2_permutations(N, layer_sequence, [&](const auto& chooser) {
    event_tuple<N> tuple;
    for (size_t i = 0, index = 0; i < chooser.size(); ++i) {
      if (chooser[i]) {
        const auto& layer_bits = layer_sequence[i];
        const auto&& layer_size = layer_bits.size();
        for (size_t j = 0; j < layer_size; ++j)
          if (layer_bits[j]) tuple[index++] = layers[i][j];
      }
    }

    /*
    for (const auto& t : tuple)
      std::cout << t << " ";
    std::cout << "\n";
    */

    out.push_back(tuple);
  });

  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS_HH */

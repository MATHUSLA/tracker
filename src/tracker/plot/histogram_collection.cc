/*
 * src/tracker/plot/histogram_collection.cc
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

#include <tracker/plot.hh>

namespace MATHUSLA { namespace TRACKER {

namespace plot { ///////////////////////////////////////////////////////////////////////////////

//__Collection Prefix Constructor_______________________________________________________________
histogram_collection::histogram_collection(const std::string& prefix) : _prefix(prefix) {}
//----------------------------------------------------------------------------------------------

//__Collection Constructor______________________________________________________________________
histogram_collection::histogram_collection(std::initializer_list<histogram>&& histograms) {
  for (auto&& hist : histograms)
    _histograms.insert(std::make_pair(hist.name(), hist));
}
//----------------------------------------------------------------------------------------------

//__Collection Constructor with Prefix__________________________________________________________
histogram_collection::histogram_collection(const std::string& prefix,
                                           std::initializer_list<histogram>&& histograms) : _prefix(prefix) {
  for (auto&& hist : histograms) {
    const auto full_name = _prefix + hist.name();
    _histograms.insert(std::make_pair(full_name, hist));
    _histograms[full_name].name(full_name);
  }
}
//----------------------------------------------------------------------------------------------

//__Emplacement Construction____________________________________________________________________
histogram& histogram_collection::emplace(const histogram& hist) {
  const auto full_name = _prefix + hist.name();
  _histograms.emplace(std::piecewise_construct, std::forward_as_tuple(full_name), std::forward_as_tuple(hist));
  auto& reference = _histograms[full_name];
  reference.name(full_name);
  return reference;
}
//----------------------------------------------------------------------------------------------

//__Emplacement Construction____________________________________________________________________
histogram& histogram_collection::emplace(histogram&& hist) {
  const auto full_name = _prefix + hist.name();
  _histograms.emplace(std::piecewise_construct, std::forward_as_tuple(full_name), std::forward_as_tuple(hist));
  auto& reference = _histograms[full_name];
  reference.name(full_name);
  return reference;
}
//----------------------------------------------------------------------------------------------

//__Emplacement Construction____________________________________________________________________
histogram& histogram_collection::emplace(const std::string& name,
                                         const size_t bins,
                                         const real min,
                                         const real max) {
  return _histograms.emplace(std::piecewise_construct, std::forward_as_tuple(_prefix + name),
                             std::forward_as_tuple(_prefix + name, bins, min, max)).first->second;
}
//----------------------------------------------------------------------------------------------

//__Emplacement Construction____________________________________________________________________
histogram& histogram_collection::emplace(const std::string& name,
                                         const std::string& title,
                                         const size_t bins,
                                         const real min,
                                         const real max) {
  return _histograms.emplace(std::piecewise_construct, std::forward_as_tuple(_prefix + name),
                             std::forward_as_tuple(_prefix + name, title, bins, min, max)).first->second;
}
//----------------------------------------------------------------------------------------------

//__Emplacement Construction____________________________________________________________________
histogram& histogram_collection::emplace(const std::string& name,
                                         const std::string& title,
                                         const std::string& x_title,
                                         const std::string& y_title,
                                         const size_t bins,
                                         const real min,
                                         const real max) {
  return _histograms.emplace(std::piecewise_construct, std::forward_as_tuple(_prefix + name),
                             std::forward_as_tuple(_prefix + name, title, x_title, y_title, bins, min, max)).first->second;
}
//----------------------------------------------------------------------------------------------

//__Indexing Operator Overload__________________________________________________________________
histogram& histogram_collection::operator[](const std::string& name) {
  const auto first_search = _histograms.find(name);
  if (first_search == _histograms.cend()) {
    return _histograms.find(_prefix + name)->second;
  }
  return first_search->second;
}
//----------------------------------------------------------------------------------------------

//__Load Histogram From File into Histogram Collection__________________________________________
histogram& histogram_collection::load(const std::string& path,
                                      const std::string& name) {
  auto& reference = _histograms.insert(std::make_pair(_prefix + name, histogram::load(path, name))).first->second;
  reference.name(_prefix + reference.name());
  return reference;
}
//----------------------------------------------------------------------------------------------

//__Load Histogram From File into Histogram Collection__________________________________________
histogram& histogram_collection::load_with_prefix(const std::string& path,
                                                  const std::string& name) {
  auto& reference = _histograms.insert(std::make_pair(_prefix + name, histogram::load(path, _prefix + name))).first->second;
  reference.name(_prefix + reference.name());
  return reference;
}
//----------------------------------------------------------------------------------------------

//__Draw all Histograms in Collection___________________________________________________________
void histogram_collection::draw_all() {
  std::for_each(_histograms.begin(), _histograms.end(), [](auto& elem) { elem.second.draw(); });
}
//----------------------------------------------------------------------------------------------

//__Clear all Histograms in Collection__________________________________________________________
void histogram_collection::clear_all() {
  std::for_each(_histograms.begin(), _histograms.end(), [](auto& elem) { elem.second.clear(); });
}
//----------------------------------------------------------------------------------------------

//__Save all Histograms in Collection___________________________________________________________
bool histogram_collection::save_all(const std::string& path) const {
  return std::accumulate(_histograms.cbegin(), _histograms.cend(), true,
    [&](const auto saved, const auto& elem) { return saved && elem.second.save(path); });
}
//----------------------------------------------------------------------------------------------

} /* namespace plot */  ////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */
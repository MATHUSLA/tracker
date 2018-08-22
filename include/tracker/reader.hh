/*
 * include/tracker/reader.hh
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

#ifndef TRACKER__READER_HH
#define TRACKER__READER_HH
#pragma once

#include <tracker/analysis/monte_carlo.hh>
#include <tracker/geometry.hh>
#include <tracker/script.hh>

namespace MATHUSLA { namespace TRACKER {

namespace reader { /////////////////////////////////////////////////////////////////////////////

using namespace type;

//__Detector Map Import_________________________________________________________________________
const geometry::detector_map import_detector_map(const script::path_type& path);
//----------------------------------------------------------------------------------------------

//__Detector Time Resolution Map Import_________________________________________________________
const geometry::time_resolution_map import_time_resolution_map(const script::path_type& path);
//----------------------------------------------------------------------------------------------

namespace root { ///////////////////////////////////////////////////////////////////////////////

//__Search and Collect File Paths_______________________________________________________________
const script::path_vector search_directory(const script::path_type& path,
                                           const std::string& ext="root");
//----------------------------------------------------------------------------------------------

//__Search and Collect File Paths_______________________________________________________________
const std::vector<script::path_vector> search_directories(const script::path_vector& paths,
                                                          const std::string& ext="root");
//----------------------------------------------------------------------------------------------

//__Transpose Search Directories________________________________________________________________
const std::vector<script::path_vector> transpose_search_directories(const script::path_vector& paths,
                                                                    const std::string& ext="root");
//----------------------------------------------------------------------------------------------

//__Import Mode Type____________________________________________________________________________
enum class ImportMode { Detector, Widths };
//----------------------------------------------------------------------------------------------

//__ROOT Event Import___________________________________________________________________________
const analysis::event_vector import_events(const script::path_type& path,
                                           const std::string& t_key,
                                           const std::string& x_key,
                                           const std::string& y_key,
                                           const std::string& z_key);
const analysis::event_vector import_events(const script::path_type& path,
                                           const std::string& t_key,
                                           const std::string& detector_key,
                                           const geometry::detector_map& map);
const analysis::event_vector import_events(const script::path_type& path,
                                           const script::tracking_options& options,
                                           const geometry::detector_map& map);
const analysis::event_vector import_events(const script::path_type& path,
                                           const script::tracking_options& options,
                                           const ImportMode mode);
//----------------------------------------------------------------------------------------------

//__ROOT Event Parallel Import__________________________________________________________________
const analysis::event_vector import_events(const script::path_vector& paths,
                                           const real_vector& timing_offsets,
                                           const std::string& t_key,
                                           const std::string& x_key,
                                           const std::string& y_key,
                                           const std::string& z_key);
const analysis::event_vector import_events(const script::path_vector& paths,
                                           const real_vector& timing_offsets,
                                           const std::string& t_key,
                                           const std::string& detector_key,
                                           const geometry::detector_map& map);
const analysis::event_vector import_events(const script::path_vector& paths,
                                           const script::tracking_options& options,
                                           const geometry::detector_map& map);
const analysis::event_vector import_events(const script::path_vector& paths,
                                           const script::tracking_options& options,
                                           const ImportMode mode);
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Import______________________________________________________________________
const analysis::full_event_vector import_full_events(const script::path_type& path,
                                                     const std::string& t_key,
                                                     const std::string& x_key,
                                                     const std::string& y_key,
                                                     const std::string& z_key,
                                                     const std::string& dt_key,
                                                     const std::string& dx_key,
                                                     const std::string& dy_key,
                                                     const std::string& dz_key);
const analysis::full_event_vector import_full_events(const script::path_type& path,
                                                     const std::string& t_key,
                                                     const std::string& dt_key,
                                                     const std::string& detector_key,
                                                     const geometry::detector_map& map);
const analysis::full_event_vector import_full_events(const script::path_type& path,
                                                     const script::tracking_options& options,
                                                     const geometry::detector_map& map);
const analysis::full_event_vector import_full_events(const script::path_type& path,
                                                     const script::tracking_options& options,
                                                     const ImportMode mode);
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Parallel Import_____________________________________________________________
const analysis::full_event_vector import_full_events(const script::path_vector& paths,
                                                     const real_vector& timing_offsets,
                                                     const std::string& t_key,
                                                     const std::string& x_key,
                                                     const std::string& y_key,
                                                     const std::string& z_key,
                                                     const std::string& dt_key,
                                                     const std::string& dx_key,
                                                     const std::string& dy_key,
                                                     const std::string& dz_key);
const analysis::full_event_vector import_full_events(const script::path_vector& paths,
                                                     const real_vector& timing_offsets,
                                                     const std::string& t_key,
                                                     const std::string& dt_key,
                                                     const std::string& detector_key,
                                                     const geometry::detector_map& map);
const analysis::full_event_vector import_full_events(const script::path_vector& paths,
                                                     const script::tracking_options& options,
                                                     const geometry::detector_map& map);
const analysis::full_event_vector import_full_events(const script::path_vector& paths,
                                                     const script::tracking_options& options,
                                                     const ImportMode mode);
//----------------------------------------------------------------------------------------------

//__ROOT Event Import with Monte-Carlo Tracks___________________________________________________
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_type& path,
                                                               const std::string& track_key,
                                                               const std::string& t_key,
                                                               const std::string& x_key,
                                                               const std::string& y_key,
                                                               const std::string& z_key);
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_type& path,
                                                               const std::string& track_key,
                                                               const std::string& t_key,
                                                               const std::string& x_key,
                                                               const std::string& y_key,
                                                               const std::string& z_key,
                                                               const std::string& detector_key,
                                                               const geometry::detector_map& map);
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_type& path,
                                                               const script::tracking_options& options,
                                                               const geometry::detector_map& map);
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_type& path,
                                                               const script::tracking_options& options,
                                                               const ImportMode mode);
//----------------------------------------------------------------------------------------------

//__ROOT Event Parallel Import with Monte-Carlo Tracks__________________________________________
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_vector& paths,
                                                               const real_vector& timing_offsets,
                                                               const std::string& track_key,
                                                               const std::string& t_key,
                                                               const std::string& x_key,
                                                               const std::string& y_key,
                                                               const std::string& z_key);
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_vector& paths,
                                                               const real_vector& timing_offsets,
                                                               const std::string& track_key,
                                                               const std::string& t_key,
                                                               const std::string& x_key,
                                                               const std::string& y_key,
                                                               const std::string& z_key,
                                                               const std::string& detector_key,
                                                               const geometry::detector_map& map);
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_vector& paths,
                                                               const script::tracking_options& options,
                                                               const geometry::detector_map& map);
const analysis::mc::event_vector_bundle import_event_mc_bundle(const script::path_vector& paths,
                                                               const script::tracking_options& options,
                                                               const ImportMode mode);
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Import with Monte-Carlo Tracks______________________________________________
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_type& path,
                                                                         const std::string& track_key,
                                                                         const std::string& t_key,
                                                                         const std::string& x_key,
                                                                         const std::string& y_key,
                                                                         const std::string& z_key,
                                                                         const std::string& dt_key,
                                                                         const std::string& dx_key,
                                                                         const std::string& dy_key,
                                                                         const std::string& dz_key);
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_type& path,
                                                                         const std::string& track_key,
                                                                         const std::string& t_key,
                                                                         const std::string& x_key,
                                                                         const std::string& y_key,
                                                                         const std::string& z_key,
                                                                         const std::string& dt_key,
                                                                         const std::string& detector_key,
                                                                         const geometry::detector_map& map);
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_type& path,
                                                                         const script::tracking_options& options,
                                                                         const geometry::detector_map& map);
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_type& path,
                                                                         const script::tracking_options& options,
                                                                         const ImportMode mode);
//----------------------------------------------------------------------------------------------

//__ROOT Full Event Parallel Import with Monte-Carlo Tracks_____________________________________
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_vector& paths,
                                                                         const real_vector& timing_offsets,
                                                                         const std::string& track_key,
                                                                         const std::string& t_key,
                                                                         const std::string& x_key,
                                                                         const std::string& y_key,
                                                                         const std::string& z_key,
                                                                         const std::string& dt_key,
                                                                         const std::string& dx_key,
                                                                         const std::string& dy_key,
                                                                         const std::string& dz_key);
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_vector& paths,
                                                                         const real_vector& timing_offsets,
                                                                         const std::string& track_key,
                                                                         const std::string& t_key,
                                                                         const std::string& x_key,
                                                                         const std::string& y_key,
                                                                         const std::string& z_key,
                                                                         const std::string& dt_key,
                                                                         const std::string& detector_key,
                                                                         const geometry::detector_map& map);
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_vector& paths,
                                                                         const script::tracking_options& options,
                                                                         const geometry::detector_map& map);
const analysis::mc::full_event_vector_bundle import_full_event_mc_bundle(const script::path_vector& paths,
                                                                         const script::tracking_options& options,
                                                                         const ImportMode mode);
//----------------------------------------------------------------------------------------------

//__Merge Save Files____________________________________________________________________________
void merge_save(const script::path_type& home_path,
                const script::path_vector& paths,
                const std::vector<std::string>& prefixes={});
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

} /* namespace reader */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__READER_HH */

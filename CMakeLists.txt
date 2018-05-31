#
# CMakeLists.txt
#
# Copyright 2018 Brandon Gomes
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(MATHUSLA_TRACKER VERSION 0.0.1 DESCRIPTION "MATHUSLA Particle Tracker")

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

find_package(Geant4 REQUIRED gdml)
find_package(ROOT   REQUIRED Minuit)

add_executable(tracker
    src/tracker.cc
    src/analysis.cc
    src/geometry.cc
    src/plot.cc
    src/reader.cc
    src/util/command_line_parser.cc
)

# TODO: make into library
#add_library(tracker SHARED
#    src/analysis.cc
#    src/geometry.cc
#    src/plot.cc
#    src/reader.cc
#    src/util/command_line_parser.cc
#)

# set_target_properties(tracker PROPERTIES VERSION ${PROJECT_VERSION})

target_link_libraries(tracker
    ${Geant4_LIBRARIES}
    ${ROOT_LIBRARIES})

target_compile_features(tracker PUBLIC cxx_std_14)

target_include_directories(tracker PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

install(TARGETS tracker DESTINATION bin)

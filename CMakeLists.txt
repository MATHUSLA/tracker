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

find_package(Geant4 REQUIRED COMPONENTS gdml)
find_package(ROOT   REQUIRED COMPONENTS MathCore Minuit)

add_library(tracker SHARED
    src/tracker/analysis/event.cc
    src/tracker/analysis/join.cc
    src/tracker/analysis/monte_carlo.cc
    src/tracker/analysis/track.cc
    src/tracker/analysis/tree.cc
    src/tracker/analysis/vertex.cc

    src/tracker/geometry.cc

    src/tracker/plot/canvas.cc
    src/tracker/plot/histogram.cc
    src/tracker/plot/histogram_collection.cc
    src/tracker/plot/plot.cc
    src/tracker/plot/value_tag.cc

    src/tracker/reader.cc
    src/tracker/script.cc

    src/tracker/core/stat.cc

    src/tracker/util/command_line_parser.cc
)

target_compile_features(tracker PUBLIC cxx_std_14)
target_compile_options(tracker PUBLIC -Wno-unknown-warning-option
                                      -Wno-missing-braces
                                      -Wunused-function
                                      -Wno-unused-value
                                      -Wno-unused-parameter
                                      -Wno-unused-but-set-parameter
                                      -Wno-switch)

set_target_properties(tracker PROPERTIES VERSION ${PROJECT_VERSION})

target_link_libraries(tracker PUBLIC
    ${Geant4_LIBRARIES}
    ${ROOT_LIBRARIES})

target_include_directories(tracker PUBLIC
    ${Geant4_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS})

target_include_directories(tracker PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

install(TARGETS tracker DESTINATION lib)

add_subdirectory(demo)

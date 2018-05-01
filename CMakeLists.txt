cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(MATHUSLA_TRACKER)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

find_package(Geant4 REQUIRED gdml ui_all vis_all)
find_package(ROOT   REQUIRED)

include(${Geant4_USE_FILE})

add_executable(tracker
    src/tracker.cc

    src/analysis.cc
    src/geometry.cc
    src/root_helper.cc

    src/util/CommandLineParser.cc
)

target_link_libraries(tracker
    ${Geant4_LIBRARIES}
    ${ROOT_LIBRARIES})

target_compile_features(tracker PUBLIC cxx_std_14)

target_include_directories(tracker PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

install(TARGETS tracker DESTINATION bin)

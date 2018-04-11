cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(MATHUSLA_TRACKER)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

find_package(Geant4 REQUIRED gdml ui_all vis_all)
find_package(ROOT   REQUIRED)

include(${Geant4_USE_FILE})

add_executable(tracker
    src/tracker.cc

    src/geometry/Construction.cc
    src/geometry/Reader.cc

    src/util/CommandLineParser.cc
)

target_link_libraries(tracker
    ${Geant4_LIBRARIES}
    ${ROOT_LIBRARIES})

target_include_directories(tracker PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

install(TARGETS tracker DESTINATION bin)

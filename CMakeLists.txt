cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(MATHUSLA_TRACKER)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

add_executable(tracker
    src/tracker.cc
)

# target_link_libraries(tracker)

target_include_directories(tracker PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_SOURCE_DIR}/include>)

install(TARGETS tracker DESTINATION bin)

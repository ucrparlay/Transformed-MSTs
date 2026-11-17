cmake_minimum_required(VERSION 3.18)

project(
  TMST-code
  VERSION 0.1.0
  DESCRIPTION "The family of balanced T-MSTs for incremental minimum spanning trees"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

enable_testing()

add_subdirectory(bench)
add_subdirectory(test)

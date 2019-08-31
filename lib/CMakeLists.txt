cmake_minimum_required(VERSION 3.1...3.14)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

include(FetchContent)

project(libs LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Formatting library from fmtlib for python-like i/o syntax.
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git)
FetchContent_GetProperties(fmt)
if(NOT fmt_POPULATED)
  FetchContent_Populate(fmt)
  add_subdirectory(
    ${fmt_SOURCE_DIR}
    ${fmt_BINARY_DIR}
    EXCLUDE_FROM_ALL)
endif()

add_subdirectory(am)
add_subdirectory(mcutils)
add_subdirectory(basis)

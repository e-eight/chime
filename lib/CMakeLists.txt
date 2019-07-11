cmake_minimum_required(VERSION 3.1...3.14)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(libs LANGUAGES CXX)

add_subdirectory(fmt)
add_subdirectory(am)
add_subdirectory(mcutils)
add_subdirectory(basis)

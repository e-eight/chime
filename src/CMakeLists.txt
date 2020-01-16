project (chime_src LANGUAGES CXX)

set (TARGET chime.exe)
set (TEST test.exe)

set (TARGET_SOURCES main.cpp m1_new.cpp gt.cpp orbital_am.cpp spin_am.cpp)
set (TEST_SOURCES test.cpp m1_new.cpp)

set (HEADERS constants.h rme_extras.h utility.h io_new.h chiral.h factory.h
             threedho.h quadrature.h m1.h gt.h orbital_am.h spin_am.h)

add_executable(${TARGET} ${TARGET_SOURCES} ${HEADERS})
add_executable(${TEST} ${TEST_SOURCES} ${HEADERS})

target_link_libraries(${TARGET} PUBLIC m basis am mcutils fmt::fmt CLI11::CLI11)
target_link_libraries(${TEST} PUBLIC m am fmt::fmt)

find_package(GSL REQUIRED)
target_link_libraries(${TARGET} PUBLIC GSL::gsl GSL::gslcblas)
target_link_libraries(${TEST} PUBLIC GSL::gsl GSL::gslcblas)

target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(${TEST} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

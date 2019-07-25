project (chime_src LANGUAGES CXX)

set (TARGET chime.exe)

set (SOURCES main.cpp chiral.cpp m1.cpp charge_radius.cpp)

set (HEADERS cli.h constants.h rme_extended.h utility.h io.h chiral.h factory.h
             operators.h threedho.h quadrature.h m1.h charge_radius.h)

add_executable(${TARGET} ${SOURCES} ${HEADERS})

target_link_libraries(${TARGET} PUBLIC m basis am mcutils fmt::fmt)

find_package(GSL REQUIRED)
target_link_libraries(${TARGET} PUBLIC GSL::gsl GSL::gslcblas)

target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

project (chime_src LANGUAGES CXX)

set (TARGET chime.exe)
set (TEST test.exe)

set (TARGET_SOURCES main.cpp m1_new.cpp gt.cpp)
set (TEST_SOURCES test.cpp)

set (HEADERS cli.h constants.h rme_extras.h utility.h io.h chiral.h factory.h
             operators.h threedho.h quadrature.h m1.h gt.h)

add_executable(${TARGET} ${TARGET_SOURCES} ${HEADERS})
add_executable(${TEST} ${TEST_SOURCES} ${HEADERS})

target_link_libraries(${TARGET} PUBLIC m basis am mcutils fmt::fmt)
target_link_libraries(${TEST} PUBLIC m am fmt::fmt)

find_package(GSL REQUIRED)
target_link_libraries(${TARGET} PUBLIC GSL::gsl GSL::gslcblas)

target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(${TEST} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

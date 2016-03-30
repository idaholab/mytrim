cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
project(mytrim)

set(DEBUG_FLAGS "-std=c++11 -g -O1 -Wall -Wextra -Werror -pedantic")
set(RELEASE_FLAGS "-std=c++11 -O3")

set(CMAKE_CXX_FLAGS ${RELEASE_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})
set(CMAKE_CONFIGURATION_TYPES Debug Release)

#
# build the mytrim library
# put install prefix somehow into config file
#
set( DATA_INSTALL /share/mytrim)
set( MYTRIM_DATA_DIR ${CMAKE_INSTALL_PREFIX}${DATA_INSTALL} )
configure_file( ${CMAKE_SOURCE_DIR}/config.h.cmake ${CMAKE_SOURCE_DIR}/config.h )
add_definitions(-DHAVE_CONFIG_H)
add_library(mytrim invert.C shim/cxx11random.C shim/point.C element.C material.C sample.C ion.C trim.C simconf.C sample_clusters.C sample_solid.C sample_wire.C sample_burried_wire.C sample_layers.C sample_dynamic.C)


install(TARGETS mytrim ARCHIVE DESTINATION lib)
install(FILES simconf.h shim/cxx11random.h shim/point.h element.h material.h ion.h trim.h invert.h sample.h functions.h DESTINATION include/mytrim)
install(FILES sample_clusters.h DESTINATION include/mytrim)
install(FILES sample_solid.h    DESTINATION include/mytrim)
install(FILES sample_wire.h     DESTINATION include/mytrim)
install(FILES sample_burried_wire.h     DESTINATION include/mytrim)
install(FILES sample_layers.h     DESTINATION include/mytrim)
install(FILES sample_dynamic.h     DESTINATION include/mytrim)
install(FILES data/SCOEF.95A data/SCOEF.95B data/SLFCTR.dat data/SNUC03.dat data/ELNAME.dat DESTINATION ${MYTRIM_DATA_DIR})

#
# build the test trim programs
#
include_directories(PUBLIC .)
set(PROGRAMS mytrim_ODS mytrim_clusters mytrim_layers mytrim_solid mytrim_solid2 mytrim_uo2 mytrim_wire mytrim_dynamic mytrim_wire2 mytrim_bobmsq mytrim_moose_verification)
foreach(program ${PROGRAMS})
  add_executable(${program} apps/${program}.C)
  target_link_libraries(${program}  mytrim m)
endforeach(program)

find_package(JsonCpp)

if (JSONCPP_FOUND)
  include_directories(${JSONCPP_INCLUDE_DIRS})
  add_executable(runmytrim apps/runmytrim.C
    apps/src/TrimVacCount.C
    apps/src/TrimVacEnergyCount.C
  )
  add_executable(runstopping apps/runstopping.C)
  target_link_libraries(runmytrim mytrim m ${JSONCPP_LIBRARIES})
  target_link_libraries(runstopping mytrim m ${JSONCPP_LIBRARIES})
endif (JSONCPP_FOUND)

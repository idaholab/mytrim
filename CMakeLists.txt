cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
project(mytrim C CXX)

option(ENABLE_GPERFTOOLS "enable the google performance tools (with GCC)" OFF)
option(ENABLE_LLVM_PROFILE "enable LLVM profiling instrumentation (with clang)" OFF)

find_package(Threads REQUIRED)

set(DEBUG_FLAGS "-std=c++11 -g -O1 -Wall -Wextra -Werror -pedantic")
set(RELEASE_FLAGS "-std=c++11")

if(ENABLE_GPERFTOOLS)
  # find_package(Gperftools REQUIRED)
  set(CMAKE_CXX_FLAGS "${RELEASE_FLAGS} -g -O2")
  set(CMAKE_EXE_LINKER_FLAGS "-g")
  set(CMAKE_SHARED_LINKER_FLAGS "-g")
  set(OPTIONAL_LIBRARIES ${OPTIONAL_LIBRARIES} profiler)
elseif(ENABLE_LLVM_PROFILE)
  set(CMAKE_CXX_FLAGS "${RELEASE_FLAGS} -fprofile-instr-generate -O2")
  set(CMAKE_EXE_LINKER_FLAGS "-fprofile-instr-generate")
  set(CMAKE_SHARED_LINKER_FLAGS "-fprofile-instr-generate")
else()
  set(CMAKE_CXX_FLAGS "${RELEASE_FLAGS} -O3")
endif()
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
add_library(mytrim invert.C shim/point.C element.C material.C sample.C ion.C trim.C simconf.C sample_clusters.C sample_solid.C sample_wire.C sample_burried_wire.C sample_layers.C sample_dynamic.C)


install(TARGETS mytrim ARCHIVE DESTINATION lib)
install(FILES simconf.h shim/point.h element.h material.h ion.h trim.h invert.h sample.h functions.h DESTINATION include/mytrim)
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
set(PROGRAMS mytrim_ODS mytrim_clusters mytrim_layers mytrim_solid mytrim_solid2 mytrim_uo2 mytrim_wire mytrim_dynamic mytrim_wire2 mytrim_bobmsq mytrim_moose_verification mytrim_distance)
foreach(program ${PROGRAMS})
  add_executable(${program} apps/${program}.C)
  target_link_libraries(${program}  mytrim m ${OPTIONAL_LIBRARIES})
endforeach(program)

find_package(JsonCpp)

if (JSONCPP_FOUND)
message(" include json: ${JSONCPP_INCLUDE_DIRS}")
message(" lib json: ${JSONCPP_LIBRARIES}")
message(" prefix json: ${JSONCPP_INCLUDE_PREFIX}") 

  include_directories(${JSONCPP_INCLUDE_DIRS})
  add_executable(runmytrim apps/runmytrim.C
    apps/src/TrimRange.C
    apps/src/TrimVacCount.C
    apps/src/TrimVacEnergyCount.C
  )
  add_executable(runstopping apps/runstopping.C)
  target_link_libraries(runmytrim mytrim m ${JSONCPP_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${OPTIONAL_LIBRARIES})
  target_link_libraries(runstopping mytrim m ${JSONCPP_LIBRARIES} ${OPTIONAL_LIBRARIES})
endif (JSONCPP_FOUND)

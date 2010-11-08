cmake_minimum_required(VERSION 2.6)
project(mytrim)

find_library(R250C_LIB NAMES r250c PATHS /opt/lib)

# build the mytrim library
# put install prefix somehow into config file
set( DATA_INSTALL /share/mytrim)
set( DATA_DIR ${CMAKE_INSTALL_PREFIX}${DATA_INSTALL} )
configure_file( ${CMAKE_SOURCE_DIR}/config.h.cmake ${CMAKE_SOURCE_DIR}/config.h )
add_library(mytrim invert.cc element.cc material.cc sample.cc ion.cc trim.cc simconf.cc sample_clusters.cc sample_solid.cc sample_wire.cc sample_burried_wire.cc sample_layers.cc sample_dynamic.cc)
install(TARGETS mytrim ARCHIVE DESTINATION lib)
install(FILES simconf.h element.h material.h ion.h trim.h invert.h sample.h DESTINATION include/mytrim)
install(FILES sample_clusters.h DESTINATION include/mytrim)
install(FILES sample_solid.h    DESTINATION include/mytrim)
install(FILES sample_wire.h     DESTINATION include/mytrim)
install(FILES sample_burried_wire.h     DESTINATION include/mytrim)
install(FILES sample_layers.h     DESTINATION include/mytrim)
install(FILES sample_dynamic.h     DESTINATION include/mytrim)
install(FILES data/SCOEF.95A data/SLFCTR.dat data/SNUC03.dat DESTINATION ${DATA_DIR})

# build the test trim programs
set(PROGRAMS mytrim_ODS mytrim_clusters mytrim_layers mytrim_solid mytrim_solid2 mytrim_uo2 mytrim_wire mytrim_dynamic mytrim_wire2)
foreach(program ${PROGRAMS})
  add_executable(${program} ${program}.cc)
  target_link_libraries(${program} ${R250C_LIB} mytrim m)
endforeach(program)
  

# might use boost sometime
#find_package(Boost REQUIRED program_options)
#include_directories(${Boost_INCLUDE_DIRS})
#target_link_libraries(firstpass ${Boost_LIBRARIES} ${R250C_LIB} m)


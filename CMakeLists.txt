cmake_minimum_required(VERSION 2.6)
project(mytrim)

# build the mytrim library
# put install prefix somehow into config file
set( DATA_INSTALL /share/mytrim)
set( DATA_DIR ${CMAKE_INSTALL_PREFIX}${DATA_INSTALL} )
configure_file( ${CMAKE_SOURCE_DIR}/config.h.cmake ${CMAKE_SOURCE_DIR}/config.h )
add_library(mytrim invert.C element.C material.C sample.C ion.C trim.C simconf.C sample_clusters.C sample_solid.C sample_wire.C sample_burried_wire.C sample_layers.C sample_dynamic.C)
install(TARGETS mytrim ARCHIVE DESTINATION lib)
install(FILES simconf.h element.h material.h ion.h trim.h invert.h sample.h functions.h DESTINATION include/mytrim)
install(FILES sample_clusters.h DESTINATION include/mytrim)
install(FILES sample_solid.h    DESTINATION include/mytrim)
install(FILES sample_wire.h     DESTINATION include/mytrim)
install(FILES sample_burried_wire.h     DESTINATION include/mytrim)
install(FILES sample_layers.h     DESTINATION include/mytrim)
install(FILES sample_dynamic.h     DESTINATION include/mytrim)
install(FILES data/SCOEF.95A data/SLFCTR.dat data/SNUC03.dat data/ELNAME.dat DESTINATION ${DATA_DIR})

# build the test trim programs
set(PROGRAMS mytrim_ODS mytrim_clusters mytrim_layers mytrim_solid mytrim_solid2 mytrim_uo2 mytrim_wire mytrim_dynamic mytrim_wire2 mytrim_bobmsq)
foreach(program ${PROGRAMS})
  add_executable(${program} ${program}.C)
  target_link_libraries(${program}  mytrim m)
endforeach(program)


# might use boost sometime
#find_package(Boost REQUIRED program_options)
#include_directories(${Boost_INCLUDE_DIRS})
#target_link_libraries(firstpass ${Boost_LIBRARIES} ${R250C_LIB} m)

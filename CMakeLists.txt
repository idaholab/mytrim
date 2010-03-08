cmake_minimum_required(VERSION 2.6)
project(mytrim)

find_library(R250C_LIB NAMES r250c PATHS /opt/lib)

# build the mytrim library
add_library(mytrim invert.cc element.cc material.cc sample.cc ion.cc trim.cc simconf.cc sample_clusters.cc sample_solid.cc sample_wire.cc)
install(TARGETS mytrim ARCHIVE DESTINATION lib)
install(FILES simconf.h element.h material.h ion.h trim.h invert.h sample.h sample_solid.h sample_wire.h DESTINATION include/mytrim)
install(FILES sample_clusters.h DESTINATION include/mytrim)
install(FILES sample_solid.h    DESTINATION include/mytrim)
install(FILES sample_wire.h     DESTINATION include/mytrim)

# build the test trim programs
set(PROGRAMS mytrim_clusters mytrim_solid mytrim_uo2 mytrim_wire)
foreach(program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(firstpass ${R250C_LIB} mytrim m)
endforeach(program)
  

# might use boost sometime
#find_package(Boost REQUIRED program_options)
#include_directories(${Boost_INCLUDE_DIRS})
#target_link_libraries(firstpass ${Boost_LIBRARIES} ${R250C_LIB} m)


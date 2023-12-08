cmake_minimum_required(VERSION 3.15)
project(clk VERSION 0.0.1 LANGUAGES C)

set(FLAGS_COMMON "-Wall -Wextra")
set(CMAKE_C_FLAGS "${FLAGS_COMMON} -O3 -g0")
set(CMAKE_C_FLAGS_DEBUG "${FLAGS_COMMON} -O0 -g2 -ggdb")

set(source_files
    edgegen.c  xnear.c delaunay.c mlinkern.c
    fmatch.c
    kdbuild.c kdnear.c kdspan.c kdtwoopt.c
    allocrus.c bgetopt.c  dheaps_i.c edgelen.c  edgemap.c  edgeutil.c
    fastread.c getdata.c  safe_io.c  sortrus.c
    urandom.c  util.c     zeit.c
    linkern.c linkern_fixed.c linkern_path.c flip_two.c
)

set(library_name clk)
set(target_name lk_main)

link_libraries(m)

add_library(${library_name} STATIC ${source_files})

add_executable(${target_name} ${target_name}.c)

target_link_libraries(
    ${target_name} PRIVATE
    ${library_name}
)

add_custom_target(clk-test
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_name} pr2392.tsp
)

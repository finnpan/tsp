cmake_minimum_required(VERSION 3.15)
project(clk VERSION 0.0.1 LANGUAGES C)

set(FLAGS_COMMON "-Wall -Werror -Wextra -Wno-unused-function -Wno-unused-but-set-variable -Wno-maybe-uninitialized")
set(CMAKE_C_FLAGS "${FLAGS_COMMON} -O3 -g0")
set(CMAKE_C_FLAGS_DEBUG "${FLAGS_COMMON} -O0 -g2 -ggdb")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
link_libraries(m)

add_subdirectory(edgegen)
add_subdirectory(fmatch)
add_subdirectory(kdtree)
add_subdirectory(linkern)
add_subdirectory(util)

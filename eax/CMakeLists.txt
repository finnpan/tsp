cmake_minimum_required(VERSION 3.15)
project(eax VERSION 0.0.1 LANGUAGES CXX)

set(FLAGS_COMMON "-Wall -Werror -Wextra -Wno-maybe-uninitialized")
set(CMAKE_CXX_FLAGS "${FLAGS_COMMON} -O3 -g0")
set(CMAKE_CXX_FLAGS_DEBUG "${FLAGS_COMMON} -O0 -g2 -ggdb")

set(library_name eax_framework)
set(source_path ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(test)


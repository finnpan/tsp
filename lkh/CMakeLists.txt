cmake_minimum_required(VERSION 3.15)
project(lkh VERSION 0.0.1 LANGUAGES C)

set(FLAGS_COMMON "-fcommon -Wall -Werror -Wextra -Wno-unused-parameter")
set(CMAKE_C_FLAGS "${FLAGS_COMMON} -O3 -g0")
set(CMAKE_C_FLAGS_DEBUG "${FLAGS_COMMON} -O0 -g2 -ggdb")

set(library_name lkh)
include_directories(include)
link_libraries(m)

add_subdirectory(src)
add_subdirectory(test_lkh)

cmake_minimum_required(VERSION 3.15)
project(clk)

set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -Wall -Werror -Wextra -O0 -g2 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -Wall -Werror -Wextra -O3")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
link_libraries(m)

add_subdirectory(edgegen)
add_subdirectory(fmatch)
add_subdirectory(kdtree)
add_subdirectory(linkern)
add_subdirectory(util)

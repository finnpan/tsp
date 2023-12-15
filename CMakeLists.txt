cmake_minimum_required(VERSION 3.15)
project(tsp)

set(base_lib "base")
add_subdirectory(${base_lib})

include_directories(${base_lib})
link_directories(./${base_lib})

add_subdirectory(eax)
add_subdirectory(clk)

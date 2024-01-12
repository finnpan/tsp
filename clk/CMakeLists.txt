cmake_minimum_required(VERSION 3.15)
project(clk VERSION 0.0.1 LANGUAGES CXX)

set(FLAGS_COMMON "-Wall -Werror -Wextra")
set(CMAKE_CXX_FLAGS "${FLAGS_COMMON} -O3 -g0")
set(CMAKE_CXX_FLAGS_DEBUG "${FLAGS_COMMON} -O0 -g2 -ggdb")


set(target_name clk)
add_executable(${target_name} ${target_name}.cpp util.cpp linkern.cpp)
target_link_libraries(${target_name} PRIVATE m ${base_lib})


add_custom_target(clk-test
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
    COMMAND ./${target_name} ${CMAKE_SOURCE_DIR}/instances/pr2392.tsp
)

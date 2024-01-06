cmake_minimum_required(VERSION 3.15)
project(clk VERSION 0.0.1 LANGUAGES C)

set(FLAGS_COMMON "-Wall -Werror -Wextra")
set(CMAKE_C_FLAGS "${FLAGS_COMMON} -O3 -g0")
set(CMAKE_C_FLAGS_DEBUG "${FLAGS_COMMON} -O0 -g2 -ggdb")


set(target_name clk)
add_executable(${target_name} ${target_name}.c util.c linkern.c)
target_link_libraries(${target_name} PRIVATE m ${base_lib})


add_custom_target(clk-test
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_name} -r 10 pr2392.tsp
)

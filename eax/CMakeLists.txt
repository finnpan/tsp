cmake_minimum_required(VERSION 3.15)
project(eax VERSION 0.0.1 LANGUAGES CXX)

set(FLAGS_COMMON "-Wall -Werror -Wextra")
set(CMAKE_CXX_FLAGS "${FLAGS_COMMON} -O3 -g0")
set(CMAKE_CXX_FLAGS_DEBUG "${FLAGS_COMMON} -O0 -g2 -ggdb -fsanitize=address")

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    link_libraries(asan)
endif()


link_libraries(${base_lib})


set(library_name eax)
add_library(${library_name} STATIC eax.cpp util.cpp)
link_libraries(${library_name})


set(target_gatest gatest)
add_executable(${target_gatest} ${target_gatest}.cpp)
add_custom_target(ga-test
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
)


set(target_kopt kopttest)
add_executable(${target_kopt} ${target_kopt}.cpp)
add_custom_target(kopt-test
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_kopt} fnl4461.tsp
)


add_custom_target(eax-test
    COMMAND make kopt-test
    COMMAND make ga-test
)


add_custom_target(heavy-eax-test
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} fnl4461.tsp
)


add_custom_target(heap-overflow-eax-test
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_gatest} eil101.tsp
)

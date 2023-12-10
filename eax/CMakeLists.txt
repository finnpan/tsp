cmake_minimum_required(VERSION 3.15)
project(eax VERSION 0.0.1 LANGUAGES CXX)

set(FLAGS_COMMON "-Wall -Werror -Wextra")
set(CMAKE_CXX_FLAGS "${FLAGS_COMMON} -O3 -g0")
set(CMAKE_CXX_FLAGS_DEBUG "${FLAGS_COMMON} -O0 -g2 -ggdb -fsanitize=address ")

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    link_libraries(asan)
endif()


set(library_name eax)
add_library(${library_name} STATIC eax.cpp tree.cpp tsp.cpp)
link_libraries(${library_name})


set(target_eaxtest eaxtest)
add_executable(${target_eaxtest} ${target_eaxtest}.cpp)
add_custom_target(eax-test
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_eaxtest} eil101.tsp
)


set(target_tree treetest)
add_executable(${target_tree} ${target_tree}.cpp)
target_link_libraries(${target_tree} PRIVATE gtest gtest_main pthread)
add_custom_target(tree-test COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_tree})


set(target_kopt kopttest)
add_executable(${target_kopt} ${target_kopt}.cpp)
add_custom_target(kopt-test
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_kopt} fnl4461.tsp
)


add_custom_target(test
    COMMAND make tree-test
    COMMAND make kopt-test
    COMMAND make eax-test
)


add_custom_target(heap-overflow-test
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/heap-overflow-test.sh
)


add_custom_target(heavy-eax-test
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_eaxtest} fnl4461.tsp
)

include_directories(${source_path})
link_libraries(${library_name})

set(target_name eaxtest)
add_executable(${target_name} ${target_name}.cpp)
add_custom_target(test-eax
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/../instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_name} eil101.tsp
)

set(target_name treetest)
add_executable(${target_name} ${target_name}.cpp)
target_link_libraries(${target_name} PRIVATE gtest gtest_main pthread)
add_custom_target(test-tree COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_name})

set(target_name kopttest)
add_executable(${target_name} ${target_name}.cpp)
add_custom_target(test-kopt
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/../instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${target_name} fnl4461.tsp
)

add_custom_target(test
    COMMAND make test-tree
    COMMAND make test-kopt
    COMMAND make test-eax
)

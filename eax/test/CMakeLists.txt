set(target_name eaxtest)
add_executable(${target_name} ${target_name}.cpp)
target_include_directories(${target_name} PRIVATE ${source_path})
target_link_libraries(${target_name} PRIVATE ${library_name})


set(target_name twoleveltreetest)
add_executable(${target_name} ${target_name}.cpp)
target_include_directories(${target_name} PRIVATE ${source_path})
target_link_libraries(${target_name} PRIVATE gtest gtest_main pthread ${library_name})

cmake_minimum_required(VERSION 3.15)
project(lkh VERSION 0.0.1 LANGUAGES C)

set(FLAGS_COMMON "-fcommon -Wall -Werror -Wextra")
set(CMAKE_C_FLAGS "${FLAGS_COMMON} -O3 -g0")
set(CMAKE_C_FLAGS_DEBUG "${FLAGS_COMMON} -O0 -g2 -ggdb")

include_directories(include)
link_libraries(m)

aux_source_directory(src source_files)

add_executable(${PROJECT_NAME} ${source_files})

add_custom_target(
    test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../instances/
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME} pr2392.par
)

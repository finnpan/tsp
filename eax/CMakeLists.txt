cmake_minimum_required(VERSION 3.15)
project(eax)

set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -Wall -Werror -Wextra -O0 -g2 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -Wall -Werror -Wextra -O3")

set(library_name eax_framework)
set(source_path ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(test_eax)
add_subdirectory(test_twoleveltree)

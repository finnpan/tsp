cmake_minimum_required(VERSION 3.15)
project(lkh)

set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -Wall -Werror -Wextra -O0 -g2 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -Wall -Werror -Wextra -O3")

set(library_name lkh)
include_directories(include)
link_libraries(m)

add_subdirectory(src)
add_subdirectory(test_lkh)

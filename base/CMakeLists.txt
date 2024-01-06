cmake_minimum_required(VERSION 3.15)
project(base VERSION 0.0.1 LANGUAGES C)

set(FLAGS_COMMON "-Wall -Werror -Wextra")
set(CMAKE_C_FLAGS "${FLAGS_COMMON} -O3 -g0")
set(CMAKE_C_FLAGS_DEBUG "${FLAGS_COMMON} -O0 -g2 -ggdb")

add_library(${base_lib} STATIC alloc.c flipper.c mpool.c)

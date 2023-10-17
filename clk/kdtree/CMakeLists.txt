set(
    source_files
    kdbuild.c kdnear.c kdspan.c kdtwoopt.c
)

set(library_name kdtree)
set(target_name kd_main)

add_library(${library_name} STATIC ${source_files})

add_executable(${target_name} ${target_name}.c)

target_link_libraries(
    ${target_name} PRIVATE
    ${library_name} util
)

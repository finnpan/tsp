set(source_files fmatch.c)

set(library_name fmatch)
set(target_name fm_main)

add_library(${library_name} STATIC ${source_files})

add_executable(${target_name} ${target_name}.c)

target_link_libraries(
    ${target_name} PRIVATE
    ${library_name} edgegen kdtree linkern util
)

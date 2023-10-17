set(
    source_files
    edgegen.c  xnear.c delaunay.c mlinkern.c
)

set(library_name edgegen)
set(target_name eg_main)

add_library(${library_name} STATIC ${source_files})

add_executable(${target_name} ${target_name}.c)

target_link_libraries(
    ${target_name} PRIVATE
    ${library_name} kdtree fmatch linkern util
)

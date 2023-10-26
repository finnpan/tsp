set(source_files
    linkern.c linkern_fixed.c linkern_path.c flip_two.c
)

set(library_name linkern)
set(target_name lk_main)

add_library(${library_name} STATIC ${source_files})

add_executable(${target_name} ${target_name}.c)

target_link_libraries(
    ${target_name} PRIVATE
    ${library_name} edgegen kdtree fmatch util
)

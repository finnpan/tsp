set(source_files
    eax.cpp
    indi.cpp
    kopt.cpp
    twoleveltree.cpp
)

add_library(${library_name} STATIC ${source_files})

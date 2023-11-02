set(target_name test_twoleveltree)

find_package(GTest REQUIRED)

add_executable(${target_name} ${target_name}.cpp)

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(gtest_repo_path "C:/Users/HX00711/Desktop/repo/googletest-1.14.0")
    message("On Windows! gtest_repo_path is " ${gtest_repo_path})
    target_link_directories(${target_name} PRIVATE
        "${gtest_repo_path}/build/lib"
    )
    list(APPEND source_path "${gtest_repo_path}/googletest/include")
endif()

target_include_directories(${target_name} PRIVATE ${source_path})

target_link_libraries(
    ${target_name} PRIVATE
    gtest gtest_main pthread ${library_name}
)

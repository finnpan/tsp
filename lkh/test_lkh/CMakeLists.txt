set(target_name LKHmain)

add_executable(${target_name} ${target_name}.c)

target_link_libraries(
    ${target_name} PRIVATE
    ${library_name}
)

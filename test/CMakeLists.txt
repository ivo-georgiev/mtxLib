cmake_minimum_required(VERSION 3.3)
project(matrix_test)

# Compiler configuration
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Werror -std=c99 -pedantic)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
    add_compile_options(/W4 /WX)
endif()

# Header include paths
include_directories(..)

# Executable targets
add_executable(matrix_test mtxTest.c ../mtxLib.c)
target_link_libraries(matrix_test m)

enable_testing()
add_test(testname matrix_test)



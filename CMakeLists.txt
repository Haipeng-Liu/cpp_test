cmake_minimum_required(VERSION 3.16)
project(test)

set(CMAKE_CXX_STANDARD 17)
find_package(Thread)
add_executable(test
        main.cpp math_test.h thread_test.h)
target_link_libraries(test pthread  ${CMAKE_THREAD_LIBS_INIT})
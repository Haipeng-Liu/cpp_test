cmake_minimum_required(VERSION 3.16)
project(test)

set(CMAKE_CXX_STANDARD 20)
find_package(Thread)

file(GLOB_RECURSE COMP_SRC ./*.cpp ./*.c ./*.cc)
file(GLOB_RECURSE COMP_HDR ./*.h ./*.hpp )
add_executable(test
#        ${COMP_SRC} ${COMP_HDR}
        main.cpp math_test.h thread_test.h math_common/data_type.h kalman.h
        math_test.h math_common/iou.h math_common/iou.cpp math_common/common.h math_common/common.cpp
        )
target_link_libraries(test pthread  ${CMAKE_THREAD_LIBS_INIT})
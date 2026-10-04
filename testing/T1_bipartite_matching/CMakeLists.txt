cmake_minimum_required(VERSION 3.22)
project(T1_bipartite_matching_grader LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 正式测评必须由总控脚本传入考生提交仓库的绝对路径；本地默认用当前仓库。
if(NOT DEFINED SUBMISSION_ROOT OR SUBMISSION_ROOT STREQUAL "")
    set(SUBMISSION_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
endif()

get_filename_component(SUBMISSION_ROOT "${SUBMISSION_ROOT}" ABSOLUTE)
set(STUDENT_SOURCE_ROOT
    "${SUBMISSION_ROOT}/problems/T1_bipartite_matching/src")
set(STUDENT_INCLUDE_ROOT
    "${SUBMISSION_ROOT}/problems/T1_bipartite_matching/include")

if(NOT IS_DIRECTORY "${STUDENT_SOURCE_ROOT}")
    message(FATAL_ERROR "找不到考生源码目录: ${STUDENT_SOURCE_ROOT}")
endif()

# 与正式测评一致：递归收集本题 src/ 下全部 .cpp，使用官方接口头文件编译。
file(GLOB_RECURSE STUDENT_SOURCES CONFIGURE_DEPENDS
    "${STUDENT_SOURCE_ROOT}/*.cpp")
list(SORT STUDENT_SOURCES)

if(NOT STUDENT_SOURCES)
    message(FATAL_ERROR "考生 src/ 目录中没有找到 .cpp 文件")
endif()

add_library(candidate OBJECT ${STUDENT_SOURCES})
target_compile_features(candidate PRIVATE cxx_std_17)
target_compile_options(candidate PRIVATE -O2 -Wall -Wextra)
target_include_directories(candidate PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/include"
    "${STUDENT_INCLUDE_ROOT}")

add_executable(grader
    grader.cpp
    $<TARGET_OBJECTS:candidate>)
target_compile_features(grader PRIVATE cxx_std_17)
target_compile_options(grader PRIVATE -O2 -Wall -Wextra)
target_include_directories(grader PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/include")

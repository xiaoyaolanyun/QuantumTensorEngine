# CMakeList.txt: TensorEngine 的 CMake 项目，包括源和定义
# 此处特定于项目的逻辑。
#
cmake_minimum_required (VERSION 3.8)

file(GLOB_RECURSE TensorEngine_SRCS_H ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB_RECURSE TensorEngine_SRCS_CPP ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp )
# 将源添加到此项目的可执行文件。
add_executable (TensorEngine ${TensorEngine_SRCS_H} ${TensorEngine_SRCS_CPP})

# TODO: 如有需要，请添加测试并安装目标。

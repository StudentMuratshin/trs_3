# CMakeList.txt : CMake project for trs_3, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

include_directories("C:\\Users\\PETA4\\Desktop\\3\\6_sem\\trs\\labs\\trs_3\\trs_3\\eigen-3.4.0")

# Add source to this project's executable.
add_executable (task_1 "task_1.cpp")
add_executable (task_2 "task_2.cpp")
add_executable (task_3 "task_3.cpp")
add_executable (task_4 "task_4.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET task_1 PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.

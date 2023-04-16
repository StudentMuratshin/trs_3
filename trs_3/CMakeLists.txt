# CMakeList.txt : CMake project for trs_3, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (task_1 "task_1.cpp")
add_executable (task_2 "task_2.cpp")
add_executable (task_3 "task_3.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET task_1 PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.

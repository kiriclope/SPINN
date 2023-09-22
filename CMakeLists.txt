cmake_minimum_required(VERSION 3.1)
project(LifNet)

# Ensure we have C++11 or later
set (CMAKE_CXX_STANDARD 17)

# Try to find the yaml-cpp package
find_package(yaml-cpp REQUIRED)

add_executable(LifNet src/main.cpp)

# Link yaml-cpp to the executable
target_link_libraries(LifNet ${YAML_CPP_LIBRARIES})

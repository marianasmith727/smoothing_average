cmake_minimum_required(VERSION 2.8.3)
project(smoothing_average)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  geometry_msgs
)

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
add_compile_options(-std=c++14)

set(EIGEN3_INCLUDE_DIR /usr/include/eigen3)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS roscpp rospy std_msgs geometry_msgs
#  DEPENDS system_lib
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)


add_executable(smoothing_average src/smoothing_average.cpp)
target_link_libraries(smoothing_average ${catkin_LIBRARIES})


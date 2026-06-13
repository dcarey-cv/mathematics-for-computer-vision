cmake_minimum_required(VERSION 3.20)
project(mathematics-for-computer-vision)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(fmt REQUIRED)
find_package(Eigen3 REQUIRED)

include(FetchContent)
FetchContent_Declare(
  pangolin
  GIT_REPOSITORY https://github.com/stevenlovegrove/Pangolin.git
  GIT_TAG v0.9.1
)
FetchContent_MakeAvailable(pangolin)

option(STRICT_WARNINGS "Enable -Wall -Wextra -Werror -Wpedantic" ON)

add_subdirectory(strang-lingear-algebra)

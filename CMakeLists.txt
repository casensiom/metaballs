cmake_minimum_required(VERSION 3.16)

project(metaball)
set(CMAKE_C_STANDARD 99)
# set(CMAKE_C_STANDARD_REQUIRED TRUE)

# raylib
find_package(raylib QUIET)
if (NOT raylib_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        raylib
        GIT_REPOSITORY https://github.com/raysan5/raylib.git
        GIT_TAG 4.5.0
    )
    FetchContent_MakeAvailable(raylib)
endif()

message(STATUS "PROJECT_NAME: ${PROJECT_NAME}")
add_executable(${PROJECT_NAME} metaball.c marching_cube.h)
target_link_libraries(${PROJECT_NAME} PUBLIC raylib)
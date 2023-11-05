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
add_executable(${PROJECT_NAME} main.c metaballs.h)
target_link_libraries(${PROJECT_NAME} PUBLIC raylib)

add_custom_target(copy_assets
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/assets ${CMAKE_CURRENT_BINARY_DIR}/assets
)
add_dependencies(${PROJECT_NAME} copy_assets)
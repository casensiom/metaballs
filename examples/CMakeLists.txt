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

add_library(metaball_dependency 
    raylib_wrapper.h 
    raylib_wrapper.c 
    ball.h 
    ball.c
)
target_link_libraries(metaball_dependency PUBLIC raylib)

add_custom_target(copy_assets
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/assets ${CMAKE_CURRENT_BINARY_DIR}/assets
)
add_dependencies(metaball_dependency copy_assets)

add_executable(example0 example0.c)
target_link_libraries(example0 PUBLIC metaball_dependency)

add_executable(example1 example1.c)
target_link_libraries(example1 PUBLIC metaball_dependency)

add_executable(example2 example2.c)
target_link_libraries(example2 PUBLIC metaball_dependency)


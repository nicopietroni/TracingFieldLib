cmake_minimum_required(VERSION 3.5)
project(quad_gui)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# Build type: change this word and run make (cmake reconfigures by itself)
#   Release : fast, asserts off (default)
#   Debug   : asserts on, -g -O0 (for lldb)
#   Sanity  : Debug + AddressSanitizer/UBSan: on a crash it prints file and
#             line of the bad access (much slower)
#   RelSanity : Release (same -O3 numbers, asserts off) + AddressSanitizer
#             and -g: for a crash of the Release only (Debug computes
#             slightly different numbers and may take another path)
set(MY_BUILD_TYPE "Release")

if(MY_BUILD_TYPE STREQUAL "Sanity")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined)
elseif(MY_BUILD_TYPE STREQUAL "RelSanity")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    add_compile_options(-g -fsanitize=address -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address)
elseif(MY_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
else()
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "MY_BUILD_TYPE = ${MY_BUILD_TYPE}, CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# as the command line quadrangulator: the LUT of the patch quadrangulation
# (relative to the build folder, as there)
set(LUT_PATH "../../../include/Lut.txt")
add_compile_definitions(LUT_PATH="${LUT_PATH}")

set(IMGUI_DIR "/Users/134750/Desktop/imgui/" CACHE INTERNAL "")
set(TRACEFIELD_DIR = "../../include" CACHE INTERNAL "")
set(GLDRAW_DIR = "../include/GLFunctions/" CACHE INTERNAL "")
set(EIGENPATH  = "../../external_lib/eigen-3.4.0/" CACHE INTERNAL "")

set(IGLLIB  = "../../external_lib/libigl/include/" CACHE INTERNAL "")

include_directories(${TRACEFIELD_DIR})
include_directories(${IMGUI_DIR})
include_directories(${IMGUI_DIR}/backends)
include_directories(/opt/homebrew/include)
include_directories(${EIGENPATH})

include_directories(${IGLLIB})

find_library(GLEW_LIB NAMES GLEW glew PATHS /opt/homebrew/lib /usr/local/lib)
find_path(GLEW_INCLUDE_DIR GL/glew.h PATHS /opt/homebrew/include /usr/local/include)
if(NOT GLEW_LIB OR NOT GLEW_INCLUDE_DIR)
  message(FATAL_ERROR "GLEW not found")
endif()

add_executable(${PROJECT_NAME} main.cpp
               ${IMGUI_DIR}/imgui.cpp
               ${IMGUI_DIR}/imgui_demo.cpp
               ${IMGUI_DIR}/imgui_draw.cpp
               ${IMGUI_DIR}/imgui_tables.cpp
               ${IMGUI_DIR}/imgui_widgets.cpp
               ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
               ${IMGUI_DIR}/backends/imgui_impl_opengl2.cpp
               )

target_include_directories(${PROJECT_NAME} PRIVATE ${GLEW_INCLUDE_DIR})

# as the command line quadrangulator (the flow graph uses threads)
find_package(Threads REQUIRED)
target_compile_options(${PROJECT_NAME} PRIVATE
  $<$<CONFIG:Release>:-O3 -DNDEBUG>)

target_link_libraries(${PROJECT_NAME}
    ${GLEW_LIB}
    Threads::Threads
    "-framework OpenGL"
    "-framework Cocoa"
    "-framework IOKit"
    "-framework CoreVideo"
    "-L/opt/homebrew/lib"
    "-lglfw"
)
